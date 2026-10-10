import os
import tempfile
import unittest
from unittest.mock import patch

import mfem.ser as mfem
import numpy as np

from bfp.assemble import BoundaryConditions, FESpace, Bilinear_Form, Linear_Form
from bfp.coeff import VectorConstCoefficient
from bfp.input import problem1_input, problem2_input, problem3_input, problem4_input, problem5_input
from bfp.mesh import create_2D_mesh
from bfp.solver import Solve_Psi, Solve_Phi, GMRES_solver, Initial_Guess


class AnalyticCoefficient(mfem.PyCoefficient):
    def __init__(self, function):
        super().__init__()
        self.function = function

    def EvalValue(self, x):
        return self.function(x)


class TestProblems(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        previous = os.getcwd()
        os.chdir(directory.name)
        self.addCleanup(os.chdir, previous)
        display = patch('bfp.vis.GlVis_2D')
        display.start()
        self.addCleanup(display.stop)

    def test_optional_output(self):
        mesh = mfem.Mesh(2, 3, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        mu, weights = np.polynomial.legendre.leggauss(2)
        xs = mfem.ConstantCoefficient(5.0)
        args = (5, mu, weights, mesh, fes, xs, 5.0, 5.0, 2.0,
                1.0, 0.5, *BoundaryConditions(mesh), 5.0, 10.0)

        with patch('bfp.vis.GlVis_2D') as display:
            angular = Solve_Psi(*args, p_level=0, visualize=False, save=False)
            phi = Solve_Phi(fes, angular, save=False)
            display.assert_not_called()
            self.assertEqual(os.listdir(), [])

            with_output = Solve_Psi(*args, p_level=0)
            saved_phi = Solve_Phi(fes, with_output)
            self.assertEqual(display.call_count, len(mu))
        self.assertEqual(set(os.listdir()), {
            'phi_dg.gf', *('psi_hs_mu_{:.3f}.gf'.format(m) for m in mu)})
        np.testing.assert_array_equal(phi.GetDataArray(), saved_phi.GetDataArray())
        for (_, _, psi), (_, _, saved) in zip(angular, with_output):
            np.testing.assert_array_equal(psi.GetDataArray(), saved.GetDataArray())

    def test_manufactured_problems(self):
        cases = [
            (problem1_input, lambda x: 40.0),
            (problem3_input, lambda x: 2.0 * (5.0 + 10.0 * x[0])),
            (problem4_input, lambda x: 2.0 * (5.0 + 10.0 * x[1])),
            (problem5_input, lambda x: 2.0 * (5.0 + 10.0 * x[0] * x[1])),
        ]
        for problem, exact in cases:
            with self.subTest(problem=problem.__name__):
                phi = problem(nx=4, nE=4, N_ang=2, tol=1e-12, p_level=0)
                self.assertLess(phi.ComputeL2Error(AnalyticCoefficient(exact)), 1e-8)

    def solve_energy(self, nE, stopping, sigma, inflow):
        mesh = create_2D_mesh(2, nE, 0.0, 1.0, 0.25, 2.25)
        fes = FESpace(1, mesh)
        markers = BoundaryConditions(mesh)
        xs = mfem.ConstantCoefficient(sigma)
        psi = Solve_Psi(1, [0.0], [1.0], mesh, fes, xs, sigma, inflow, 1.0,
                        1.0, 0.5, *markers, 0.0, 0.0, p_level=0,
                        S_coeff=stopping)[0][2]
        return mesh, fes, markers, xs, psi

    def test_zero_stopping_power(self):
        exact = AnalyticCoefficient(lambda x: 4.0 + 8.0 * (
            np.exp(-5.0 * np.sqrt(3.0) * x[0]) +
            np.exp(-5.0 * np.sqrt(3.0) * (1.0 - x[0]))))
        errors = []
        for nx in (16, 32):
            phi = problem2_input(nx=nx, nE=2, p_level=0)
            errors.append(phi.ComputeL2Error(exact))
        self.assertLess(errors[1], 0.02)
        self.assertGreater(errors[0] / errors[1], 3.5)

    def test_high_energy_inflow(self):
        # The deliberately different low-E value must not enter the solution.
        inflow = AnalyticCoefficient(lambda x: 1.0 if x[1] > 2.0 else 17.0)
        exact = AnalyticCoefficient(lambda x: np.exp(x[1] - 2.25))
        errors = []
        for nE in (8, 16):
            mesh, fes, markers, xs, psi = self.solve_energy(
                nE, mfem.ConstantCoefficient(1.0), 1.0, inflow)
            errors.append(psi.ComputeL2Error(exact))
        self.assertLess(errors[1], 1e-3)
        self.assertGreater(errors[0] / errors[1], 3.5)

    def test_variable_stopping_power(self):
        # -d[(1+E)*psi]/dE = 0, with unit incoming energy-space flux.
        stopping = AnalyticCoefficient(lambda x: 1.0 + x[1])
        exact = AnalyticCoefficient(lambda x: 1.0 / (1.0 + x[1]))
        errors = []
        for nE in (8, 16):
            mesh, fes, markers, xs, psi = self.solve_energy(nE, stopping, 0.0, exact)
            errors.append(psi.ComputeL2Error(exact))
            _, low_energy = BoundaryConditions(mesh, y_min=1, y_max=0)
            outgoing = mfem.LinearForm(fes)
            outgoing.AddBdrFaceIntegrator(mfem.BoundaryLFIntegrator(stopping), low_energy)
            outgoing.Assemble()
            self.assertAlmostEqual(outgoing.GetDataArray() @ psi.GetDataArray(), 1.0, delta=1e-8)

            velocity_x = VectorConstCoefficient(0.0, 0.0)
            velocity_E = mfem.ScalarVectorProductCoefficient(
                stopping, VectorConstCoefficient(0.0, -1.0))
            A = Bilinear_Form(fes, velocity_x, velocity_E, xs, *markers)
            rhs = Linear_Form(fes, mfem.ConstantCoefficient(0.0), exact,
                              velocity_x, velocity_E, *markers)
            assembled = Initial_Guess(fes, 0.0)
            GMRES_solver(A, rhs, assembled, p_level=0)
            np.testing.assert_allclose(assembled.GetDataArray(), psi.GetDataArray(),
                                       rtol=1e-8, atol=1e-10)
            residual = mfem.Vector(fes.GetVSize())
            A.Mult(psi, residual)
            self.assertLess(np.linalg.norm(residual.GetDataArray() - rhs.GetDataArray()), 1e-8)
        self.assertLess(errors[1], 1e-3)
        self.assertGreater(errors[0] / errors[1], 3.5)


if __name__ == '__main__':
    unittest.main()
