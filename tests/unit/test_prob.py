import contextlib
import io
import os
import tempfile
import unittest
import warnings
from unittest.mock import patch

import mfem.ser as mfem
import numpy as np

from bfp import gauss_legendre_dirs, Solve_Transport
from bfp.angular import (
    gauss_legendre_dirs as angular_dirs,
    integrate_angular_flux,
    morel_edge_coefficients,
    build_morel_angular_operator,
)
from bfp.assemble import (
    BoundaryConditions,
    FESpace,
    Bilinear_Form,
    Linear_Form,
    Mass_Matrix,
    GFP2_Matrix,
)
from bfp.coeff import VectorConstCoefficient, gfp2_coefficients
from bfp.input import (
    problem1_input,
    problem2_input,
    problem3_input,
    problem4_input,
    problem5_input,
    slab_input,
)
from bfp.mesh import create_2D_mesh
from bfp.solver import (
    Solve_Psi,
    Solve_Phi,
    Solve_FP,
    Solve_GFP2,
    GMRES_solver,
    Initial_Guess,
)
from bfp.utils import gauss_legendre_dirs as legacy_dirs
from main import main


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
        display = patch("bfp.vis.GlVis_2D")
        display.start()
        self.addCleanup(display.stop)

    def test_transport_dispatch(self):
        mesh = mfem.Mesh(2, 3, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        mu, weights = angular_dirs(4)
        c = mfem.ConstantCoefficient
        markers = BoundaryConditions(mesh)
        common = dict(
            mu_vals=mu,
            w_vals=weights,
            mesh=mesh,
            fes=fes,
            S_coeff=AnalyticCoefficient(lambda x: 1.0 + x[1]),
            q_coeff=c(1.0),
            inflow=c(0.0),
            dir_bdr1=markers[0],
            dir_bdr2=markers[1],
            visualize=False,
            save=False,
            p_level=0,
        )
        cases = [
            (
                "csda",
                Solve_Psi,
                dict(
                    pn=1,
                    xs_t_coeff=c(1.0),
                    xs_t_const=1.0,
                    S_const=1.0,
                    alpha=1.0,
                    beta=0.5,
                    a_const=0.0,
                    b_const=0.0,
                ),
            ),
            ("fp", Solve_FP, dict(absorption_coeff=c(1.0), transport_coeff=c(1.0))),
            (
                "gfp2",
                Solve_GFP2,
                dict(
                    absorption_coeff=c(1.0),
                    alpha_coeff=c(0.6),
                    beta_coeff=c(0.1),
                ),
            ),
        ]
        for model, solve, coefficients in cases:
            for limit in (1000, 1):
                with (
                    self.subTest(model=model, limit=limit),
                    warnings.catch_warnings(record=True),
                ):
                    args = dict(common, **coefficients, iter_=limit)
                    direct, expected = solve(**args, return_info=True)
                    dispatched, info = Solve_Transport(model, **args, return_info=True)
                    default = Solve_Transport(model, **args)
                    self.assertEqual(info, expected)
                    diagnostics = info if model == "csda" else [info]
                    self.assertTrue(
                        all(d["converged"] == (limit == 1000) for d in diagnostics)
                    )
                    for a, b, d in zip(direct, dispatched, default):
                        self.assertEqual(a[:2], b[:2])
                        np.testing.assert_array_equal(
                            a[2].GetDataArray(), b[2].GetDataArray()
                        )
                        np.testing.assert_array_equal(
                            a[2].GetDataArray(), d[2].GetDataArray()
                        )
                    direct_phi = Solve_Phi(fes, direct, save=False)
                    dispatched_phi = Solve_Phi(fes, dispatched, save=False)
                    np.testing.assert_array_equal(
                        direct_phi.GetDataArray(), dispatched_phi.GetDataArray()
                    )
        self.assertEqual(os.listdir(), [])

    def test_transport_argument_validation(self):
        c = mfem.ConstantCoefficient(1.0)
        common = dict(
            mu_vals=[0.0],
            w_vals=[2.0],
            mesh=None,
            fes=None,
            absorption_coeff=c,
            S_coeff=c,
            q_coeff=c,
            inflow=c,
            dir_bdr1=None,
            dir_bdr2=None,
        )
        with patch("bfp.assemble.Mass_Matrix") as mass:
            with self.assertRaisesRegex(ValueError, "Unknown transport model"):
                Solve_Transport("unknown")
            for model, scattering in [
                ("fp", dict(transport_coeff=c)),
                ("gfp2", dict(alpha_coeff=c, beta_coeff=c)),
            ]:
                for name in scattering:
                    args = dict(common, **scattering)
                    del args[name]
                    with self.assertRaisesRegex(TypeError, name):
                        Solve_Transport(model, **args)
                    args[name] = None
                    with self.assertRaisesRegex(TypeError, name):
                        Solve_Transport(model, **args)
                with self.assertRaisesRegex(TypeError, "xs_t_coeff"):
                    Solve_Transport(model, **common, **scattering, xs_t_coeff=c)
            mass.assert_not_called()

    def test_cli_models_and_slab(self):
        numerical = ["--nx", "2", "--nE", "3", "--N_ang", "4", "--p_level", "0"]
        for model in ("fp", "gfp2"):
            for limit in (1000, 1):
                output = io.StringIO()
                argv = [
                    "bfp",
                    "--model",
                    model,
                    "-p",
                    "slab",
                    "--iter_",
                    str(limit),
                    *numerical,
                ]
                with (
                    patch("sys.argv", argv),
                    contextlib.redirect_stdout(output),
                    warnings.catch_warnings(record=True),
                ):
                    main()
                expected = (
                    "GMRES converged"
                    if limit == 1000
                    else "Saved results are unconverged"
                )
                self.assertIn(expected, output.getvalue())
        for problem in range(1, 6):
            outputs = []
            for model_args in ([], ["--model", "csda"]):
                output = io.StringIO()
                with (
                    patch(
                        "sys.argv", ["bfp", "-p", str(problem), *model_args, *numerical]
                    ),
                    contextlib.redirect_stdout(output),
                ):
                    main()
                outputs.append(output.getvalue().splitlines()[-1])
            self.assertEqual(outputs[0], outputs[1])
        invalid = [
            ["--model", "unknown", "-p", "slab"],
            ["--model", "csda", "-p", "slab"],
            ["--model", "fp", "-p", "1"],
            ["--model", "gfp2", "-p", "5"],
            ["--model", "fp"],
            ["-p", "6"],
        ]
        with patch("bfp.input.create_2D_mesh") as create_mesh:
            for args in invalid:
                with (
                    patch("sys.argv", ["bfp", *args]),
                    contextlib.redirect_stderr(io.StringIO()),
                ):
                    with self.assertRaises(SystemExit) as error:
                        main()
                    self.assertEqual(error.exception.code, 2)
            with self.assertRaisesRegex(ValueError, "slab example"):
                slab_input("csda")
            create_mesh.assert_not_called()

    def test_angular_integration(self):
        self.assertIs(gauss_legendre_dirs, angular_dirs)
        self.assertIs(legacy_dirs, angular_dirs)
        mu, weights = angular_dirs(4)
        self.assertIsInstance(mu, list)
        self.assertIsInstance(weights, list)
        self.assertAlmostEqual(sum(weights), 2.0)
        self.assertEqual(mu, sorted(mu))
        mesh = mfem.Mesh(1, 1, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        angular = []
        for direction, weight in zip(mu, weights):
            psi = Initial_Guess(fes, 1.0 + direction**2)
            angular.append((direction, weight, psi))
        phi = integrate_angular_flux(fes, angular)
        wrapped = Solve_Phi(fes, angular, save=False)
        np.testing.assert_allclose(phi.GetDataArray(), 8.0 / 3.0, rtol=0.0, atol=1e-14)
        np.testing.assert_array_equal(phi.GetDataArray(), wrapped.GetDataArray())
        self.assertEqual(os.listdir(), [])

    def test_morel_angular_operator(self):
        for N in (1, 2, 3, 4, 8, 16, 32):
            with self.subTest(N=N):
                mu, weights = angular_dirs(N)
                alpha = morel_edge_coefficients(mu, weights)
                L = build_morel_angular_operator(mu, weights)
                self.assertEqual(alpha.shape, (N + 1,))
                self.assertEqual(alpha[0], 0.0)
                self.assertAlmostEqual(alpha[-1], 0.0, delta=1e-14)
                self.assertTrue(np.all(alpha[1:-1] > 0.0))
                np.testing.assert_allclose(L @ np.ones(N), 0.0, rtol=0.0, atol=1e-12)
                np.testing.assert_allclose(
                    L @ mu, -2.0 * np.asarray(mu), rtol=0.0, atol=1e-12
                )
                np.testing.assert_allclose(
                    np.asarray(weights) @ L, 0.0, rtol=0.0, atol=1e-12
                )
                eigenvalues = np.linalg.eigvals(L)
                np.testing.assert_allclose(eigenvalues.imag, 0.0, rtol=0.0, atol=1e-12)
                self.assertLessEqual(np.max(eigenvalues.real), 1e-12)

    def test_gmres_diagnostics(self):
        A = mfem.SparseMatrix(4)
        for i in range(4):
            A.Set(i, i, 4.0)
            if i:
                A.Set(i, i - 1, -1.0)
                A.Set(i - 1, i, -1.0)
        A.Finalize()
        rhs = mfem.Vector([1.0, 2.0, 3.0, 4.0])
        for limit, converged in ((1000, True), (1, False)):
            with (
                self.subTest(limit=limit),
                warnings.catch_warnings(record=True) as caught,
            ):
                warnings.simplefilter("always", RuntimeWarning)
                default = mfem.Vector(4)
                default.Assign(0.0)
                self.assertIs(
                    GMRES_solver(A, rhs, default, iter_=limit, p_level=0), default
                )
                solution = mfem.Vector(4)
                solution.Assign(0.0)
                returned, info = GMRES_solver(
                    A, rhs, solution, iter_=limit, p_level=0, return_info=True
                )
                self.assertIs(returned, solution)
                self.assertEqual(info["converged"], converged)
                self.assertEqual(info["iterations"], 4 if converged else 1)
                np.testing.assert_array_equal(
                    solution.GetDataArray(), default.GetDataArray()
                )
                residual = mfem.Vector(4)
                A.Mult(solution, residual)
                residual -= rhs
                preconditioned = mfem.Vector(4)
                mfem.GSSmoother(A).Mult(residual, preconditioned)
                self.assertAlmostEqual(
                    info["final_norm"], preconditioned.Norml2(), delta=1e-13
                )
                self.assertEqual(len(caught), 0 if converged else 2)
                for warning in caught:
                    self.assertIs(warning.category, RuntimeWarning)
                    self.assertIn("GMRES did not converge", str(warning.message))

    def test_fp_zero_scattering(self):
        mesh = mfem.Mesh(2, 3, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        mu, weights = angular_dirs(4)
        markers = BoundaryConditions(mesh)
        c = mfem.ConstantCoefficient
        legacy = Solve_Psi(
            5,
            mu,
            weights,
            mesh,
            fes,
            c(2.0),
            2.0,
            c(1.0),
            1.0,
            1.0,
            0.5,
            *markers,
            1.0,
            2.0,
            p_level=0,
            q_coeff=c(3.0),
            visualize=False,
            save=False,
        )
        flux, info = Solve_FP(
            mu,
            weights,
            mesh,
            fes,
            absorption_coeff=c(2.0),
            transport_coeff=c(0.0),
            S_coeff=c(1.0),
            q_coeff=c(3.0),
            inflow=c(1.0),
            dir_bdr1=markers[0],
            dir_bdr2=markers[1],
            p_level=0,
            visualize=False,
            save=False,
            return_info=True,
        )
        self.assertTrue(info["converged"])
        self.assertLess(info["relative_true_residual"], 1e-10)
        for (_, _, old), (_, _, new) in zip(legacy, flux):
            np.testing.assert_allclose(
                old.GetDataArray(), new.GetDataArray(), atol=1e-10, rtol=1e-10
            )

    def test_fp_gfp2_variable_coefficient_mms(self):
        mesh = mfem.Mesh(2, 2, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(2, mesh)
        mu, weights = angular_dirs(4)
        markers = BoundaryConditions(mesh)
        h = lambda x: 1.0 + x[0] + x[1]
        stopping = lambda x: 1.0 + x[1]
        absorption = lambda x: 1.0 + 0.1 * x[0]
        alpha = lambda x: 0.6 + 0.1 * x[0] + 0.2 * x[1]
        beta = lambda x: 0.1 + 0.05 * (x[0] + x[1])
        for solve in (Solve_FP, Solve_GFP2):
            with self.subTest(solve=solve.__name__):
                b = (lambda x: 0.0) if solve is Solve_FP else beta
                db = 0.0 if solve is Solve_FP else 0.05
                exact, sources = [], []
                for m in mu:
                    exact.append(
                        AnalyticCoefficient(
                            lambda x, m=m: h(x) * (2.0 + (1.0 + 2.0 * b(x)) * m)
                        )
                    )
                    sources.append(
                        AnalyticCoefficient(
                            lambda x, m=m: (m - stopping(x))
                            * (2.0 + (1.0 + 2.0 * b(x)) * m + 2.0 * db * m * h(x))
                            + (absorption(x) - 1.0)
                            * h(x)
                            * (2.0 + (1.0 + 2.0 * b(x)) * m)
                            + 2.0 * alpha(x) * m * h(x)
                        )
                    )
                scattering = (
                    dict(transport_coeff=AnalyticCoefficient(lambda x: 2.0 * alpha(x)))
                    if solve is Solve_FP
                    else dict(
                        alpha_coeff=AnalyticCoefficient(alpha),
                        beta_coeff=AnalyticCoefficient(beta),
                    )
                )
                flux, info = solve(
                    mu,
                    weights,
                    mesh,
                    fes,
                    **scattering,
                    absorption_coeff=AnalyticCoefficient(absorption),
                    S_coeff=AnalyticCoefficient(stopping),
                    q_coeff=sources,
                    inflow=exact,
                    dir_bdr1=markers[0],
                    dir_bdr2=markers[1],
                    p_level=0,
                    visualize=False,
                    save=False,
                    return_info=True,
                )
                self.assertTrue(info["converged"])
                self.assertLess(info["relative_true_residual"], 1e-10)
                for (_, _, psi), analytic in zip(flux, exact):
                    self.assertLess(psi.ComputeL2Error(analytic), 1e-9)
                phi = Solve_Phi(fes, flux, save=False)
                self.assertLess(
                    phi.ComputeL2Error(AnalyticCoefficient(lambda x: 4.0 * h(x))), 1e-9
                )

    def test_gfp2_schur_complement_and_fp_limit(self):
        mesh = mfem.Mesh(1, 1, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        mu, weights = angular_dirs(4)
        L = build_morel_angular_operator(mu, weights)
        c = mfem.ConstantCoefficient
        mass = Mass_Matrix(fes, c(1.0))
        transport = [
            Bilinear_Form(
                fes,
                VectorConstCoefficient(m, 0.0),
                VectorConstCoefficient(0.0, -1.0),
                c(1.0),
                *BoundaryConditions(mesh),
            )
            for m in mu
        ]
        a, b = gfp2_coefficients(0.8, 1.6)
        K = a * np.linalg.solve(np.eye(4) - b * L, L)
        matrix = GFP2_Matrix(
            transport, L, Mass_Matrix(fes, c(a)), Mass_Matrix(fes, c(b)), mass
        )
        dense_matrix = matrix.ToDenseMatrix()
        dense = dense_matrix.GetDataArray()
        size = 4 * fes.GetVSize()
        schur = dense[:size, :size] - dense[:size, size:] @ np.linalg.solve(
            dense[size:, size:], dense[size:, :size]
        )
        dense_mass = mass.ToDenseMatrix()
        expected = dense[:size, :size] - np.kron(K, dense_mass.GetDataArray())
        np.testing.assert_allclose(schur, expected, atol=1e-14, rtol=1e-14)

        common = dict(
            absorption_coeff=c(1.0),
            S_coeff=c(1.0),
            q_coeff=c(2.0),
            inflow=[c(1.0 + m) for m in mu],
            dir_bdr1=BoundaryConditions(mesh)[0],
            dir_bdr2=BoundaryConditions(mesh)[1],
            visualize=False,
            save=False,
            p_level=0,
        )
        fp = Solve_FP(mu, weights, mesh, fes, transport_coeff=c(0.8), **common)
        errors = []
        for b in (1e-2, 1e-4, 0.0):
            flux = Solve_GFP2(
                mu,
                weights,
                mesh,
                fes,
                alpha_coeff=c(0.4 * (1.0 + 2.0 * b)),
                beta_coeff=c(b),
                **common,
            )
            errors.append(
                max(
                    np.max(np.abs(p.GetDataArray() - f.GetDataArray()))
                    for (_, _, p), (_, _, f) in zip(fp, flux)
                )
            )
        self.assertLess(errors[1], errors[0] / 50.0)
        self.assertLess(errors[2], 1e-10)

    def test_gfp2_angular_properties_and_coefficients(self):
        a, b = gfp2_coefficients(0.8, 1.6)
        self.assertGreater(a, 0.0)
        self.assertGreaterEqual(b, 0.0)
        for N in (2, 4, 8, 16):
            mu, weights = np.polynomial.legendre.leggauss(N)
            L = build_morel_angular_operator(mu, weights)
            K = a * np.linalg.solve(np.eye(N) - b * L, L)
            np.testing.assert_allclose(K @ np.ones(N), 0.0, atol=1e-13)
            np.testing.assert_allclose(weights @ K, 0.0, atol=1e-13)
            np.testing.assert_allclose(K @ mu, -0.8 * mu, atol=1e-13)
            eig = np.linalg.eigvals(K)
            self.assertLess(np.max(eig.real), 1e-13)
            self.assertLess(np.max(np.abs(eig.imag)), 1e-13)
        for first, second in (
            (0.0, 0.0),
            (-1.0, 2.0),
            (1.0, 1.0),
            (1.0, 0.5),
            (1.0, 4.0),
            (np.nan, 2.0),
        ):
            with self.assertRaises(ValueError):
                gfp2_coefficients(first, second)
        self.assertEqual(gfp2_coefficients(1.0, 3.0), (0.5, 0.0))
        mesh = mfem.Mesh(1, 1, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        for value, positive in ((-0.1, False), (0.0, True), (np.inf, True)):
            with self.assertRaises(ValueError):
                Mass_Matrix(fes, mfem.ConstantCoefficient(value), positive=positive)
        c = mfem.ConstantCoefficient
        with self.assertRaisesRegex(ValueError, "GSSmoother requires nonzero diagonal"):
            Solve_GFP2(
                *angular_dirs(3),
                mesh,
                fes,
                absorption_coeff=c(0.0),
                alpha_coeff=c(0.5),
                beta_coeff=c(0.1),
                S_coeff=c(0.0),
                q_coeff=c(1.0),
                inflow=c(0.0),
                dir_bdr1=BoundaryConditions(mesh)[0],
                dir_bdr2=BoundaryConditions(mesh)[1],
                visualize=False,
                save=False,
                p_level=0,
            )

    def test_angle_diagnostics_and_cli(self):
        for limit, converged in ((1000, True), (1, False)):
            with self.subTest(limit=limit), warnings.catch_warnings(record=True):
                phi, info = problem5_input(
                    nx=3, nE=4, N_ang=2, iter_=limit, p_level=0, return_info=True
                )
                mu, _ = np.polynomial.legendre.leggauss(2)
                np.testing.assert_array_equal([entry["mu"] for entry in info], mu)
                self.assertTrue(all(entry["converged"] == converged for entry in info))
                self.assertTrue(np.isfinite(phi.GetDataArray()).all())
                output = io.StringIO()
                argv = [
                    "main.py",
                    "-p",
                    "5",
                    "--nx",
                    "3",
                    "--nE",
                    "4",
                    "--N_ang",
                    "2",
                    "--iter_",
                    str(limit),
                    "--p_level",
                    "0",
                ]
                with patch("sys.argv", argv), contextlib.redirect_stdout(output):
                    main()
                if converged:
                    self.assertIn(
                        "GMRES converged for all angular directions", output.getvalue()
                    )
                else:
                    self.assertIn("Saved results are unconverged", output.getvalue())
                    self.assertNotIn(
                        "GMRES converged for all angular directions", output.getvalue()
                    )

    def test_optional_output(self):
        mesh = mfem.Mesh(2, 3, "QUADRILATERAL", True, 1.0, 1.0)
        fes = FESpace(1, mesh)
        mu, weights = np.polynomial.legendre.leggauss(2)
        xs = mfem.ConstantCoefficient(5.0)
        args = (
            5,
            mu,
            weights,
            mesh,
            fes,
            xs,
            5.0,
            5.0,
            2.0,
            1.0,
            0.5,
            *BoundaryConditions(mesh),
            5.0,
            10.0,
        )

        with patch("bfp.vis.GlVis_2D") as display:
            angular = Solve_Psi(*args, p_level=0, visualize=False, save=False)
            phi = Solve_Phi(fes, angular, save=False)
            display.assert_not_called()
            self.assertEqual(os.listdir(), [])

            with_output = Solve_Psi(*args, p_level=0)
            saved_phi = Solve_Phi(fes, with_output)
            self.assertEqual(display.call_count, len(mu))
        self.assertEqual(
            set(os.listdir()),
            {"phi_dg.gf", *("psi_hs_mu_{:.3f}.gf".format(m) for m in mu)},
        )
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
        psi = Solve_Psi(
            1,
            [0.0],
            [1.0],
            mesh,
            fes,
            xs,
            sigma,
            inflow,
            1.0,
            1.0,
            0.5,
            *markers,
            0.0,
            0.0,
            p_level=0,
            S_coeff=stopping,
        )[0][2]
        return mesh, fes, markers, xs, psi

    def test_zero_stopping_power(self):
        exact = AnalyticCoefficient(
            lambda x: 4.0
            + 8.0
            * (
                np.exp(-5.0 * np.sqrt(3.0) * x[0])
                + np.exp(-5.0 * np.sqrt(3.0) * (1.0 - x[0]))
            )
        )
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
                nE, mfem.ConstantCoefficient(1.0), 1.0, inflow
            )
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
            outgoing.AddBdrFaceIntegrator(
                mfem.BoundaryLFIntegrator(stopping), low_energy
            )
            outgoing.Assemble()
            self.assertAlmostEqual(
                outgoing.GetDataArray() @ psi.GetDataArray(), 1.0, delta=1e-8
            )

            velocity_x = VectorConstCoefficient(0.0, 0.0)
            velocity_E = mfem.ScalarVectorProductCoefficient(
                stopping, VectorConstCoefficient(0.0, -1.0)
            )
            A = Bilinear_Form(fes, velocity_x, velocity_E, xs, *markers)
            rhs = Linear_Form(
                fes,
                mfem.ConstantCoefficient(0.0),
                exact,
                velocity_x,
                velocity_E,
                *markers,
            )
            assembled = Initial_Guess(fes, 0.0)
            GMRES_solver(A, rhs, assembled, p_level=0)
            np.testing.assert_allclose(
                assembled.GetDataArray(), psi.GetDataArray(), rtol=1e-8, atol=1e-10
            )
            residual = mfem.Vector(fes.GetVSize())
            A.Mult(psi, residual)
            self.assertLess(
                np.linalg.norm(residual.GetDataArray() - rhs.GetDataArray()), 1e-8
            )
        self.assertLess(errors[1], 1e-3)
        self.assertGreater(errors[0] / errors[1], 3.5)


if __name__ == "__main__":
    unittest.main()
