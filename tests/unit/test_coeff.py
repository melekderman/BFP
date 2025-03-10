import unittest
import numpy as np
from bfp.coeff import (
    TotalXSCoefficient,
    ScatteringXSCoefficient,
    StoppingPowerCoefficient,
    StoppingPowerDerivativeCoefficient,
    InflowCoefficientSN,
    QCoefficient,
    EnergyDependentCoefficient,
    XDependentCoefficient,
    VelocityCoefficientOld,
    VelocityCoefficient,
    ConstantCoefficient
)

class TestCoefficients(unittest.TestCase):

    def test_TotalXSCoefficient(self):
        xs_t_data = [1.0, 2.0, 3.0]
        coeff = TotalXSCoefficient(xs_t_data, 0.0, 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.0]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.5]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.9]), 3.0)

    def test_ScatteringXSCoefficient(self):
        xs_s_data = [0.5, 1.5, 2.5]
        coeff = ScatteringXSCoefficient(xs_s_data, 1.0, 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 1.0]), 0.5)
        self.assertAlmostEqual(coeff.EvalValue([0, 2.5]), 1.5)

    def test_StoppingPowerCoefficient(self):
        S_data = [10, 20, 30]
        coeff = StoppingPowerCoefficient(S_data, 2.0, 5.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 2.0]), 10)
        self.assertAlmostEqual(coeff.EvalValue([0, 4.9]), 30)

    def test_StoppingPowerDerivativeCoefficient(self):
        dS_data = [-1, -2, -3]
        coeff = StoppingPowerDerivativeCoefficient(dS_data, E_start=0.0, E_end=3.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.0]), -1)
        self.assertAlmostEqual(coeff.EvalValue([0, 2.9]), -3)

    def test_InflowCoefficientSN_positive_mu(self):
        coeff = InflowCoefficientSN(in_flux=5.0, mu=0.7)
        self.assertEqual(coeff.EvalValue([0, 0]), 5.0)
        self.assertEqual(coeff.EvalValue([10, 100]), 5.0)

    def test_InflowCoefficientSN_negative_mu(self):
        coeff = InflowCoefficientSN(in_flux=5.0, mu=-0.7)
        self.assertEqual(coeff.EvalValue([0, 0]), 0.0)
        self.assertEqual(coeff.EvalValue([10, 100]), 0.0)

    def test_QCoefficient_constant(self):
        coeff = QCoefficient(7.0)
        self.assertEqual(coeff.EvalValue([0, 0.5]), 7.0)

    def test_QCoefficient_energy_dependent(self):
        Q_data = [1.0, 2.0, 3.0, 4.0]
        coeff = QCoefficient(Q_data, 0.0, 4.0)
        self.assertEqual(coeff.EvalValue([0, 0.0]), 1.0)
        self.assertEqual(coeff.EvalValue([0, 3.9]), 4.0)

    def test_EnergyDependentCoefficient(self):
        data = [10.0, 20.0, 30.0]
        coeff = EnergyDependentCoefficient(data, 0, 3)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.0]), 10.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 2.9]), 30.0)

    def test_VelocityCoefficient(self):
        mu = 0.8
        S_arr = [1.0, 2.0, 3.0]
        E_start, E_end = 0.0, 3.0
        coeff = VelocityCoefficient(mu, S_arr, E_start, E_end)
        self.assertEqual(coeff.EvalValue([0, 0.0]), [0.8, 1.0])
        self.assertEqual(coeff.EvalValue([0, 2.9]), [0.8, 3.0])

    def test_ConstantCoefficient(self):
        const_value = 7.0
        coeff = ConstantCoefficient(const_value)
        self.assertEqual(coeff.EvalValue([0, 0.5]), 7.0)
        self.assertEqual(coeff.EvalValue([10, 100]), 7.0)

    def test_XDependentCoefficient(self):
        data = [1.0, 3.0, 5.0]
        coeff = XDependentCoefficient(data, 0.0, 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([1.5, 0]), 3.0)

    def test_VelocityCoefficientOld(self):
        mu = 0.5
        S_coeff = ConstantCoefficient(2.0)
        coeff = VelocityCoefficientOld(mu, S_coeff)
        self.assertEqual(coeff._EvalPy([0, 0]), [0.5, 2.0])
        self.assertEqual(coeff._EvalPy([1, 1]), [0.5, 2.0])

if __name__ == '__main__':
    unittest.main()
