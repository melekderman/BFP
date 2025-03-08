import unittest
import numpy as np
from bfp.coeff import (
    TotalXSCoefficient,
    ScatteringXSCoefficient,
    StoppingPowerCoefficient,
    StoppingPowerDerivativeCoefficient,
    InflowCoefficient,
    QCoefficient,
    EnergyDependentCoefficient,
    VelocityCoefficient
)

class IntegrationPointMock:
    def __init__(self, y):
        self.y = y
    def __getitem__(self, idx):
        return [0, self.y][idx]


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
        self.assertAlmostEqual(coeff.EvalValue([0, 0]), 0.5)
        self.assertAlmostEqual(coeff.EvalValue([0, 0.5]), 1.5)

    def test_StoppingPowerCoefficient(self):
        S_data = [10, 20, 30]
        coeff = StoppingPowerCoefficient(S_data, 2.0, 5.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0]), 10)
        self.assertAlmostEqual(coeff.EvalValue([0, 1]), 30)

    def test_StoppingPowerDerivativeCoefficient(self):
        dS_data = [-1, -2, -3]
        coeff = StoppingPowerDerivativeCoefficient(dS_data=dS_data, E_start=0.0, E_end=3.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 0]), -1)
        self.assertAlmostEqual(coeff.EvalValue([0, 1]), -3)

    def test_InflowCoefficient(self):
        coeff = InflowCoefficient(5.0)
        self.assertEqual(coeff.EvalValue([0, 0]), 5.0)
        self.assertEqual(coeff.EvalValue([10, 100]), 5.0)

    def test_QCoefficient_constant(self):
        coeff = QCoefficient(7.0)
        self.assertEqual(coeff.EvalValue([0, 0.5]), 7.0)

    def test_QCoefficient_energy_dependent(self):
        Q_data = [1.0, 2.0, 3.0, 4.0]
        coeff = QCoefficient(Q_data, 0.0, 4.0)
        self.assertEqual(coeff.EvalValue([0, 0]), 1.0)
        self.assertEqual(coeff.EvalValue([0, 1]), 4.0)

    def test_EnergyDependentCoefficient(self):
        data = [10.0, 20.0, 30.0]
        coeff = EnergyDependentCoefficient(data, 0, 3)
        self.assertAlmostEqual(coeff.EvalValue([0, 0]), 10.0)
        self.assertAlmostEqual(coeff.EvalValue([0, 1]), 30.0)

    def test_VelocityCoefficient(self):
        class MockIP:
            x = 0
            y = 0.5

            def __getitem__(self, idx):
                return [0, self.y][idx]

        mu = 1.0
        S_data = [2.0, 4.0]
        s_coeff = StoppingPowerCoefficient(S_data, 0.0, 2.0)
        coeff = VelocityCoefficient(mu, s_coeff)

        V = np.zeros(2)
        coeff._EvalPy(V, [0, 0.5])
        self.assertAlmostEqual(V[0], mu)
        self.assertAlmostEqual(V[1], 4.0)


if __name__ == '__main__':
    unittest.main()
