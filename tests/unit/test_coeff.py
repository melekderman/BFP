import unittest
from bfp.coeff import (
    TotalXSCoefficientE,
    ScatteringXSCoefficientE,
    StoppingPowerCoefficientE,
    QCoefficientE,
    EDependentCoefficient,
    XDependentCoefficient,
    ConstantCoefficient,
)

class TestCoefficients(unittest.TestCase):

    # Test the TotalXSCoefficientE class.
    def test_TotalXSCoefficientE(self):
        data = [0,1,2,3,4]
        coeff = TotalXSCoefficientE(data, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 4.0)

    def test_TotalXSCoefficientE_constant(self):
        coeff = TotalXSCoefficientE(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 40.0)

    # Test the ScatteringXSCoefficientE class.
    def test_ScatteringXSCoefficientE(self):
        data = [0,1,2,3,4]
        coeff = ScatteringXSCoefficientE(data, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 4.0)

    def test_ScatteringXSCoefficientE_constant(self):
        coeff = ScatteringXSCoefficientE(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 40.0)

    # Test the StoppingPowerCoefficientE class.
    def test_StoppingPowerCoefficientE(self):
        data = [0,1,2,3,4]
        coeff = StoppingPowerCoefficientE(data, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 4.0)
    
    def test_StoppingPowerCoefficientE_constant(self):
        coeff = StoppingPowerCoefficientE(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 40.0)

    # Test the QCoefficientE class.
    def test_QCoefficientE(self):
        data = [0,1,2,3,4]
        coeff = QCoefficientE(data, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 4.0)
    
    def test_QCoefficientE_constant(self):
        coeff = QCoefficientE(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 40.0)

    # Test the EDependentCoefficient class.
    def test_EDependentCoefficient(self):
        data = [0,1,2,3,4]
        coeff = EDependentCoefficient(data, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 4.0)
    
    def test_EDependentCoefficient_constant(self):
        coeff = EDependentCoefficient(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.802]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.604]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.406]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.208]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.01]), 40.0)

    # Test the XDependentCoefficient class.
    def test_XDependentCoefficient(self):
        data = [0,1,2,3,4]
        coeff = XDependentCoefficient(data, 0.0, 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.0]), 0.0)
        self.assertAlmostEqual(coeff.EvalValue([0.208, 0.0]), 1.0)
        self.assertAlmostEqual(coeff.EvalValue([0.406, 0.0]), 2.0)
        self.assertAlmostEqual(coeff.EvalValue([0.604, 0.0]), 3.0)
        self.assertAlmostEqual(coeff.EvalValue([0.802, 0.0]), 4.0)
        self.assertAlmostEqual(coeff.EvalValue([1.0, 0.0]), 4.0)

    def test_XDependentCoefficient_constant(self):
        coeff = XDependentCoefficient(40, 1, 0.01)
        self.assertAlmostEqual(coeff.EvalValue([0.0, 0.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.208, 0.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.406, 0.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.604, 0.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([0.802, 0.0]), 40.0)
        self.assertAlmostEqual(coeff.EvalValue([1.0, 0.0]), 40.0)

    # Test the ConstantCoefficient class.
    def test_ConstantCoefficient(self):
        coeff = ConstantCoefficient(40.0)
        self.assertEqual(coeff.EvalValue([0.0, 1.0]), 40.0)
        self.assertEqual(coeff.EvalValue([1.0, 2.0]), 40.0)
        self.assertEqual(coeff.EvalValue([2.0, 3.0]), 40.0)

if __name__ == '__main__':
    unittest.main()
