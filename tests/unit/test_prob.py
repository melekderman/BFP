import unittest
from unittest.mock import patch

from bfp.input import problem1_input, problem2_input, problem3_input, problem4_input, problem5_input


class TestProblems(unittest.TestCase):

    def test_problem1(self):

        phi = problem1_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-8, p_level=1)
        self.assertIsNotNone(phi, "Problem 1 did not return a valid solution.")

    def test_problem2(self):
        phi = problem2_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-8, p_level=1)
        self.assertIsNotNone(phi, "Problem 2 did not return a valid solution.")

    def test_problem3(self):
        phi = problem3_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-8, p_level=1)
        self.assertIsNotNone(phi, "Problem 3 did not return a valid solution.")

    def test_problem4(self):
        phi = problem4_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-8, p_level=1)
        self.assertIsNotNone(phi, "Problem 4 did not return a valid solution.")

    def test_problem5(self):
        phi = problem5_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-8, p_level=1)
        self.assertIsNotNone(phi, "Problem 5 did not return a valid solution.")

if __name__ == "__main__":
    unittest.main()
