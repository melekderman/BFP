import mfem.ser as mfem
import numpy as np

__all__ = [
    "gauss_legendre_dirs",
    "integrate_angular_flux",
    "morel_edge_coefficients",
    "build_morel_angular_operator",
]


def gauss_legendre_dirs(N_dir):
    """Return Gauss-Legendre directions and weights on [-1, 1] as lists."""
    mu, w = np.polynomial.legendre.leggauss(N_dir)
    return mu.tolist(), w.tolist()


def integrate_angular_flux(fes, psi_mu_list):
    """Sum (mu, weight, GridFunction) entries in their supplied order."""
    phi = mfem.GridFunction(fes)
    # PyMFEM 4.7 GridFunction does not retain its finite element space.
    phi._fes = fes
    phi.Assign(0.0)
    for _, w, X in psi_mu_list:
        phi.Add(w, X)
    return phi


def morel_edge_coefficients(mu, weights):
    """Return recursive edge coefficients for full-range Gauss-Legendre data."""
    mu = np.asarray(mu, dtype=float)
    weights = np.asarray(weights, dtype=float)
    alpha = np.zeros(mu.size + 1)
    for n in range(mu.size):
        alpha[n + 1] = alpha[n] - 2.0 * mu[n] * weights[n]
    return alpha


def build_morel_angular_operator(mu, weights):
    """Return Morel's Laplace-Beltrami matrix for increasing Gauss-Legendre nodes."""
    mu = np.asarray(mu, dtype=float)
    weights = np.asarray(weights, dtype=float)
    alpha = morel_edge_coefficients(mu, weights)
    N = mu.size
    L = np.zeros((N, N))
    for n in range(N):
        if n + 1 < N:
            c_plus = alpha[n + 1] / (weights[n] * (mu[n + 1] - mu[n]))
            L[n, n + 1] += c_plus
            L[n, n] -= c_plus
        if n - 1 >= 0:
            c_minus = alpha[n] / (weights[n] * (mu[n] - mu[n - 1]))
            L[n, n - 1] += c_minus
            L[n, n] -= c_minus
    return L
