import mfem.ser as mfem
import numpy as np

__all__ = ['gauss_legendre_dirs', 'integrate_angular_flux']


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
