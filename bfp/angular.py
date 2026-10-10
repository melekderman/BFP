import mfem.ser as mfem
import numpy as np

__all__ = [
    "gauss_legendre_dirs",
    "integrate_angular_flux",
    "morel_edge_coefficients",
    "build_morel_angular_operator",
]


def gauss_legendre_dirs(N_dir):
    """Return ``N_dir`` Gauss–Legendre directions and weights as two lists.

    Directions increase on [-1, 1]; weights are positive and sum to 2.
    ``N_dir`` must be a positive integer. No angular-average normalization
    is applied. Also available as ``bfp.gauss_legendre_dirs`` and through utils.
    """
    mu, w = np.polynomial.legendre.leggauss(N_dir)
    return mu.tolist(), w.tolist()


def integrate_angular_flux(fes, psi_mu_list):
    """Return scalar flux as an MFEM GridFunction on ``fes``.

    Sum ``weight * psi`` from ``(mu, weight, GridFunction)`` entries in their
    supplied order, with no factor of 1/2. The function writes no files.
    """
    phi = mfem.GridFunction(fes)
    # PyMFEM 4.7 GridFunction does not retain its finite element space.
    phi._fes = fes
    phi.Assign(0.0)
    for _, w, X in psi_mu_list:
        phi.Add(w, X)
    return phi


def morel_edge_coefficients(mu, weights):
    """Return N+1 Morel edge coefficients as a NumPy array.

    Supply matching increasing full-range Gauss–Legendre ``mu`` and ``weights``.
    The recursion is ``alpha[0]=0``, ``alpha[n+1]=alpha[n]-2*mu[n]*weights[n]``.
    These dimensionless edge coefficients are unrelated to GFP2's material alpha.
    """
    mu = np.asarray(mu, dtype=float)
    weights = np.asarray(weights, dtype=float)
    alpha = np.zeros(mu.size + 1)
    for n in range(mu.size):
        alpha[n + 1] = alpha[n] - 2.0 * mu[n] * weights[n]
    return alpha


def build_morel_angular_operator(mu, weights):
    """Return the dense N-by-N Morel angular Laplace–Beltrami matrix.

    ``mu`` and ``weights`` are matching full-range Gauss–Legendre data in
    increasing direction order, with positive weights summing to 2. The
    dimensionless tridiagonal matrix approximates ``d_mu((1-mu**2)*d_mu(psi))``.
    Within floating-point tolerance, ``L*1=0``, ``L*mu=-2*mu``, ``weights.T*L=0``
    and the spectrum is nonpositive. Its coefficients already include the
    quadrature weights; do not multiply angular coupling blocks by them again.
    """
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
