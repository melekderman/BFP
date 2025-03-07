import mfem.ser as mfem
import numpy as np
import h5py

__all__ = [
    'read_data',
    'gauss_legendre_dirs',
    'compute_S_derivative']

###############################################################################
# Utilities:
###############################################################################

def read_data(nE):
    
    """Read cross-section and stopping power data from an HDF5 file.

    Args:
        nE (int): Number of energy groups.

    Returns:
        tuple: A tuple containing:
            - E_arr (numpy.ndarray): Energy array.
            - E_grid_arr (numpy.ndarray): Energy grid array.
            - xs_t_arr (numpy.ndarray): Total cross-section array.
            - xs_s_arr (numpy.ndarray): Scattering cross-section array.
            - S_arr (numpy.ndarray): Stopping power array.
    """
    
    data_file = "data.h5"

    with h5py.File(data_file, "r") as f:
        E_arr = f[f"E_{nE}"][:]
        E_grid_arr = f[f"E_grid_{nE}"][:]
        xs_t_arr = f[f"xs_t_{nE}"][:]
        xs_s_arr = f[f"xs_s_{nE}"][:]
        S_arr    = f[f"S_{nE}"][:]
        
    return E_arr, E_grid_arr, xs_t_arr, xs_s_arr, S_arr


def gauss_legendre_dirs(N_dir):

    """Compute Gauss-Legendre quadrature points and weights for discrete ordinates.

    This function computes the Gauss-Legendre quadrature points (mu_i) and weights (w_i)
    for the specified number of discrete ordinates using NumPy's polynomial module.

    Args:
        N_dir (int): The number of discrete ordinates.

    Returns:
        tuple: A tuple containing two lists:
            - mu (list of float): The Gauss-Legendre quadrature points.
            - w (list of float): The corresponding Gauss-Legendre quadrature weights.
    """

    mu, w = np.polynomial.legendre.leggauss(N_dir)
    return mu.tolist(), w.tolist()


def compute_S_derivative(E, S):
    """Compute the derivative dS/dE using finite differences.

    This function approximates the derivative of S with respect to E by applying finite differences.
    It uses a forward difference at the first point, a backward difference at the last point, and
    central differences for all interior points.

    Args:
        E (numpy.ndarray): 1D array of energy values.
        S (numpy.ndarray): 1D array of S values corresponding to the energy values.

    Returns:
        numpy.ndarray: Array of the finite difference approximations of dS/dE.
    """
    n = len(S)
    dS_dE = np.zeros(n)
    
    # Forward difference for the first point
    dS_dE[0] = (S[1] - S[0]) / (E[1] - E[0])
    
    # Central difference for interior points
    for i in range(1, n - 1):
        dS_dE[i] = (S[i + 1] - S[i - 1]) / (E[i + 1] - E[i - 1])
    
    # Backward difference for the last point
    dS_dE[-1] = (S[-1] - S[-2]) / (E[-1] - E[-2])
    
    return dS_dE

