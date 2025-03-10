import mfem.ser as mfem
import numpy as np
import h5py

__all__ = [
    'read_data',
    'gauss_legendre_dirs',
    'compute_S_derivative',
    'gridfunction_to_array',
    'assign_boundary_attributes',
    'get_marker_for_mu'
    ]

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
        E (np.ndarray): 1D array of energy values.
        S (np.ndarray): 1D array of S values corresponding to the energy values.

    Returns:
        np.ndarray: Array of the finite difference approximations of dS/dE.
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


def gridfunction_to_array(phi):
    """Converts an MFEM GridFunction into a NumPy array.

    Args:
        phi_new: An MFEM GridFunction-like object with methods Size() and __getitem__.

    Returns:
        A NumPy array containing the grid function values.
    """
    return np.array([phi[i] for i in range(phi.Size())])


def assign_boundary_attributes(mesh, x_start, E_start, x_end, tol=1e-6):
    """
    Assigns boundary attributes to the mesh:
      - If any vertex of a boundary element is near (x_start, E_start), assign attribute 1.
      - Else if any vertex of a boundary element is near x_end, assign attribute 2.
      - Otherwise, leave the attribute as 0.
    
    Parameters:
      mesh    : The mfem.Mesh object.
      x_start : The x-coordinate for the left (inflow) boundary.
      E_start : The energy coordinate for the left boundary.
      x_end   : The x-coordinate for the right (outflow) boundary.
      tol     : Tolerance for comparing coordinates.
    """
    vertex_array = mesh.GetVertexArray()

    for i in range(mesh.GetNBE()):
        mesh.SetBdrAttribute(i, 0)

    for i in range(mesh.GetNBE()):
        v_indices = mesh.GetBdrElementVertices(i)
        for idx in v_indices:
            x_coord = vertex_array[idx][0]
            E_coord = vertex_array[idx][1]
            if abs(x_coord - x_start) < tol and abs(E_coord - E_start) < tol:
                mesh.SetBdrAttribute(i, 1)
                break
            elif abs(x_coord - x_end) < tol:
                mesh.SetBdrAttribute(i, 2)
                break


def get_marker_for_mu(mesh, mu):
    """
    Returns a marker array based on the mesh boundary attributes and the sign of mu.
      - If mu > 0, only boundaries with attribute 1 (left/inflow) are marked.
      - If mu < 0, only boundaries with attribute 2 (right/inflow) are marked.
    
    Parameters:
      mesh : The mfem.Mesh object (its boundary attributes should already be set).
      mu   : The angular variable.
    
    Returns:
      marker (mfem.intArray): An integer array with 1's where the appropriate boundary
                              (inflow) is active and 0's elsewhere.
    """
    num_bdr = mesh.bdr_attributes.Size()
    marker = mfem.intArray(num_bdr)
    marker.Assign(0)
    for i in range(num_bdr):
        attr = mesh.bdr_attributes[i]
        if mu > 0 and attr == 1:
            marker[i] = 1
        elif mu < 0 and attr == 2:
            marker[i] = 1
    return marker