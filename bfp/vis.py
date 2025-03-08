import mfem.ser as mfem
import matplotlib.pyplot as plt
import numpy as np
from glvis import glvis

__all__ = [
    'glvis_visualize',
    'matplotlib_visualize',
    ]

###############################################################################
# Visualizers:
###############################################################################


def glvis_visualize(mesh, solution, title="CSDA Transport Solution"):
    """Visualize the solution using GLVis.

    Args:
        mesh: The computational mesh.
        solution: The solution grid function.
        title (str): Title for the GLVis window.
    """
    ss = mfem.socketstream("localhost", 19916)
    if not ss.good():
        print("Unable to open GLVis socket connection.")
        return
    ss.send_text("solution\n" + title + "\n")
    mesh.Print(ss)
    solution.Save(ss)
    ss.send_text("\n")


def matplotlib_visualize(mesh, solution, E_range, nx, nE, title):
    """Visualize the solution using matplotlib.

    Assumes the mesh is a structured Cartesian mesh.
    The y-coordinate is mapped to energy using:
         E = E_start - y*(E_start - E_end).

    Args:
        mesh: The computational mesh.
        solution: The solution grid function.
        E_range (tuple): (E_end, E_start) for energy mapping.
        nx (int): Number of spatial elements.
        nE (int): Number of energy groups.
    """

    # Get the mesh nodes.
    nodes = mesh.GetNodes()
    if nodes.Size() == 0:
        print("Mesh has no nodes; cannot plot.")
        return
    
    # nodes is a mfem.Matrix with dimensions (dim x num_nodes).
    dim = nodes.Width()
    num_nodes = nodes.Height()

    # Extract x and y coordinates.
    x_coords = np.array([nodes.Get(i, 0) for i in range(num_nodes)])
    y_coords = np.array([nodes.Get(i, 1) for i in range(num_nodes)])

    # Map y-coordinate to energy.
    E_end, E_start = E_range
    E_values = E_start - y_coords * (E_start - E_end)

    # Get solution values at nodes.
    sol = np.array([solution.GetValue(i) for i in range(num_nodes)])
    
    # Create a scatter plot.
    plt.figure()
    sc = plt.scatter(x_coords, E_values, c=sol, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('Energy (MeV)')
    plt.title(title)
    plt.colorbar(sc, label='Flux')
    plt.show()