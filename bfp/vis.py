import mfem.ser as mfem
import matplotlib.pyplot as plt
import numpy as np
from glvis import glvis
from mfem.ser import socketstream

__all__ = [
    "GlVis_Visualizer",
    "Matplotlib_Visualizer",
    "Mesh_Report",
    "GlVis_2D",
    "GlVis_3D",
    "visualize_sol",
]


def GlVis_Visualizer(mesh, solution, title="CSDA Transport Solution"):
    """Send an MFEM ``mesh`` and solution GridFunction to GLVis on localhost:19916.

    Requires a running external GLVis server. ``title`` is the window label.
    Prints a message if the socket connection fails; returns None.
    """
    ss = mfem.socketstream("localhost", 19916)
    if not ss.good():
        print("Unable to open GLVis socket connection.")
        return
    ss.send_text("solution\n" + title + "\n")
    mesh.Print(ss)
    solution.Save(ss)
    ss.send_text("\n")


def Matplotlib_Visualizer(mesh, solution, E_range, nx, nE, title):
    """Legacy nodal plotting helper; not verified for general DG GridFunctions.

    Assumes a structured mesh with accessible nodal coordinates and matching
    solution values. The second coordinate is physical kinetic energy in MeV.
    ``E_range``, ``nx`` and ``nE`` are retained but unused; ``title`` labels the
    plot. Displays a Matplotlib figure and returns None. For current DG output,
    prefer :func:`GlVis_2D`; this helper may fail on meshes without nodal data.
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

    E_values = y_coords

    # Get solution values at nodes.
    sol = np.array([solution.GetValue(i) for i in range(num_nodes)])

    # Create a scatter plot.
    plt.figure()
    sc = plt.scatter(x_coords, E_values, c=sol, cmap="viridis")
    plt.xlabel("x")
    plt.ylabel("Energy (MeV)")
    plt.title(title)
    plt.colorbar(sc, label="Flux")
    plt.show()


def Mesh_Report(mesh, width=600, height=600, num_vertices_to_print=None):
    """Print vertex/element counts and display an MFEM mesh with PyGLVis.

    ``width`` and ``height`` are display dimensions in pixels.
    ``num_vertices_to_print`` defaults to the first ten vertices. Rendering
    requires a compatible notebook frontend. Returns None.
    """
    try:
        total_vertices = mesh.GetNV()
        print(f"Total number of vertices: {total_vertices}")
    except AttributeError as e:
        print("Error retrieving number of vertices:", e)

    try:
        total_elements = mesh.GetNE()
        print(f"Total number of elements: {total_elements}")
    except AttributeError as e:
        print("Error retrieving number of elements:", e)

    count = 10 if num_vertices_to_print is None else int(num_vertices_to_print)
    print(f"Number of vertices to display: {count}")

    try:
        vertices = mesh.GetVertexArray()
        print("Sample vertex coordinates:")
        for idx, vertex in enumerate(vertices):
            print(f"Vertex {idx}: {vertex}")
            if idx + 1 >= count:
                print("...")
                break
    except AttributeError as e:
        print("Error retrieving vertex coordinates:", e)

    try:
        print("Launching glvis visualization...")
        glvis(mesh, width, height)
    except ImportError as e:
        print("glvis module is not available:", e)
    except RuntimeError as e:
        print("An error occurred during visualization:", e)


def GlVis_2D(mesh, solution):
    """Render a two-dimensional MFEM ``mesh`` and solution GridFunction with PyGLVis.

    Requires a compatible notebook frontend. Uses a 1000-by-1000 display and
    returns None. This is the visualization invoked by the transport solvers.
    """
    g = glvis((mesh, solution), 500, 500, keys="aRjLGmcbp")
    g.set_size(1000, 1000)
    g.render()


def GlVis_3D(mesh, solution):
    """Render a three-dimensional MFEM ``mesh`` and solution GridFunction with PyGLVis.

    Requires a compatible notebook frontend; returns None. This visualization
    helper does not provide a three-dimensional transport solver.
    """
    g = glvis((mesh, solution), 500, 500, keys="ArljmGac//0./")
    g.set_size(1000, 1000)
    g.render()

    # This function was taken from the PyMFEM


# https://github.com/mfem/PyMFEM/blob/06e1efa9b0bc5bc5a13246e2e73d85badecf60ad/miniapps/gslib/field_diff.py#L174
def visualize_sol(
    sol,
    mesh,
    dim,
    title: str = "Mesh and solution",
    vishost: str = "localhost",
    visport: int = 19916,
    geometry: str = "0 0 600 600",
    keys: str = "ArljmGac//0./",
):
    """Legacy incomplete socket helper; use GlVis_2D or GlVis_Visualizer instead.

    The current body references undefined variables and does not correctly use
    its supplied solution/connection arguments. It is retained for compatibility
    but is not a functioning visualization entry point.
    """
    sol_sock = socketstream("localhost", 19916)
    sol_sock.precision(8)
    sol_sock << "solution\n" << mesh << x
    sol_sock << "window_title 'Solution'\n"
    sol.sock.flush()
    flux_sock = socketstream("localhost", 19916)
    flux_sock.precision(8)
    flux_sock << "solution\n" << mesh << flux
    flux_sock << "keys vvv\n"
    flux_sock << "window_geometry 402 0 400 350\n"
    flux_sock << "window_title 'Flux'\n"
    flux_sock.flush()
