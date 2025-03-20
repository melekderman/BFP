import mfem.ser as mfem
import matplotlib.pyplot as plt
import numpy as np
from glvis import glvis
import seaborn as sns
import pandas as pd


__all__ = [
    'GlVis_Visualizer',
    'Matplotlib_Visualizer',
    'Mesh_Report',
    'GlVis_2D',
    'GlVis_3D',
    ]

def GlVis_Visualizer(mesh, solution, title="CSDA Transport Solution"):
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


def Matplotlib_Visualizer(mesh, solution, E_range, nx, nE, title):
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


def Mesh_Report(mesh, width=600, height=600, num_vertices_to_print=None):
    """
    Print mesh information and visualize the mesh using glvis.

    This function prints the total number of vertices and elements,
    displays a specified number of vertex coordinates (defaulting to the first 10 if not provided),
    and then launches the glvis visualization window with the specified dimensions.

    Parameters:
        mesh : mfem Mesh object
            The mesh object containing vertices and elements.
        width : int, optional
            Width of the glvis window (default is 600).
        height : int, optional
            Height of the glvis window (default is 600).
        num_vertices_to_print : int, optional
            The number of vertex coordinates to print. If not provided, the first 10 vertices will be printed.
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
    g = glvis((mesh, solution), 500,500, keys='aRjLGmcbp')
    g.set_size(1000, 1000)
    g.render()

def GlVis_3D(mesh, solution):
    g = glvis((mesh, solution), 500,500, keys='ArljmGac//0./')
    g.set_size(1000, 1000)
    g.render()