import mfem.ser as mfem
import matplotlib.pyplot as plt
import numpy as np
from glvis import glvis
import seaborn as sns
import pandas as pd


__all__ = [
    'glvis_visualize',
    'matplotlib_visualize',
    'mesh_report',
    'HeatmapPlot',
    'ScatterPlot'
    ]

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


def mesh_report(mesh, width=600, height=600, num_vertices_to_print=None):
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


class HeatmapPlot:
    """Generates heatmaps from MFEM GridFunction data for 2D and 3D problems."""

    def __init__(self, phi_new, mesh, fes):
        """Initializes HeatmapPlot with MFEM objects.

        Args:
            phi_new: MFEM GridFunction object.
            mesh: MFEM mesh object.
            fes: MFEM FiniteElementSpace object.
        """
        self.phi_new = phi_new
        self.mesh = mesh
        self.fes = fes
        self.phi_arr = np.array([phi_new[i] for i in range(phi_new.Size())])

    def _compute_cell_centers(self, start, end, n_cells):
        """Computes cell-center coordinates.

        Args:
            start (float): Start coordinate.
            end (float): End coordinate.
            n_cells (int): Number of cells.

        Returns:
            np.ndarray: Cell-center coordinates.
        """
        nodes = np.linspace(start, end, n_cells + 1)
        return 0.5 * (nodes[:-1] + nodes[1:])

    def plot_heatmap(self, x_start, x_end, nx, E_start, E_end, nE,
                     y_start=None, y_end=None, ny=None, y_slice_idx=None,
                     vmin=None, vmax=None):
        """Plots a heatmap for either a 2D problem or a 2D slice of a 3D problem.

        Args:
            x_start, x_end (float): x-axis boundaries.
            nx (int): Number of cells in x-axis.
            E_start, E_end (float): Energy boundaries.
            nE (int): Number of cells in energy axis.
            y_start, y_end (float, optional): y-axis boundaries (for 3D).
            ny (int, optional): Number of cells in y-axis (for 3D).
            y_slice_idx (int, optional): y-slice index (for 3D visualization).
            vmin, vmax (float, optional): Color scale limits.
        """
        cell_x = self._compute_cell_centers(x_start, x_end, nx)
        cell_E = self._compute_cell_centers(E_start, E_end, nE)

        if y_start is not None and y_end is not None and ny is not None:
            cell_y = self._compute_cell_centers(y_start, y_end, ny)
            cell_avg = np.zeros((nE, ny, nx))
            cell_counter = 0
            for iE in range(nE):
                for iy in range(ny):
                    for ix in range(nx):
                        dof_indices = self.fes.GetElementDofs(cell_counter)
                        vals = [self.phi_arr[i] for i in dof_indices]
                        cell_avg[iE, iy, ix] = np.mean(vals)
                        cell_counter += 1
            if y_slice_idx is None:
                raise ValueError("y_slice_idx must be provided for 3D plots.")
            plot_data = cell_avg[:, y_slice_idx, :]
            title = f"Scalar Flux Φ(x,E) at y = {cell_y[y_slice_idx]:.4f}"
        else:
            cell_avg = np.zeros((nE, nx))
            cell_counter = 0
            for iE in range(nE):
                for ix in range(nx):
                    dof_indices = self.fes.GetElementDofs(cell_counter)
                    vals = [self.phi_arr[i] for i in dof_indices]
                    cell_avg[iE, ix] = np.mean(vals)
                    cell_counter += 1
            plot_data = cell_avg
            title = "Scalar Flux Φ(x,E)"

        df = pd.DataFrame(plot_data, index=np.round(cell_E, 4), columns=np.round(cell_x, 4))

        plt.figure(figsize=(8, 6))
        sns.heatmap(df, cmap="jet", cbar_kws={'label': 'Scalar Flux Φ'},
                    vmin=vmin, vmax=vmax, annot=False, fmt=".2f",
                    linewidths=0.01, linecolor='grey')
        plt.xlabel("x (cm)")
        plt.ylabel("Energy (MeV)")
        plt.title(title)
        plt.grid(True)
        plt.show()


class ScatterPlot:
    """Generates scatter plots from MFEM GridFunction data."""

    def __init__(self, phi_new, mesh, fes):
        """Initializes ScatterPlot with MFEM objects.

        Args:
            phi_new: MFEM GridFunction object.
            mesh: MFEM mesh object.
            fes: MFEM FiniteElementSpace object.
        """
        self.phi_new = phi_new
        self.mesh = mesh
        self.fes = fes
        self.dim = mesh.Dimension()

    def plot_scatter(self, cmap='jet', point_size=50, vmin=None, vmax=None):
        """Plots a scatter plot showing element-averaged scalar flux at cell centers.

        Args:
            cmap (str, optional): Matplotlib colormap name.
            point_size (int, optional): Size of scatter plot points.
            vmin, vmax (float, optional): Color scale limits.
        """
        num_cells = self.mesh.GetNE()
        cell_centers = np.zeros((num_cells, self.dim))
        phi_cell_avg = np.zeros(num_cells)
        center = mfem.Vector(self.dim)

        phi_arr = np.array([self.phi_new[i] for i in range(self.phi_new.Size())])

        for i in range(num_cells):
            self.mesh.GetElementCenter(i, center)
            cell_centers[i, :] = [center[j] for j in range(self.dim)]
            dof_indices = self.fes.GetElementDofs(i)
            phi_vals = [phi_arr[dof] for dof in dof_indices]
            phi_cell_avg[i] = np.mean(phi_vals)

        plt.figure(figsize=(8, 6))
        sc = plt.scatter(cell_centers[:, 0], cell_centers[:, 1], c=phi_cell_avg,
                         cmap=cmap, s=point_size, vmin=vmin, vmax=vmax)
        plt.colorbar(sc, label='Scalar Flux Φ')
        plt.xlabel("x (cm)")
        plt.ylabel("Energy (MeV)")
        plt.title("Element-Averaged Scalar Flux Distribution Φ(x,E)")
        plt.grid(True)
        plt.show()