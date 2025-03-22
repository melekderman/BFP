#!/usr/bin/env python
"""
CSDA Transport Equation Solver with Angular Discretization

Solves the equation:
    Ω ∂ψ(x,E)/∂x + ∂[S(E) ψ(x,E)]/∂E + Σ_t(E) ψ(x,E) = Q(x,E)
in a 2D (x, E) domain using discrete ordinates.
The angular variable (μ) is discretized using a 6‑point Gauss–Legendre rule
over μ ∈ [–1, 1] (yielding both right- and left-moving directions).
The spatial domain is [0, 0.3] (20 cells) and the energy domain is represented by 50 cells,
with energy decreasing from 1.0 MeV to 0.01 MeV.
Cross-section and stopping power data are read from an HDF5 file.
"""

import mfem.ser as mfem
import numpy as np
import h5py
import math
import matplotlib.pyplot as plt
from glvis import glvis

###############################################################################
# MeshGenerator: Create a 2D mesh in (x, E)
###############################################################################
def create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end):
    """
    Create a 2D Cartesian mesh for the (x,E) domain.

    The x-coordinate spans [x_start, x_end] (nx cells).
    The second coordinate is a normalized variable in [0,1] that is mapped to energy
    via E = E_end + y*(E_start - E_end) (nE cells).

    Args:
        nx (int): Number of cells in x-direction.
        nE (int): Number of cells in energy direction.
        x_start (float): Starting x value.
        x_end (float): Ending x value.
        E_start (float): Highest energy value.
        E_end (float): Lowest energy value.

    Returns:
        mfem.Mesh: Generated 2D mesh.
    """
    mesh = mfem.Mesh(nx, nE, "QUADRILATERAL", True, 0.0, 0.0)
    x_coords = np.linspace(x_start, x_end, nx + 1)
    y_coords = np.linspace(0, 1, nE + 1)  # normalized energy coordinate

    verts = mesh.GetVertexArray()
    k = 0
    for j in range(nE + 1):
        for i in range(nx + 1):
            verts[k][0] = x_coords[i]
            verts[k][1] = y_coords[j]
            k += 1
    return mesh

###############################################################################
# DataReader: Read cross-section and stopping power data from HDF5
###############################################################################
def read_data(nE):
    """
    Read cross-section and stopping power data from an HDF5 file.

    Expects datasets "xs_t_{nE}", "xs_s_{nE}", and "S_{nE}".

    Args:
        nE (int): Number of energy groups.

    Returns:
        tuple: (xs_t_arr, xs_s_arr, S_arr) as numpy arrays.
    """
    data_file = "data.h5"
    with h5py.File(data_file, "r") as f:
        xs_t_arr = f[f"xs_t_{nE}"][:]
        xs_s_arr = f[f"xs_s_{nE}"][:]
        S_arr    = f[f"S_{nE}"][:]
    return xs_t_arr, xs_s_arr, S_arr

###############################################################################
# Coefficient Wrappers
###############################################################################
class CrossSectionCoefficient(mfem.PyCoefficient):
    """
    Coefficient for the total cross-section Σ_t(E).

    Maps the normalized energy coordinate (x[1]) to the appropriate cross-section value.
    """
    def __init__(self, xs_t_data, E_start, E_end):
        super(CrossSectionCoefficient, self).__init__()
        self.xs_t_data = xs_t_data
        self.E_start = E_start
        self.E_end = E_end

    def EvalValue(self, x):
        y = x[1]  # normalized energy coordinate in [0,1]
        E = self.E_start + y * (self.E_start - self.E_end)
        n_groups = len(self.xs_t_data)
        group = min(n_groups - 1, int((E - self.E_end) / (self.E_start - self.E_end) * n_groups))
        return float(self.xs_t_data[group])

class StoppingPowerCoefficient(mfem.PyCoefficient):
    """
    Coefficient for the stopping power S(E).

    Maps the normalized energy coordinate (x[1]) to the corresponding S(E) value.
    """
    def __init__(self, S_data, E_start, E_end):
        super(StoppingPowerCoefficient, self).__init__()
        self.S_data = S_data
        self.E_start = E_start
        self.E_end = E_end

    def EvalValue(self, x):
        y = x[1]
        E = self.E_start + y * (self.E_start - self.E_end)
        n_groups = len(self.S_data)
        group = min(n_groups - 1, int((E - self.E_end) / (self.E_start - self.E_end) * n_groups))
        return float(self.S_data[group])

def Q_function(x):
    """
    Source term Q(x,E). Returns 0.0 in this example.

    Args:
        x (mfem.Vector): The coordinate vector.
        
    Returns:
        float: 0.0.
    """
    return 0.0

def inflow_function(x):
    """
    Inflow boundary condition function.

    Imposes unit flux at the spatial boundary (x=0) when the physical energy is near 1.0 MeV.

    Args:
        x (mfem.Vector): The coordinate vector.
        
    Returns:
        float: Inflow flux value.
    """
    E_start = 1.0
    E_end = 0.01
    y = x[1]
    E = E_end + y * (E_start - E_end)
    if abs(x[0]) < 1e-3 and abs(E - 1.0) < 1e-3:
        return 1.0
    return 0.0

###############################################################################
# Vector Coefficient Wrappers for Convection Integrators
###############################################################################
class AngularVectorCoefficient(mfem.VectorPyCoefficient):
    """
    A vector coefficient for the angular term.

    Returns a constant vector [mu, 0] representing the angular direction
    in the (x,E) space.
    """
    def __init__(self, mu, dim=2):
        super(AngularVectorCoefficient, self).__init__(dim)
        self.mu = mu
        self.dim = dim

    def EvalValue(self, x):
        """
        Evaluate the angular vector coefficient at a point x.

        Args:
            x (mfem.Vector): The coordinate vector.
        
        Returns:
            list: A list representing the vector [mu, 0] in 2D (or [mu, 0, 0] in 3D).
        """
        if self.dim == 1:
            return [self.mu]
        elif self.dim == 2:
            return [self.mu, 0.0]
        elif self.dim == 3:
            return [self.mu, 0.0, 0.0]

class EnergyVectorCoefficient(mfem.VectorPyCoefficient):
    """
    A vector coefficient that wraps a scalar coefficient for energy convection.

    Given a scalar coefficient (for example, stopping power S(E)),
    returns a vector [0, S(E)] in 2D.
    """
    def __init__(self, scalar_coeff, dim=2):
        super(EnergyVectorCoefficient, self).__init__(dim)
        self.scalar_coeff = scalar_coeff
        self.dim = dim

    def EvalValue(self, x):
        """
        Evaluate the energy vector coefficient at a point x.

        Args:
            x (mfem.Vector): The coordinate vector.
        
        Returns:
            list: A list representing the vector [0, S(E)] in 2D (or [0, S(E), 0] in 3D).
        """
        S_val = self.scalar_coeff.EvalValue(x)
        if self.dim == 1:
            return [S_val]
        elif self.dim == 2:
            return [0.0, S_val]
        elif self.dim == 3:
            return [0.0, S_val, 0.0]


###############################################################################
# Assembler: Assemble the DG weak form for one discrete angle
###############################################################################
def assemble_system_for_angle(fes, mesh, E_range, xs_t_arr, S_arr, mu, inflow_coeff):
    """
    Assemble the DG system for the CSDA transport equation for a given angular direction.

    For a fixed discrete angle μ, the weak form is:
        μ ∂ψ/∂x + ∂(S(E)ψ)/∂E + Σ_t(E)ψ = Q(x,E)
    with DG interface integrals for upwinding.

    Args:
        fes (mfem.FiniteElementSpace): Finite element space.
        mesh (mfem.Mesh): The computational mesh.
        E_range (tuple): (E_start, E_end) for energy.
        xs_t_arr (numpy.array): Total cross-section data.
        S_arr (numpy.array): Stopping power data.
        mu (float): Discrete angular cosine.
        inflow_coeff (mfem.Coefficient): Inflow boundary coefficient.

    Returns:
        tuple: (A, b) where A is the system matrix and b is the RHS vector.
    """
    E_start, E_end = E_range

    # Create physical scalar coefficients.
    xs_t_coeff = CrossSectionCoefficient(xs_t_arr, E_start, E_end)
    S_scalar_coeff = StoppingPowerCoefficient(S_arr, E_start, E_end)
    # Use constant coefficients for Q (source = 0) and for inflow, use the provided inflow_coeff.
    Q_coeff = mfem.ConstantCoefficient(0.0)

    # Create vector coefficients for convection integrators.
    # For the spatial convection term, wrap mu in a vector [mu, 0].
    mu_vector_coeff = AngularVectorCoefficient(mu, dim=2)
    # For the energy convection term, wrap the scalar S(E) into a vector [0, S(E)].
    S_vector_coeff = EnergyVectorCoefficient(S_scalar_coeff, dim=1)

    # Build the bilinear form.
    a = mfem.BilinearForm(fes)
    # Spatial convection term: μ ∂ψ/∂x.
    a.AddDomainIntegrator(mfem.ConvectionIntegrator(mu_vector_coeff, 1.0))
    # Energy convection term: ∂(S(E)ψ)/∂E.
    a.AddDomainIntegrator(mfem.ConvectionIntegrator(S_vector_coeff, 1.0))
    # Absorption term: Σ_t(E)ψ.
    a.AddDomainIntegrator(mfem.MassIntegrator(xs_t_coeff))
    # Add DG interior face integrators for upwind flux.
    a.AddInteriorFaceIntegrator(mfem.DGTraceIntegrator(mu_vector_coeff, 1.0, -0.5))
    a.AddInteriorFaceIntegrator(mfem.DGTraceIntegrator(S_vector_coeff, 1.0, -0.5))
    a.Assemble()
    A = a.SpMat()

    # Assemble the right-hand side.
    b = mfem.LinearForm(fes)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(Q_coeff))
    b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(inflow_coeff))
    b.Assemble()

    return A, b

###############################################################################
# TransportSolver: Solve the system for a given discrete angle
###############################################################################
class TransportSolver:
    """
    Solver for the CSDA transport equation for a given discrete angle.

    Provides a direct LU solver.
    """
    def __init__(self, A, b, fes):
        """
        Initialize the solver.

        Args:
            A (mfem.SparseMatrix): The system matrix.
            b (mfem.Vector): The right-hand side.
            fes (mfem.FiniteElementSpace): The finite element space.
        """
        self.A = A
        self.b = b
        self.fes = fes

    def solve_lu(self):
        """
        Solve the system using a direct LU solver.

        Returns:
            mfem.Vector: The solution vector.
        """
        x = mfem.Vector(self.b.Size())
        lu_solver = mfem.SuperLUSolver()
        lu_solver.SetOperator(self.A)
        lu_solver.Mult(self.b, x)
        return x

###############################################################################
# Angular Discretization: Get discrete angles from Gauss-Legendre quadrature
###############################################################################
def get_angular_directions():
    """
    Retrieve discrete angular directions (μ) and weights using a 6-point Gauss-Legendre rule.

    Returns:
        list of tuples: Each tuple is (mu, weight) with μ in [-1, 1].
    """
    IR = mfem.IntRules.Get(mfem.Geometry.SEGMENT, 6)
    directions = []
    for i in range(IR.GetNPoints()):
        ip = IR.IntPoint(i)
        directions.append((ip.x, ip.weight))
    return directions

###############################################################################
# Visualizer: Plot the solution using GLVis and Matplotlib
###############################################################################
def glvis_visualize(mesh, solution, title="CSDA Transport Scalar Flux"):
    """
    Visualize the solution using GLVis.

    Args:
        mesh (mfem.Mesh): The mesh.
        solution (mfem.GridFunction): The solution.
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

def matplotlib_visualize(mesh, solution, E_range, nx, nE):
    """
    Visualize the solution using Matplotlib.

    Maps the normalized energy coordinate y to physical energy:
      E = E_end + y*(E_start - E_end)
    
    Args:
        mesh (mfem.Mesh): The mesh.
        solution (mfem.GridFunction): The solution.
        E_range (tuple): (E_start, E_end).
        nx (int): Number of spatial cells.
        nE (int): Number of energy cells.
    """
    nodes = mesh.GetNodes()
    if nodes.Size() == 0:
        print("Mesh has no nodes; cannot plot.")
        return
    num_nodes = nodes.Height()
    x_coords = np.array([nodes.Get(i, 0) for i in range(num_nodes)])
    y_coords = np.array([nodes.Get(i, 1) for i in range(num_nodes)])
    E_start, E_end = E_range
    E_values = E_end + y_coords * (E_start - E_end)
    sol = np.array([solution.GetValue(i) for i in range(num_nodes)])
    plt.figure()
    sc = plt.scatter(x_coords, E_values, c=sol, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('Energy (MeV)')
    plt.title('CSDA Transport Scalar Flux')
    plt.colorbar(sc, label='Flux')
    plt.show()

###############################################################################
# Main driver: Loop over discrete angles, solve, and integrate the solution
###############################################################################
def main():
    # Problem parameters
    nx = 20
    nE = 50
    x_start = 0.0
    x_end = 0.3
    E_start = 1.0   # highest energy
    E_end = 0.01    # lowest energy
    E_range = (E_start, E_end)

    # Read data from HDF5.
    xs_t_arr, xs_s_arr, S_arr = read_data(nE)

    # Create mesh.
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)
    
    # Define finite element space (L2, order=1).
    order = 1
    fec = mfem.L2_FECollection(order, mesh.Dimension())
    fes = mfem.FiniteElementSpace(mesh, fec)
    print("Number of unknowns:", fes.GetVSize())

    # Define inflow coefficient using a PyCoefficient subclass.
    class InflowCoefficient(mfem.PyCoefficient):
        """
        Inflow boundary coefficient that returns 1.0 when x[0] ≈ 0 and E ≈ 1.0 MeV.
        """
        def __init__(self):
            super(InflowCoefficient, self).__init__()

        def EvalValue(self, x):
            E_start = 1.0
            E_end = 0.01
            y = x[1]
            E = E_end + y * (E_start - E_end)
            if abs(x[0]) < 1e-3 and abs(E - 1.0) < 1e-3:
                return 1.0
            return 0.0

    inflow_coeff = InflowCoefficient()

    # Get discrete angular directions (μ) and weights.
    directions = get_angular_directions()  # list of (mu, weight)
    
    # Initialize a GridFunction to accumulate the angular-integrated (scalar) flux.
    scalar_flux = mfem.GridFunction(fes)
    scalar_flux.Assign(0.0)
    
    # Loop over each discrete angle.
    for (mu, weight) in directions:
        print("Solving for mu =", mu, "with weight =", weight)
        # Assemble the system for this angle.
        A, b = assemble_system_for_angle(fes, mesh, E_range, xs_t_arr, S_arr, mu, inflow_coeff)
        # Solve the system using direct LU.
        solver = TransportSolver(A, b, fes)
        sol = solver.solve_lu()
        # Convert the solution vector to a GridFunction.
        sol_gf = mfem.GridFunction(fes)
        sol_gf.Assign(sol)
        # Accumulate the weighted solution.
        scalar_flux.Add(weight, sol_gf)
    
    # Optionally normalize the scalar flux by the total weight.
    total_weight = sum(w for (_, w) in directions)
    if total_weight > 0:
        scalar_flux.Scale(1.0 / total_weight)
    
    # Visualize the angular-integrated (scalar) flux.
    glvis_visualize(mesh, scalar_flux, "CSDA Transport Scalar Flux")
    matplotlib_visualize(mesh, scalar_flux, E_range, nx, nE)

if __name__ == "__main__":
    main()



