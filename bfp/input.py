# input.py
import os
os.environ['DISPLAY'] = ':0'

from mfem.ser import mfem
from mfem.ser import socketstream
from glvis import glvis
from glvis.widget import GlvisWidget
from ipywidgets.embed import embed_minimal_html

from bfp.mesh import create_2D_mesh
from bfp.assemble import FESpace, BoundaryConditions
from bfp.coeff import TotalXSCoefficientE
from bfp.utils import gauss_legendre_dirs
from bfp.solver import Solve_Psi, Solve_Phi
from bfp.vis import visualize_sol

__all__ = ['problem1_input',
        'problem2_input',
        'problem3_input',
        'problem4_input',
        'problem5_input',
        ]


def problem1_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-12, p_level=1):
    """
    Executes Problem 1: Infinite Medium Problem, ψ = Q/σₜ

    Args:
        nx (int): Number of cells in the x-direction.
        nE (int): Number of energy cells.
        N_ang (int): Number of angles for the SN method.
        iter_ (int): Maximum number of solver iterations.
        tol (float): Solver tolerance.
        p_level (int): Print level (1 for verbose, 0 for silent).

    Returns:
        phi: The computed solution.
    """

    # Set the problem parameters
    x_start = 0.0
    x_end = 1.0
    E_start = 0.0
    E_end = 1.0
    order = 1
    a = 0
    b = 0
    alpha = 1.0
    beta = 0.5

    # Problem 1: ψ = Q/σₜ, 
    pn = 1
    inflow = 20.0
    S_const = 5.0
    xs_t_const = 5.0
    q_const = 100.0

    # Create mesh
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)
    dim = mesh.Dimension()
    # Define finite element space
    fes = FESpace(order, mesh)

    # Define boundary attributes
    dir_bdr1, dir_bdr2 = BoundaryConditions(mesh, x_min=0, x_max=0, y_min=1, y_max=1)

    # Set coefficients and SN quadratures
    xs_t_coeff = TotalXSCoefficientE(xs_t_const, E_start, E_end)
    mu_vals, w_vals = gauss_legendre_dirs(N_ang)

    # Compute solution: Calculate psi then phi
    psi_mu_pos_list = Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const,
                                 inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2,
                                 a, b, q_const, E_start, E_end, iter_, tol, p_level)
    phi = Solve_Phi(fes, psi_mu_pos_list)

    # Visualize the solution
    flux_sock = socketstream("localhost", 19916)
    flux_sock.precision(8)
    flux_sock.send_text("solution\n")
    flux_sock.send_solution(mesh, phi)
    flux_sock.send_text("\n")
    flux_sock.flush()
    widget = GlvisWidget((mesh, phi), 600, 600, keys="ArljmGac//0")
    widget.render()

    return phi

def problem2_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-12, p_level=1):
    """
    Executes Problem 2: S=0 Problem,  ψ(x) = Q/σₜ + ψₗ * exp(-σₜ * x / μ).

    Args:
        nx (int): Number of cells in the x-direction.
        nE (int): Number of energy cells.
        N_ang (int): Number of angles for the SN method.
        iter_ (int): Maximum number of solver iterations.
        tol (float): Solver tolerance.
        p_level (int): Print level (1 for verbose, 0 for silent).

    Returns:
        phi: The computed solution.
    """
    # Set the problem parameters
    x_start = 0.0
    x_end = 1.0
    E_start = 0.0
    E_end = 1.0
    order = 1
    a = 0
    b = 0
    alpha = 1.0
    beta = 0.5

    # Problem 2: ψ(x) = Q/σₜ + ψₗ * exp(-σₜ * x / μ) 
    pn = 2
    inflow = 10.0
    S_const = 0.0
    xs_t_const = 5.0
    q_const = 10.0

    # Create mesh and apply uniform refinement
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)

    # Define finite element space
    fes = FESpace(order, mesh)

    # Define boundary attributes (for Problem 2, boundaries are set differently)
    dir_bdr1, dir_bdr2 = BoundaryConditions(mesh, x_min=1, x_max=1, y_min=1, y_max=1)

    # Set coefficients and SN quadratures
    xs_t_coeff = TotalXSCoefficientE(xs_t_const, E_start, E_end)
    mu_vals, w_vals = gauss_legendre_dirs(N_ang)

    # Compute solution: Calculate psi then phi
    psi_mu_pos_list = Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const,
                                 inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2,
                                 a, b, q_const, E_start, E_end, iter_, tol, p_level)
    phi = Solve_Phi(fes, psi_mu_pos_list)

    return phi

def problem3_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-12, p_level=1):
    """
    Executes Problem 3: MMS Problem,  ψ = a + b * x.
    
    Parameters:
        nx (int): Number of cells in the x-direction.
        nE (int): Number of energy cells.
        N_ang (int): Number of angles for the SN method.
        iter_ (int): Maximum number of solver iterations.
        tol (float): Solver tolerance.
        p_level (int): Print level (1 for verbose output, 0 for silent operation).
    
    Returns:
        phi: The computed solution for Problem 3.
    """
    # Set the problem parameters
    x_start = 0.0
    x_end = 1.0
    E_start = 0.0
    E_end = 1.0
    order = 1
    a = 5
    b = 10
    alpha = 1.0
    beta = 0.5

    # Problem 3: ψ = a+bx
    pn=3
    inflow = 5.0
    S_const = 0.0
    xs_t_const = 5.0
    q_const = 0.0

    # Create mesh and apply uniform refinement.
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)

    # Define the finite element space.
    fes = FESpace(order, mesh)

    # Define boundary attributes.
    dir_bdr1, dir_bdr2 = BoundaryConditions(mesh, x_min=0, x_max=0, y_min=1, y_max=1)

    # Set coefficients and compute SN quadratures.
    xs_t_coeff = TotalXSCoefficientE(xs_t_const, E_start, E_end)
    mu_vals, w_vals = gauss_legendre_dirs(N_ang)

    # Compute solution: first calculate psi, then calculate phi.
    psi_mu_pos_list = Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const,
                                 inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2,
                                 a, b, q_const, E_start, E_end, iter_, tol, p_level)
    phi = Solve_Phi(fes, psi_mu_pos_list)

    return phi

def problem4_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-12, p_level=1):
    """
    Executes Problem 4: MMS Problem, ψ = a + b * E.
    
    Parameters:
        nx (int): Number of cells in the x-direction.
        nE (int): Number of energy cells.
        N_ang (int): Number of angles for the SN method.
        iter_ (int): Maximum number of solver iterations.
        tol (float): Solver tolerance.
        p_level (int): Print level (1 for verbose output, 0 for silent operation).
    
    Returns:
        phi: The computed solution for Problem 4.
    """
    # Set the problem parameters
    x_start = 0.0
    x_end = 1.0
    E_start = 0.0
    E_end = 1.0
    order = 1
    a = 5
    b = 10
    alpha = 1.0
    beta = 0.5

    # Problem 4: ψ = a+bE
    pn=4
    inflow = 5.0
    S_const = 2.0
    xs_t_const = 5.0
    q_const = 0.0

    # Create mesh and apply uniform refinement.
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)

    # Define the finite element space.
    fes = FESpace(order, mesh)

    # Define boundary attributes.
    dir_bdr1, dir_bdr2 = BoundaryConditions(mesh, x_min=0, x_max=0, y_min=1, y_max=1)

    # Set coefficients and compute SN quadratures.
    xs_t_coeff = TotalXSCoefficientE(xs_t_const, E_start, E_end)
    mu_vals, w_vals = gauss_legendre_dirs(N_ang)

    # Compute solution: first calculate psi, then calculate phi.
    psi_mu_pos_list = Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const,
                                 inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2,
                                 a, b, q_const, E_start, E_end, iter_, tol, p_level)
    phi = Solve_Phi(fes, psi_mu_pos_list)

    return phi

def problem5_input(nx=10, nE=10, N_ang=2, iter_=1000, tol=1e-12, p_level=1):
    """
    Executes Problem 5: ψ = a + b * xE.
    
    Args:
        nx (int): Number of cells in the x-direction.
        nE (int): Number of energy cells.
        N_ang (int): Number of angles for the SN method.
        iter_ (int): Maximum number of solver iterations.
        tol (float): Solver tolerance.
        p_level (int): Print level (1 for verbose, 0 for silent).
    
    Returns:
        phi: The computed solution.
    """
    # Set the problem parameters
    x_start = 0.0
    x_end = 1.0
    E_start = 0.0
    E_end = 1.0
    order = 1
    a = 5
    b = 10
    alpha = 1.0
    beta = 0.5

    # Problem 5: ψ = a+bxE
    pn=5
    inflow = 5.0
    S_const = 2.0
    xs_t_const = 5.0

    # Create mesh and apply uniform refinement
    mesh = create_2D_mesh(nx, nE, x_start, x_end, E_start, E_end)

    # Define finite element space
    fes = FESpace(order, mesh)

    # Define boundary attributes
    dir_bdr1, dir_bdr2 = BoundaryConditions(mesh, x_min=0, x_max=0, y_min=1, y_max=1)

    # Set coefficients and SN quadratures
    xs_t_coeff = TotalXSCoefficientE(xs_t_const, E_start, E_end)
    mu_vals, w_vals = gauss_legendre_dirs(N_ang)

    # Compute solution: Calculate psi and then phi
    psi_mu_pos_list = Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const,
                                 inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2,
                                 a, b, iter_, tol, p_level)
    phi = Solve_Phi(fes, psi_mu_pos_list)

    return phi



