import warnings

import mfem.ser as mfem
import numpy as np
from bfp import angular, assemble, vis, coeff

__all__ = ['GMRES_solver',
            'Initial_Guess',
            'Solve_Phi',
            'Solve_Psi',
            'Solve_FP',
            'Solve_GFP2',
           ]

def GMRES_solver(A, b, psi, iter_=1000, tol=1e-12, p_level=1, *, return_info=False):
    """Optionally return (psi, info), with MFEM's preconditioned residual norm."""
    prec = mfem.GSSmoother(A)
    solver = mfem.GMRESSolver()
    solver.SetOperator(A)
    solver.SetPreconditioner(prec)
    solver.SetRelTol(tol)
    solver.SetAbsTol(tol)
    solver.SetMaxIter(iter_)
    solver.SetKDim(30)
    solver.SetPrintLevel(p_level)
    solver.Mult(b, psi)

    info = dict(converged=solver.GetConverged(),
                iterations=solver.GetNumIterations(),
                final_norm=solver.GetFinalNorm())
    if not info['converged']:
        warnings.warn(
            f"GMRES did not converge after {info['iterations']} iterations "
            f"(reported residual norm: {info['final_norm']:.6g}).",
            RuntimeWarning, stacklevel=2)
    if return_info:
        return psi, info
    return psi


def Initial_Guess(fes, psi_const):
    psi = mfem.GridFunction(fes)
    psi._fes = fes
    psi.Assign(psi_const)
    
    return psi


def Solve_Phi(fes, psi_mu_list, *, save=True):
    """Integrate angular flux; set save=False to skip the scalar flux file."""
    phi = angular.integrate_angular_flux(fes, psi_mu_list)
    if save:
        phi.Save("phi_dg.gf")

    return phi


def Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const, inflow, S_const, 
              alpha, beta, dir_bdr1, dir_bdr2, a_const, b_const, q_const=0, E_start=0, E_end=0, iter_=1000, tol=1e-12, p_level=1,
              *, S_coeff=None, q_coeff=None, visualize=True, save=True, return_info=False):
    """Solve in physical (x, E); S_coeff overrides positive S_const.

    S_coeff is a continuous MFEM scalar coefficient. For variable-S MMS,
    supply q_coeff including -d(S*psi)/dE and the exact inflow coefficient.
    Both energy faces must be marked: high-E inflow and low-E outflow.
    Set visualize=False and save=False to disable GLVis and flux files.
    With return_info=True, return (angular fluxes, per-angle diagnostics).
    """
    if S_coeff is not None and pn in (3, 4, 5) and q_coeff is None:
        raise ValueError("Variable-S manufactured solutions require q_coeff.")
    source_coeff = q_coeff
    if S_coeff is None:
        S_coeff = mfem.ConstantCoefficient(S_const)
    energy_direction = coeff.VectorConstCoefficient(0.0, -1.0)
    v_coeff2 = mfem.ScalarVectorProductCoefficient(S_coeff, energy_direction)

    psi_mu_list = []
    info_list = []

    for mu, w in zip(mu_vals, w_vals):
        a = a_const
        b = b_const
        print("  Solving for mu =", mu)
        inflow_coeff = inflow if isinstance(inflow, mfem.Coefficient) else coeff.InflowCoefficient(
            inflow, mu, pn=pn, a_val=a, b_val=b)
        v_coeff1 = coeff.VectorConstCoefficient(mu, 0.0)
        if source_coeff is not None:
            q_coeff = source_coeff
        elif pn == 1 or pn == 2:
            q_coeff = coeff.QCoefficientE(q_const, E_start, E_end)
        else:
            q_coeff = coeff.QFuncCoefficient(pn, a_val=a, b_val=b, xs_t_val=xs_t_const,
                                         mu_val=mu, S_val=S_const)

        A = assemble.Bilinear_Form(fes, v_coeff1, v_coeff2, xs_t_coeff,
                                  dir_bdr1, dir_bdr2, alpha=alpha, beta=beta)
        b = assemble.Linear_Form(fes, q_coeff, inflow_coeff, v_coeff1, v_coeff2,
                                 dir_bdr1, dir_bdr2, alpha=alpha, beta=beta)

        psi = Initial_Guess(fes, 1.0)
        psi, info = GMRES_solver(A, b, psi, iter_, tol, p_level, return_info=True)
        if visualize:
            vis.GlVis_2D(mesh, psi)
        if save:
            psi.Save("psi_hs_mu_{:.3f}.gf".format(mu))
        psi_mu_list.append((mu, w, psi))
        info_list.append(dict(mu=float(mu), **info))

    if return_info:
        return psi_mu_list, info_list
    return psi_mu_list


def Solve_FP(mu_vals, w_vals, mesh, fes, *, absorption_coeff, transport_coeff,
             S_coeff, q_coeff, inflow, dir_bdr1, dir_bdr2,
             iter_=1000, tol=1e-12, p_level=1,
             visualize=True, save=True, return_info=False):
    """Solve FP-CSDA with genuine absorption and elastic transport cross section.

    Supply MFEM coefficients; q_coeff and inflow may also be per-angle lists.
    S_coeff uses physical energy; mark both energy faces as in Solve_Psi.
    return_info adds global GMRES diagnostics and the true coupled residual.
    """
    mass_tr = assemble.Mass_Matrix(fes, transport_coeff)
    return _solve_coupled(mu_vals, w_vals, mesh, fes, absorption_coeff, S_coeff,
                          q_coeff, inflow, dir_bdr1, dir_bdr2, mass_tr, None, None,
                          iter_, tol, p_level, visualize, save, return_info)


def Solve_GFP2(mu_vals, w_vals, mesh, fes, *, absorption_coeff, alpha_coeff,
               beta_coeff, S_coeff, q_coeff, inflow, dir_bdr1, dir_bdr2,
               iter_=1000, tol=1e-12, p_level=1,
               visualize=True, save=True, return_info=False):
    """Solve GFP2-CSDA with alpha > 0, beta >= 0 and an auxiliary DG variable.

    Arguments/output follow Solve_FP. Coefficients may vary in (x,E); align
    material interfaces with the mesh. chi has no independent inflow data.
    gfp2_coefficients converts elastic relaxation cross sections to alpha,beta.
    With S=absorption=0, use even angular order for the GSSmoother.
    """
    mass_alpha = assemble.Mass_Matrix(fes, alpha_coeff, positive=True)
    mass_beta = assemble.Mass_Matrix(fes, beta_coeff)
    mass0 = assemble.Mass_Matrix(fes, mfem.ConstantCoefficient(1.0))
    return _solve_coupled(mu_vals, w_vals, mesh, fes, absorption_coeff, S_coeff,
                          q_coeff, inflow, dir_bdr1, dir_bdr2, mass_alpha, mass_beta, mass0,
                          iter_, tol, p_level, visualize, save, return_info)


def _solve_coupled(mu_vals, w_vals, mesh, fes, absorption_coeff, S_coeff,
                   q_coeff, inflow, dir_bdr1, dir_bdr2, mass, mass_beta, mass0,
                   iter_, tol, p_level, visualize, save, return_info):
    N = len(mu_vals)
    if N == 0 or len(w_vals) != N:
        raise ValueError("Supply matching nonempty angular directions and weights.")
    sources = [q_coeff] * N if isinstance(q_coeff, mfem.Coefficient) else list(q_coeff)
    inflows = [inflow] * N if isinstance(inflow, mfem.Coefficient) else list(inflow)
    if len(sources) != N or len(inflows) != N:
        raise ValueError("Source and inflow lists must match the angular order.")
    L = angular.build_morel_angular_operator(mu_vals, w_vals)
    energy_direction = coeff.VectorConstCoefficient(0.0, -1.0)
    velocity_E = mfem.ScalarVectorProductCoefficient(S_coeff, energy_direction)
    ndof = fes.GetVSize()
    size = N * ndof if mass_beta is None else 2 * N * ndof
    rhs = mfem.Vector(size)
    rhs.Assign(0.0)
    transport = []
    for m, mu in enumerate(mu_vals):
        velocity_x = coeff.VectorConstCoefficient(mu, 0.0)
        transport.append(assemble.Bilinear_Form(
            fes, velocity_x, velocity_E, absorption_coeff, dir_bdr1, dir_bdr2))
        source = assemble.Linear_Form(fes, sources[m], inflows[m], velocity_x,
                                      velocity_E, dir_bdr1, dir_bdr2)
        rhs.GetDataArray()[m * ndof:(m + 1) * ndof] = source.GetDataArray()
    if mass_beta is None:
        A = assemble.FP_Matrix(transport, L, mass)
    else:
        A = assemble.GFP2_Matrix(transport, L, mass, mass_beta, mass0)
    diagonal = mfem.Vector(size)
    A.GetDiag(diagonal)
    if np.any(diagonal.GetDataArray() == 0.0):
        raise ValueError("GSSmoother requires nonzero diagonal entries; for GFP2 with "
                         "S=absorption=0, use an even angular order to avoid mu=0.")
    solution = mfem.Vector(size)
    solution.Assign(1.0)
    solution, info = GMRES_solver(A, rhs, solution, iter_, tol, p_level, return_info=True)
    residual = mfem.Vector(size)
    A.Mult(solution, residual)
    residual -= rhs
    info['true_residual_norm'] = residual.Norml2()
    rhs_norm = rhs.Norml2()
    info['relative_true_residual'] = residual.Norml2() / rhs_norm if rhs_norm else None

    psi_mu_list = []
    for m, (mu, w) in enumerate(zip(mu_vals, w_vals)):
        psi = Initial_Guess(fes, 0.0)
        psi.GetDataArray()[:] = solution.GetDataArray()[m * ndof:(m + 1) * ndof]
        if visualize:
            vis.GlVis_2D(mesh, psi)
        if save:
            psi.Save("psi_hs_mu_{:.3f}.gf".format(mu))
        psi_mu_list.append((mu, w, psi))
    if return_info:
        return psi_mu_list, info
    return psi_mu_list
