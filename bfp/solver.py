import mfem.ser as mfem
from bfp import assemble, vis, coeff

__all__ = ['GMRES_solver',
            'Initial_Guess',
            'Solve_Phi',
            'Solve_Psi',
           ]

def GMRES_solver(A, b, psi, iter_=1000, tol=1e-12, p_level=1):
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

    return psi


def Initial_Guess(fes, psi_const):
    psi = mfem.GridFunction(fes)
    psi._fes = fes
    psi.Assign(psi_const)
    
    return psi


def Solve_Phi(fes, psi_mu_list, *, save=True):
    """Integrate angular flux; set save=False to skip the scalar flux file."""
    phi = mfem.GridFunction(fes)
    # PyMFEM 4.7 GridFunction does not retain its finite element space.
    phi._fes = fes
    phi.Assign(0.0)
    for _, w, X in psi_mu_list:
        phi.Add(w, X)
    if save:
        phi.Save("phi_dg.gf")

    return phi


def Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const, inflow, S_const, 
              alpha, beta, dir_bdr1, dir_bdr2, a_const, b_const, q_const=0, E_start=0, E_end=0, iter_=1000, tol=1e-12, p_level=1,
              *, S_coeff=None, q_coeff=None, visualize=True, save=True):
    """Solve in physical (x, E); S_coeff overrides positive S_const.

    S_coeff is a continuous MFEM scalar coefficient. For variable-S MMS,
    supply q_coeff including -d(S*psi)/dE and the exact inflow coefficient.
    Both energy faces must be marked: high-E inflow and low-E outflow.
    Set visualize=False and save=False to disable GLVis and flux files.
    """
    if S_coeff is not None and pn in (3, 4, 5) and q_coeff is None:
        raise ValueError("Variable-S manufactured solutions require q_coeff.")
    source_coeff = q_coeff
    if S_coeff is None:
        S_coeff = mfem.ConstantCoefficient(S_const)
    energy_direction = coeff.VectorConstCoefficient(0.0, -1.0)
    v_coeff2 = mfem.ScalarVectorProductCoefficient(S_coeff, energy_direction)

    psi_mu_list = []

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
        GMRES_solver(A, b, psi, iter_, tol, p_level)
        if visualize:
            vis.GlVis_2D(mesh, psi)
        if save:
            psi.Save("psi_hs_mu_{:.3f}.gf".format(mu))
        psi_mu_list.append((mu, w, psi))

    return psi_mu_list
