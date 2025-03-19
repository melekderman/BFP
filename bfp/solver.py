import mfem.ser as mfem
from bfp import vis, coeff

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
    psi.Assign(psi_const)
    
    return psi


def Solve_Phi(fes, psi_mu_list):
    phi = mfem.GridFunction(fes)
    phi.Assign(0.0)
    for _, w, X in psi_mu_list:
        phi.Add(w, X)
    phi.Save("phi_dg.gf")

    return phi


def Solve_Psi(pn, mu_vals, w_vals, mesh, fes, xs_t_coeff, xs_t_const, inflow, S_const, alpha, beta, dir_bdr1, dir_bdr2, a_const, b_const, iter_=1000, tol=1e-12, p_level=1):

    psi_mu = []
    psi_mu_list = []

    for mu, w in zip(mu_vals, w_vals):
        a = a_const
        b = b_const
        print("  Solving for mu =", mu)
        inflow_coeff = coeff.InflowCoefficient(inflow, mu)
        v_coeff1 = coeff.VectorConstCoefficient([mu, 0.0])
        v_coeff2 = coeff.VectorConstCoefficient([0.0, S_const])
        q_coeff = coeff.QFuncCoefficient(pn, a_val=a, b_val=b, xs_t_val=xs_t_const, mu_val=mu, S_val=S_const)

        a = mfem.BilinearForm(fes)
        a.AddDomainIntegrator(mfem.ConvectionIntegrator(v_coeff1, 1.0))
        a.AddDomainIntegrator(mfem.ConvectionIntegrator(v_coeff2, 1.0))
        a.AddDomainIntegrator(mfem.MassIntegrator(xs_t_coeff))
        a.AddInteriorFaceIntegrator(mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff1, -alpha, beta)))
        a.AddBdrFaceIntegrator(mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff1, -alpha, beta)), dir_bdr1)
        a.AddInteriorFaceIntegrator(mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff2, -alpha, beta)))
        a.AddBdrFaceIntegrator(mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff2, -alpha, beta)), dir_bdr2)
        a.Assemble()
        a.Finalize()
        A = a.SpMat()

        b = mfem.LinearForm(fes)
        b.AddDomainIntegrator(mfem.DomainLFIntegrator(q_coeff))
        b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff1, -alpha, -beta), dir_bdr1)
        b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff2, -alpha, -beta), dir_bdr2)
        b.Assemble()

        psi = Initial_Guess(fes, 1.0)
        GMRES_solver(A, b, psi, iter_, tol, p_level)
        vis.GlVis_2D(mesh, psi)
        psi.Save("psi_hs_mu_{:.3f}.gf".format(mu))
        psi_mu.append(psi.GetDataArray())
        psi_mu_list.append((mu, w, psi))

    return psi_mu_list

    b = mfem.LinearForm(fes)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(q_coeff))
    b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff1, -alpha, -beta), dir_bdr1)
    b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff2, -alpha, -beta), dir_bdr2)
    b.Assemble()

    return b