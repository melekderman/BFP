import mfem.ser as mfem

__all__ = ['GMRES_solver',
           ]

def GMRES_solver(A, b, psi, iter=1000, tol=1e-12, p_level=1):
    prec = mfem.GSSmoother(A)
    solver = mfem.GMRESSolver()
    solver.SetOperator(A)
    solver.SetPreconditioner(prec)
    solver.SetRelTol(tol)
    solver.SetAbsTol(tol)
    solver.SetMaxIter(iter)
    solver.SetKDim(30)
    solver.SetPrintLevel(p_level)
    solver.Mult(b, psi)

    return psi