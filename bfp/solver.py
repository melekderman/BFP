import mfem.ser as mfem

__all__ = ['GMRES_solver',
           ]

def GMRES_solver(A, b, psi):
    prec = mfem.GSSmoother(A)
    solver = mfem.GMRESSolver()
    solver.SetOperator(A)
    solver.SetPreconditioner(prec)
    solver.SetRelTol(1e-12)
    solver.SetAbsTol(1e-12)
    solver.SetMaxIter(1000)
    solver.SetKDim(30)
    solver.SetPrintLevel(1)
    solver.Mult(b, psi)

    return psi