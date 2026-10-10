import mfem.ser as mfem
import numpy as np

__all__ = [
    "BoundaryConditions",
    "FESpace",
    "UniformRefinement",
    "Bilinear_Form",
    "Linear_Form",
    "Mass_Matrix",
    "FP_Matrix",
    "GFP2_Matrix",
]


def BoundaryConditions(mesh, x_min=1, x_max=1, y_min=1, y_max=1):

    dir_bdr1 = mfem.intArray(mesh.bdr_attributes.Max())
    dir_bdr1.Assign(0)
    dir_bdr1[3] = x_min
    dir_bdr1[1] = x_max

    dir_bdr2 = mfem.intArray(mesh.bdr_attributes.Max())
    dir_bdr2.Assign(0)
    dir_bdr2[0] = y_min
    dir_bdr2[2] = y_max

    return dir_bdr1, dir_bdr2


def FESpace(order, mesh):
    # mesh.UniformRefinement()
    dim = mesh.Dimension()
    fec = mfem.DG_FECollection(order, dim, btype=0)
    fes = mfem.FiniteElementSpace(mesh, fec)
    return fes


def UniformRefinement(mesh, ref_level=0):
    for _ in range(ref_level):
        mesh.UniformRefinement()
    dim = mesh.Dimension()
    return mesh, dim


def Bilinear_Form(
    fes, v_coeff1, v_coeff2, xs_t_coeff, dir_bdr1, dir_bdr2, alpha=1.0, beta=0.5
):
    """Assemble transport with energy velocity v_coeff2 = (0, -S(E))."""
    a = mfem.BilinearForm(fes)
    a.AddDomainIntegrator(mfem.ConvectionIntegrator(v_coeff1, 1.0))
    # Conservative energy weak form includes dS/dE without differentiating S.
    a.AddDomainIntegrator(mfem.ConservativeConvectionIntegrator(v_coeff2, 1.0))
    a.AddDomainIntegrator(mfem.MassIntegrator(xs_t_coeff))
    a.AddInteriorFaceIntegrator(
        mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff1, -alpha, beta))
    )
    a.AddBdrFaceIntegrator(
        mfem.TransposeIntegrator(mfem.DGTraceIntegrator(v_coeff1, -alpha, beta)),
        dir_bdr1,
    )
    a.AddInteriorFaceIntegrator(mfem.DGTraceIntegrator(v_coeff2, alpha, beta))
    a.AddBdrFaceIntegrator(mfem.DGTraceIntegrator(v_coeff2, alpha, beta), dir_bdr2)
    a.Assemble()
    a.Finalize()
    A = a.SpMat()
    A._bilinear_form = a

    return A


def Linear_Form(
    fes,
    q_coeff,
    inflow_coeff,
    v_coeff1,
    v_coeff2,
    dir_bdr1,
    dir_bdr2,
    alpha=1.0,
    beta=0.5,
):
    b = mfem.LinearForm(fes)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(q_coeff))
    b.AddBdrFaceIntegrator(
        mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff1, -alpha, -beta), dir_bdr1
    )
    b.AddBdrFaceIntegrator(
        mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff2, -alpha, -beta), dir_bdr2
    )
    b.Assemble()

    return b


def Mass_Matrix(fes, coefficient, *, positive=False):
    """Assemble a nonnegative mass coefficient, checked at integration points."""
    for i in range(fes.GetNE()):
        fe = fes.GetFE(i)
        transformation = fes.GetElementTransformation(i)
        rule = mfem.MassIntegrator.GetRule(fe, fe, transformation)
        for j in range(rule.GetNPoints()):
            point = rule.IntPoint(j)
            transformation.SetIntPoint(point)
            value = coefficient.Eval(transformation, point)
            if not np.isfinite(value) or value < 0.0 or (positive and value == 0.0):
                bound = "positive" if positive else "nonnegative"
                raise ValueError(
                    f"Mass coefficient must be finite and {bound} at integration points."
                )
    form = mfem.BilinearForm(fes)
    form.AddDomainIntegrator(mfem.MassIntegrator(coefficient))
    form.Assemble()
    form.Finalize()
    matrix = form.SpMat()
    matrix._bilinear_form = form
    return matrix


def _scaled_matrix(matrix, factor):
    result = mfem.SparseMatrix(matrix)
    result *= float(factor)
    return result


def FP_Matrix(transport, L, mass_tr):
    """Angle-major blocks T_m delta_mn - L_mn M_tr/2; L already includes weights."""
    N = len(transport)
    offsets = mfem.intArray([m * mass_tr.Height() for m in range(N + 1)])
    blocks = mfem.BlockMatrix(offsets)
    for m in range(N):
        blocks.SetBlock(m, m, mfem.Add(1.0, transport[m], -0.5 * L[m, m], mass_tr))
        for n in (m - 1, m + 1):
            if 0 <= n < N:
                blocks.SetBlock(m, n, _scaled_matrix(mass_tr, -0.5 * L[m, n]))
    matrix = blocks.CreateMonolithic()
    # MFEM allocates a new matrix; PyMFEM 4.7 does not mark the result as owned.
    matrix.thisown = True
    return matrix


def GFP2_Matrix(transport, L, mass_alpha, mass_beta, mass0):
    """Sparse auxiliary system, ordered (psi_0,...,psi_N-1,chi_0,...,chi_N-1)."""
    N = len(transport)
    offsets = mfem.intArray([m * mass0.Height() for m in range(2 * N + 1)])
    blocks = mfem.BlockMatrix(offsets)
    for m in range(N):
        blocks.SetBlock(m, m, transport[m])
        blocks.SetBlock(N + m, m, _scaled_matrix(mass0, -1.0))
        blocks.SetBlock(N + m, N + m, mfem.Add(1.0, mass0, -L[m, m], mass_beta))
        for n in range(max(0, m - 1), min(N, m + 2)):
            blocks.SetBlock(m, N + n, _scaled_matrix(mass_alpha, -L[m, n]))
            if n != m:
                blocks.SetBlock(N + m, N + n, _scaled_matrix(mass_beta, -L[m, n]))
    matrix = blocks.CreateMonolithic()
    matrix.thisown = True
    return matrix
