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
    """Return spatial and energy MFEM boundary-marker arrays for a rectangle.

    Use the standard quadrilateral attributes: left=4, right=2, low-E=1,
    high-E=3. Arguments set the corresponding marker values (1 enables).
    In transport, y is physical kinetic energy. Keep both energy faces marked:
    the DG operator needs low-E outflow as well as high-E inflow. Flux integrators
    select incoming traces; these arrays do not prescribe the flux values.
    """

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
    """Return a scalar MFEM DG finite element space of polynomial ``order`` on ``mesh``."""
    # mesh.UniformRefinement()
    dim = mesh.Dimension()
    fec = mfem.DG_FECollection(order, dim, btype=0)
    fes = mfem.FiniteElementSpace(mesh, fec)
    return fes


def UniformRefinement(mesh, ref_level=0):
    """Refine ``mesh`` in place ``ref_level`` times; return ``(mesh, dimension)``."""
    for _ in range(ref_level):
        mesh.UniformRefinement()
    dim = mesh.Dimension()
    return mesh, dim


def Bilinear_Form(
    fes, v_coeff1, v_coeff2, xs_t_coeff, dir_bdr1, dir_bdr2, alpha=1.0, beta=0.5
):
    """Return the assembled DG transport/removal SparseMatrix on ``fes``.

    Use vector coefficients ``v_coeff1=(mu, 0)`` and ``v_coeff2=(0, -S)``
    in physical ``(x, E)``. Energy convection is conservative, including the
    derivative of variable stopping power without explicit differentiation.
    ``xs_t_coeff`` is the scalar removal coefficient (genuine absorption in
    FP/GFP2). Boundary markers come from :func:`BoundaryConditions`.
    ``alpha=1`` and ``beta=0.5`` set the existing upwind DG flux; they are
    unrelated to GFP2 material parameters. The matrix retains its MFEM form.
    """
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
    """Return the assembled MFEM LinearForm for source and incoming boundary flux.

    ``q_coeff`` and ``inflow_coeff`` are scalar MFEM coefficients. Velocities,
    markers and DG ``alpha``/``beta`` must match :func:`Bilinear_Form`.
    With energy velocity ``(0, -S)`` and S > 0, prescribed energy inflow enters
    at the highest energy; no incoming data are imposed at low-energy outflow.
    """
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
    """Return the weighted DG mass SparseMatrix for an MFEM scalar ``coefficient``.

    Values must be finite and nonnegative at mass-integration points, or
    strictly positive when ``positive=True``; otherwise raise ValueError.
    The returned matrix retains its MFEM BilinearForm. These sampled checks
    do not establish positivity between quadrature points.
    """
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
    """Return the monolithic FP matrix with angle-major ordering.

    ``transport`` holds one DG SparseMatrix per angle; ``L`` is the dense
    Morel matrix and ``mass_tr`` is weighted by the elastic transport cross
    section. Blocks are ``T_m*delta_mn - L_mn*M_tr/2``. The angular weights
    are already in L; no additional weights enter the coupling.
    """
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
    """Return the monolithic GFP2 auxiliary system in angle-major blocks.

    Ordering is ``(psi_0,...,psi_N-1,chi_0,...,chi_N-1)``. ``transport`` holds
    the DG streaming/CSDA/absorption matrices and ``L`` is the Morel matrix.
    ``mass_alpha``, ``mass_beta`` and ``mass0`` are DG mass matrices weighted
    by GFP2 alpha, beta and unity. Block equations are
    ``T_m*psi_m - sum_n(L_mn*M_alpha*chi_n) = q_m`` and
    ``-M0*psi_m + M0*chi_m - sum_n(L_mn*M_beta*chi_n) = 0``.
    No angular quadrature weights are added to L.
    """
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
