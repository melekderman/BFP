import mfem.ser as mfem

__all__ = ['BoundaryConditions',
           'FESpace',
           'UniformRefinement',
           'Bilinear_Form',
           'Linear_Form',
           ]


def BoundaryConditions(mesh, x_min=1, x_max=1, y_min=1, y_max=1):

    dir_bdr1 = mfem.intArray(mesh.bdr_attributes.Max())
    dir_bdr1.Assign(0)
    dir_bdr1[3] = x_min
    dir_bdr1[1] = x_max

    dir_bdr2 = mfem.intArray(mesh.bdr_attributes.Max())
    dir_bdr2.Assign(0)
    dir_bdr2[0] = y_max
    dir_bdr2[2] = y_min

    return dir_bdr1, dir_bdr2


def FESpace(order, mesh):
    #mesh.UniformRefinement()
    dim = mesh.Dimension()
    fec = mfem.DG_FECollection(order, dim, btype=0)
    fes = mfem.FiniteElementSpace(mesh, fec)
    return fes


def UniformRefinement(mesh, ref_level=0):
    for _ in range(ref_level):
        mesh.UniformRefinement()
    dim = mesh.Dimension()
    return mesh, dim


def Bilinear_Form(fes, v_coeff1, v_coeff2, xs_t_coeff, dir_bdr1, dir_bdr2, alpha = 1.0, beta = 0.5):
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

    return A


def Linear_Form(fes, q_coeff, inflow_coeff, v_coeff1, v_coeff2, dir_bdr1, dir_bdr2, alpha = 1.0, beta = 0.5):
    b = mfem.LinearForm(fes)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(q_coeff))
    b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff1, -alpha, -beta), dir_bdr1)
    b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow_coeff, v_coeff2, -alpha, -beta), dir_bdr2)
    b.Assemble()

    return b

