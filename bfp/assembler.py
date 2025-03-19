import mfem.ser as mfem

__all__ = ['Bilinear_Form',
           'Linear_Form',
           ]

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