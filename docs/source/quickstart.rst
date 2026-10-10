Quick Start
===========

Command line
------------

The installed ``bfp`` command and ``python main.py`` from the checkout accept
the same arguments::

    bfp --help
    bfp -p 4 --nx 4 --nE 4 --N_ang 4 --iter_ 1000 --tol 1e-12 --p_level 0
    bfp --model fp --problem slab --nx 4 --nE 4 --N_ang 4 --p_level 0
    bfp --model gfp2 --problem slab --nx 4 --nE 4 --N_ang 4 --p_level 0

``--model csda`` is the default and supports problems 1–5: infinite medium,
spatial attenuation, and manufactured solutions linear in position, linear
in energy, and bilinear in position/energy. FP/GFP2 support only ``slab``.
Incompatible model/problem combinations are rejected.

The fixed slab has :math:`x\in[0,1]` cm, :math:`E\in[0,1]` MeV, first-order DG,
stopping power 1 MeV/cm, absorption 1 cm\ :sup:`−1`, and isotropic source 1.
Spatial inflow and high-energy inflow are zero; the low-energy boundary is
outflow. FP uses :math:`\Sigma_{\mathrm{tr}}=1` cm\ :sup:`−1`.
GFP2 uses :math:`\alpha=0.6` cm\ :sup:`−1` and :math:`\beta=0.1`, giving the
same first angular relaxation cross section. This is a simple example, not
the Leakeas–Larsen two-region benchmark.

``nx`` and ``nE`` count spatial and energy elements; ``N_ang`` counts angular
ordinates. ``p_level`` is the GMRES print level, not the DG polynomial order.
The CLI saves flux files and invokes visualization. It reports nonconvergence
but currently retains a zero exit status in that case. Use the Python output
controls below for headless solves.

Python interface
----------------

This FP example uses unit source, absorption and **unit inflow**, so the exact
solution is constant angular flux 1 and scalar flux 2. Its inflow differs from
the zero-inflow CLI slab.

.. code-block:: python

    import bfp
    import mfem.ser as mfem

    mesh = bfp.create_2D_mesh(4, 4, 0.0, 1.0, 0.0, 1.0)
    fes = bfp.FESpace(1, mesh)
    mu, weights = bfp.gauss_legendre_dirs(4)
    spatial, energy = bfp.BoundaryConditions(mesh)
    one = mfem.ConstantCoefficient(1.0)

    psi, info = bfp.Solve_Transport(
        "fp", mu_vals=mu, w_vals=weights, mesh=mesh, fes=fes,
        absorption_coeff=one, transport_coeff=one, S_coeff=one,
        q_coeff=one, inflow=one, dir_bdr1=spatial, dir_bdr2=energy,
        visualize=False, save=False, return_info=True, p_level=0,
    )
    phi = bfp.Solve_Phi(fes, psi, save=False)
    print(info["converged"], info["relative_true_residual"])

For GFP2, select ``"gfp2"`` and replace ``transport_coeff=one`` with
``alpha_coeff=mfem.ConstantCoefficient(0.6)`` and
``beta_coeff=mfem.ConstantCoefficient(0.1)``. Custom coefficients can depend
on ``(x, E)``; source and inflow may be one MFEM scalar coefficient or a list
of coefficients in quadrature order. Keep absorption separate from elastic
scattering. See :func:`~bfp.coeff.gfp2_coefficients` for converting elastic
relaxation cross sections.

For CSDA, select ``"csda"`` and supply the arguments of
:func:`~bfp.solver.Solve_Psi`. The dispatcher does not invent missing material
coefficients or convert an existing CSDA example to FP/GFP2.

Outputs and convergence
-----------------------

The default result is a list of ``(mu, weight, GridFunction)`` entries.
``return_info=True`` returns ``(fluxes, info)``; CSDA has one diagnostics
dictionary per angle, while FP/GFP2 have one dictionary for the coupled solve.
``converged``, ``iterations`` and ``final_norm`` are MFEM's GMRES statistics;
``final_norm`` is the preconditioned residual norm. FP/GFP2 also report
``true_residual_norm`` and ``relative_true_residual`` for the assembled system
(including GFP2's auxiliary unknown). The latter is ``None`` for zero RHS.
Nonconvergence emits a warning and returns the current approximation.

``visualize=False`` and ``save=False`` disable solver visualization and flux
files without changing the calculation. The mesh creation helper still writes
a mesh file under ``mesh/usr/``. PyGLVis rendering requires a compatible
notebook frontend; the socket helper requires an external GLVis server.
