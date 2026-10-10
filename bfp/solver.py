import warnings
from inspect import signature

import mfem.ser as mfem
import numpy as np
from bfp import angular, assemble, vis, coeff

__all__ = [
    "GMRES_solver",
    "Initial_Guess",
    "Solve_Phi",
    "Solve_Psi",
    "Solve_FP",
    "Solve_GFP2",
    "Solve_Transport",
]


def Solve_Transport(
    model,
    *,
    iter_=1000,
    tol=1e-12,
    p_level=1,
    visualize=True,
    save=True,
    return_info=False,
    **problem,
):
    """Select a transport model without changing its numerical solver.

    Args:
        model (str): ``"csda"``, ``"fp"`` or ``"gfp2"``, selecting
            :func:`Solve_Psi`, :func:`Solve_FP` or :func:`Solve_GFP2`.
        iter_ (int): Maximum GMRES iterations.
        tol (float): Both absolute and relative GMRES tolerance.
        p_level (int): MFEM solver print level; 0 is quiet.
        visualize (bool): Render angular fluxes with PyGLVis.
        save (bool): Write angular flux files.
        return_info (bool): Include the selected solver's diagnostics.
        **problem: The selected solver's problem arguments by keyword: mesh,
            quadrature, coefficients, source and boundary data. No missing
            scattering coefficients are inferred from a CSDA configuration.

    Returns:
        list or tuple: Angular fluxes, or ``(fluxes, info)`` if requested.
        CSDA diagnostics are per angle; FP/GFP2 diagnostics describe the
        coupled system. Direct solver return values are preserved.

    Raises:
        ValueError: The model name is unknown.
        TypeError: Arguments are missing/unexpected or required material
            coefficients are not MFEM scalar coefficients. These checks
            precede assembly; value checks remain in the selected solver.
    """
    if model == "csda":
        solve = Solve_Psi
        coefficients = ("xs_t_coeff",)
    elif model == "fp":
        solve = Solve_FP
        coefficients = ("absorption_coeff", "transport_coeff", "S_coeff")
    elif model == "gfp2":
        solve = Solve_GFP2
        coefficients = ("absorption_coeff", "alpha_coeff", "beta_coeff", "S_coeff")
    else:
        raise ValueError("Unknown transport model; choose 'csda', 'fp' or 'gfp2'.")

    problem.update(
        iter_=iter_,
        tol=tol,
        p_level=p_level,
        visualize=visualize,
        save=save,
        return_info=return_info,
    )
    signature(solve).bind(**problem)
    for name in coefficients:
        if not isinstance(problem[name], mfem.Coefficient):
            raise TypeError(f"{model}: {name} must be an MFEM scalar coefficient.")
    return solve(**problem)


def GMRES_solver(A, b, psi, iter_=1000, tol=1e-12, p_level=1, *, return_info=False):
    """Solve ``A psi = b`` in place with GMRES and GSSmoother.

    ``A`` is an MFEM sparse matrix; ``b`` and ``psi`` are compatible MFEM
    vectors (a GridFunction is also accepted for ``psi``). ``iter_`` limits
    iterations, ``tol`` sets both absolute and relative tolerances, and
    ``p_level`` controls MFEM output (0 is quiet). The restart size is 30.

    Return ``psi``, or ``(psi, info)`` when ``return_info=True``. The dictionary
    contains ``converged`` (bool), ``iterations`` (int) and ``final_norm``
    (MFEM's preconditioned residual norm, not the true residual). Failure to
    converge emits a RuntimeWarning and still returns the current iterate.
    """
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

    info = dict(
        converged=solver.GetConverged(),
        iterations=solver.GetNumIterations(),
        final_norm=solver.GetFinalNorm(),
    )
    if not info["converged"]:
        warnings.warn(
            f"GMRES did not converge after {info['iterations']} iterations "
            f"(reported residual norm: {info['final_norm']:.6g}).",
            RuntimeWarning,
            stacklevel=2,
        )
    if return_info:
        return psi, info
    return psi


def Initial_Guess(fes, psi_const):
    """Return an MFEM GridFunction on ``fes`` with every DOF set to ``psi_const``."""
    psi = mfem.GridFunction(fes)
    psi._fes = fes
    psi.Assign(psi_const)

    return psi


def Solve_Phi(fes, psi_mu_list, *, save=True):
    """Return scalar flux on ``fes`` by summing ``(mu, weight, GridFunction)`` entries.

    The integral is ``sum(weight * psi)`` without a factor of 1/2. With
    ``save=True``, write ``phi_dg.gf`` in the current working directory.
    """
    phi = angular.integrate_angular_flux(fes, psi_mu_list)
    if save:
        phi.Save("phi_dg.gf")

    return phi


def Solve_Psi(
    pn,
    mu_vals,
    w_vals,
    mesh,
    fes,
    xs_t_coeff,
    xs_t_const,
    inflow,
    S_const,
    alpha,
    beta,
    dir_bdr1,
    dir_bdr2,
    a_const,
    b_const,
    q_const=0,
    E_start=0,
    E_end=0,
    iter_=1000,
    tol=1e-12,
    p_level=1,
    *,
    S_coeff=None,
    q_coeff=None,
    visualize=True,
    save=True,
    return_info=False,
):
    """Solve CSDA streaming/removal independently at each angular direction.

    The equation is ``mu*d_x(psi) - d_E(S*psi) + Sigma_t*psi = q`` on the
    coupled physical ``(x, E)`` DG mesh; there is no angular scattering gain.

    Args:
        pn (int): Source/inflow convention: 1–2 for prescribed source,
            3–5 for linear-position, linear-energy or bilinear manufactured
            solutions. Explicit coefficient overrides take precedence.
        mu_vals, w_vals (sequence): Matching directions and quadrature weights.
        mesh: Serial MFEM mesh with position in cm and kinetic energy in MeV.
        fes: Scalar DG finite element space on ``mesh``.
        xs_t_coeff: MFEM removal coefficient in inverse cm.
        xs_t_const (float): Constant removal used to construct manufactured sources.
        inflow: MFEM scalar boundary coefficient, or a number passed to
            :class:`~bfp.coeff.InflowCoefficient` with ``pn``, ``a_const``, ``b_const``.
        S_const (float): Positive stopping power in MeV/cm (0 disables CSDA).
        alpha, beta (float): DG flux parameters, normally 1 and 0.5;
            these are not GFP2 material coefficients.
        dir_bdr1, dir_bdr2: Spatial and energy boundary markers from
            :func:`~bfp.assemble.BoundaryConditions`. Mark both energy faces:
            high-energy inflow and low-energy outflow. Spatial inflow is on
            the left for positive mu and the right for negative mu.
        a_const, b_const (float): Manufactured flux coefficients.
        q_const: Source data for problems 1–2, in angular-flux units per cm.
        E_start, E_end (float): High/low energy bounds for grouped source data.
        iter_, tol, p_level: Iteration limit, tolerance and print level;
            see :func:`GMRES_solver`.
        S_coeff: Optional continuous MFEM stopping-power coefficient overriding
            ``S_const``. The conservative weak form includes the derivative of S.
        q_coeff: Optional MFEM source coefficient overriding the built-in source.
            For variable-S manufactured solutions, supply the full conservative
            source and exact boundary inflow explicitly.
        visualize, save (bool): Render with PyGLVis and write per-angle
            ``psi_hs_mu_*.gf`` files, respectively.
        return_info (bool): Also return convergence diagnostics.

    Returns:
        list or tuple: ``(mu, weight, GridFunction)`` entries, or
        ``(fluxes, info_list)``. Each info dictionary adds ``mu`` to the
        :func:`GMRES_solver` fields. Nonconvergence warns and retains the iterate.
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
        inflow_coeff = (
            inflow
            if isinstance(inflow, mfem.Coefficient)
            else coeff.InflowCoefficient(inflow, mu, pn=pn, a_val=a, b_val=b)
        )
        v_coeff1 = coeff.VectorConstCoefficient(mu, 0.0)
        if source_coeff is not None:
            q_coeff = source_coeff
        elif pn == 1 or pn == 2:
            q_coeff = coeff.QCoefficientE(q_const, E_start, E_end)
        else:
            q_coeff = coeff.QFuncCoefficient(
                pn, a_val=a, b_val=b, xs_t_val=xs_t_const, mu_val=mu, S_val=S_const
            )

        A = assemble.Bilinear_Form(
            fes,
            v_coeff1,
            v_coeff2,
            xs_t_coeff,
            dir_bdr1,
            dir_bdr2,
            alpha=alpha,
            beta=beta,
        )
        b = assemble.Linear_Form(
            fes,
            q_coeff,
            inflow_coeff,
            v_coeff1,
            v_coeff2,
            dir_bdr1,
            dir_bdr2,
            alpha=alpha,
            beta=beta,
        )

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


def Solve_FP(
    mu_vals,
    w_vals,
    mesh,
    fes,
    *,
    absorption_coeff,
    transport_coeff,
    S_coeff,
    q_coeff,
    inflow,
    dir_bdr1,
    dir_bdr2,
    iter_=1000,
    tol=1e-12,
    p_level=1,
    visualize=True,
    save=True,
    return_info=False,
):
    """Solve coupled FP–CSDA with Morel angular differencing.

    The equation is ``mu*d_x(psi) - d_E(S*psi) + Sigma_abs*psi
    - (Sigma_tr/2)*L*psi = q``. Elastic scattering is conservative; do not
    include the total elastic cross section in ``absorption_coeff``.

    Args:
        mu_vals, w_vals (sequence): Increasing full-range Gauss–Legendre
            directions and positive weights summing to 2.
        mesh, fes: Serial MFEM mesh and scalar DG space in physical ``(x, E)``.
        absorption_coeff: Genuine absorption, in inverse cm, as an MFEM coefficient.
        transport_coeff: Elastic transport cross section, in inverse cm;
            finite and nonnegative at mass-integration points.
        S_coeff: MFEM stopping power in MeV/cm; positive for energy loss.
            Constant or continuous variable values are supported; 0 disables CSDA.
        q_coeff: MFEM source in angular-flux units per cm, shared or a list
            of one coefficient per angle in quadrature order.
        inflow: Incoming angular flux as an MFEM coefficient or per-angle list.
            Prescribe the high-energy trace and spatial inflow (left for mu > 0,
            right for mu < 0); low energy is outflow.
        dir_bdr1, dir_bdr2: Spatial and energy boundary markers from
            :func:`~bfp.assemble.BoundaryConditions`; retain both energy faces.
        iter_, tol, p_level: :func:`GMRES_solver` controls.
        visualize, save (bool): Render angular fluxes and write per-angle flux files.
        return_info (bool): Include global convergence diagnostics.

    Returns:
        list or tuple: ``(mu, weight, GridFunction)`` entries, or ``(fluxes, info)``.
        Info contains the :func:`GMRES_solver` fields plus ``true_residual_norm``
        (Euclidean norm of A*psi-b) and ``relative_true_residual`` (divided by
        the RHS norm, or None for zero RHS). Nonconvergence returns the iterate
        with a warning. Material coefficients may depend on position and energy.
    """
    mass_tr = assemble.Mass_Matrix(fes, transport_coeff)
    return _solve_coupled(
        mu_vals,
        w_vals,
        mesh,
        fes,
        absorption_coeff,
        S_coeff,
        q_coeff,
        inflow,
        dir_bdr1,
        dir_bdr2,
        mass_tr,
        None,
        None,
        iter_,
        tol,
        p_level,
        visualize,
        save,
        return_info,
    )


def Solve_GFP2(
    mu_vals,
    w_vals,
    mesh,
    fes,
    *,
    absorption_coeff,
    alpha_coeff,
    beta_coeff,
    S_coeff,
    q_coeff,
    inflow,
    dir_bdr1,
    dir_bdr2,
    iter_=1000,
    tol=1e-12,
    p_level=1,
    visualize=True,
    save=True,
    return_info=False,
):
    """Solve GFP2–CSDA using the angular operator ``alpha*L*(I-beta*L)^(-1)``.

    This operator enters with a minus sign on the left-hand side, replacing
    the FP scattering term. ``alpha_coeff`` is an MFEM coefficient in inverse
    cm, finite and strictly positive at mass-integration points. ``beta_coeff``
    is dimensionless, finite and nonnegative there. They are distinct from
    the DG flux parameters in :func:`Solve_Psi`. Use
    :func:`~bfp.coeff.gfp2_coefficients` to convert elastic relaxation cross sections.

    All other parameters, units, inflow conditions and return values follow
    :func:`Solve_FP`. Coefficients may vary in ``(x, E)``; align material
    interfaces with the mesh. The auxiliary DG variable ``chi`` satisfies
    ``(I-beta*L)*chi=psi`` and needs no independent inflow data. It is not
    returned; the true residual diagnostics include the entire auxiliary system.

    With stopping power and absorption both zero, use even angular order to
    avoid a zero transport diagonal in GSSmoother. Nonconvergence warns and
    returns the current angular flux approximation.
    """
    mass_alpha = assemble.Mass_Matrix(fes, alpha_coeff, positive=True)
    mass_beta = assemble.Mass_Matrix(fes, beta_coeff)
    mass0 = assemble.Mass_Matrix(fes, mfem.ConstantCoefficient(1.0))
    return _solve_coupled(
        mu_vals,
        w_vals,
        mesh,
        fes,
        absorption_coeff,
        S_coeff,
        q_coeff,
        inflow,
        dir_bdr1,
        dir_bdr2,
        mass_alpha,
        mass_beta,
        mass0,
        iter_,
        tol,
        p_level,
        visualize,
        save,
        return_info,
    )


def _solve_coupled(
    mu_vals,
    w_vals,
    mesh,
    fes,
    absorption_coeff,
    S_coeff,
    q_coeff,
    inflow,
    dir_bdr1,
    dir_bdr2,
    mass,
    mass_beta,
    mass0,
    iter_,
    tol,
    p_level,
    visualize,
    save,
    return_info,
):
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
        transport.append(
            assemble.Bilinear_Form(
                fes, velocity_x, velocity_E, absorption_coeff, dir_bdr1, dir_bdr2
            )
        )
        source = assemble.Linear_Form(
            fes, sources[m], inflows[m], velocity_x, velocity_E, dir_bdr1, dir_bdr2
        )
        rhs.GetDataArray()[m * ndof : (m + 1) * ndof] = source.GetDataArray()
    if mass_beta is None:
        A = assemble.FP_Matrix(transport, L, mass)
    else:
        A = assemble.GFP2_Matrix(transport, L, mass, mass_beta, mass0)
    diagonal = mfem.Vector(size)
    A.GetDiag(diagonal)
    if np.any(diagonal.GetDataArray() == 0.0):
        raise ValueError(
            "GSSmoother requires nonzero diagonal entries; for GFP2 with "
            "S=absorption=0, use an even angular order to avoid mu=0."
        )
    solution = mfem.Vector(size)
    solution.Assign(1.0)
    solution, info = GMRES_solver(
        A, rhs, solution, iter_, tol, p_level, return_info=True
    )
    residual = mfem.Vector(size)
    A.Mult(solution, residual)
    residual -= rhs
    info["true_residual_norm"] = residual.Norml2()
    rhs_norm = rhs.Norml2()
    info["relative_true_residual"] = residual.Norml2() / rhs_norm if rhs_norm else None

    psi_mu_list = []
    for m, (mu, w) in enumerate(zip(mu_vals, w_vals)):
        psi = Initial_Guess(fes, 0.0)
        psi.GetDataArray()[:] = solution.GetDataArray()[m * ndof : (m + 1) * ndof]
        if visualize:
            vis.GlVis_2D(mesh, psi)
        if save:
            psi.Save("psi_hs_mu_{:.3f}.gf".format(mu))
        psi_mu_list.append((mu, w, psi))
    if return_info:
        return psi_mu_list, info
    return psi_mu_list
