Overview
========

BFP solves steady charged-particle transport in a one-dimensional slab, with
physical position :math:`x`, kinetic energy :math:`E`, and direction cosine
:math:`\mu`. This guide explains the model choices and the Python interfaces;
the API reference is generated directly from the source docstrings.

Transport models
----------------

All models retain the coupled spatial–energy discontinuous Galerkin (DG)
discretization and Gauss–Legendre angular quadrature:

.. math::

   \mu\partial_x\psi-\partial_E(S\psi)+\Sigma\psi-\mathcal{A}\psi=Q.

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Model
     - Removal term
     - Angular scattering
   * - ``csda``
     - :math:`\Sigma_t\psi`
     - None; each angular direction is independent.
   * - ``fp``
     - :math:`\Sigma_{\mathrm{abs}}\psi`
     - :math:`\mathcal{A}=(\Sigma_{\mathrm{tr}}/2)L_\mu`
   * - ``gfp2``
     - :math:`\Sigma_{\mathrm{abs}}\psi`
     - :math:`\mathcal{A}=\alpha L_\mu(I-\beta L_\mu)^{-1}`

Here :math:`L_\mu=\partial_\mu[(1-\mu^2)\partial_\mu]` is discretized by
Morel differencing with nonpositive eigenvalues. FP/GFP2 elastic scattering
conserves the weighted angular integral. Absorption is separate from elastic
scattering: do not add the total elastic cross section to the absorption term.
The CSDA option has removal but no angular scattering gain term.

Positive stopping power :math:`S` gives energy velocity :math:`-S`.
The conservative term includes :math:`-(\partial_E S)\psi` when :math:`S`
varies. Prescribe incoming flux at the highest energy, at the left spatial
boundary for :math:`\mu>0`, and at the right boundary for :math:`\mu<0`.
The lowest energy is outflow. Both energy faces remain marked for DG assembly.
Examples use cm, MeV, cross sections in cm\ :sup:`−1`, and stopping power in
MeV/cm; other consistent unit systems are possible.

Choosing a calculation
----------------------

* **Transport model:** :func:`~bfp.solver.Solve_Transport` selects ``csda``,
  ``fp`` or ``gfp2``. The direct solver functions remain available.
* **Problem configuration:** mesh, material coefficients, source and incoming
  flux define the physical problem. The five CSDA examples and the FP/GFP2
  ``slab`` example are separate configurations.
* **Numerical solver:** all current models use GMRES with GSSmoother.
  ``iter_``, ``tol`` and ``p_level`` control this solve, not the transport model.

Current scope
-------------

BFP uses the serial ``mfem.ser`` backend. MPI, GPU transport, higher GFP orders,
and three-dimensional transport are future work. Mesh and visualization helpers
for 3D do not constitute a 3D transport solver. Large or strongly coupled systems
may need improved preconditioning. Discontinuous stopping-power interfaces have
not been verified; use constant or continuous stopping power for verified cases.

The published documentation follows successful builds of ``main``; it is not
a declaration of a new release. ``dev`` and pull requests build documentation
for validation only. The source version is |release|.

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   api
