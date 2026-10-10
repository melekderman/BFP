BFP documentation
=================

BFP provides CSDA, Fokker–Planck–CSDA and GFP2–CSDA transport in 1D slab
geometry. It uses a coupled spatial–energy DG discretization, Gauss–Legendre
angular quadrature, and the serial PyMFEM backend. MPI and GPU transport
solvers are not implemented.

The distribution is named ``bfp-transport`` and the import is ``bfp``.
Install from the repository with ``python -m pip install .`` or
``python -m pip install -e .`` using Python 3.9–3.11 and PyMFEM 4.7.0.1.
See the `README <https://github.com/melekderman/BFP/blob/main/README.md>`_
for dependencies, binary-wheel requirements, and runnable examples.

The ``bfp`` CLI runs the five existing CSDA examples. FP/GFP2 are selected
through the Python functions below.

Solver interfaces
-----------------

.. autofunction:: bfp.Solve_Psi

.. autofunction:: bfp.Solve_FP

.. autofunction:: bfp.Solve_GFP2

.. autofunction:: bfp.Solve_Phi

.. autofunction:: bfp.gfp2_coefficients

Mesh helpers
------------

.. toctree::
   :maxdepth: 2

   mesh
