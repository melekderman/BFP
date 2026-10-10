API Reference
=============

Signatures and descriptions below come from BFP's Python source. Functions
exported by ``bfp`` can also be imported from their defining modules.
Material coefficients are serial MFEM scalar coefficients unless stated
otherwise. Physical units must be consistent; examples use cm and MeV.

Transport and numerical solvers
-------------------------------

.. automodule:: bfp.solver
   :members:

Angular discretization
----------------------

.. automodule:: bfp.angular
   :members:

DG assembly and finite element spaces
-------------------------------------

.. automodule:: bfp.assemble
   :members:

Coefficients
------------

.. automodule:: bfp.coeff
   :members:

Mesh helpers
------------

.. toctree::
   :maxdepth: 1

   mesh

Visualization
-------------

.. automodule:: bfp.vis
   :members:

Example configuration
---------------------

.. autofunction:: bfp.input.slab_input
