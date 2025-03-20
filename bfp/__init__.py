# __init__.py
"""
Problem Package for the CSDA Transport Equation Solver using PyMFEM.

This package contains modules for mesh generation, system assembly,
solvers, and visualization.
"""

__version__ = "0.0.1"

from .input import (problem1_input,
                    problem2_input,
                    problem3_input,
                    problem4_input,
                    problem5_input,
)

from .mesh import (create_2D_mesh,
                   create_3D_mesh
)

from .utils import (read_data,
                    gauss_legendre_dirs,
                    FiniteDiffDerivative,
                    gridfunction_to_array,
                    assign_boundary_attributes,
                    get_marker_for_mu,
                    save_angular_flux,
)

from .coeff import (TotalXSCoefficientE,
                    ScatteringXSCoefficientE,
                    StoppingPowerCoefficientE,
                    QCoefficientE,
                    QFuncCoefficient,
                    EDependentCoefficient,
                    XDependentCoefficient,
                    InflowCoefficient,
                    ConstantCoefficient,
                    VectorConstCoefficient,
)

from .solver import (GMRES_solver,
                     Initial_Guess,
                     Solve_Phi,
                     Solve_Psi
)

from .assemble import (BoundaryConditions,
                       FESpace,
                       UniformRefinement,
                       Bilinear_Form,
                       Linear_Form,
)

from .vis import (GlVis_Visualizer,
                  Matplotlib_Visualizer,
                  Mesh_Report,
                  GlVis_2D,
                  GlVis_3D,
)

__all__ = [
    # input.py
    'problem1_input',
    'problem2_input',
    'problem3_input',
    'problem4_input',
    'problem5_input',

    # mesh.py
    'create_2D_mesh',
    'create_3D_mesh',

    # utils.py
    'read_data',
    'gauss_legendre_dirs',
    'FiniteDiffDerivative',
    'gridfunction_to_array',
    'assign_boundary_attributes',
    'get_marker_for_mu',
    'save_angular_flux',

    # coeff.py
    'TotalXSCoefficientE',
    'ScatteringXSCoefficientE',
    'StoppingPowerCoefficientE',
    'QCoefficientE',
    'QFuncCoefficient',
    'EDependentCoefficient',
    'XDependentCoefficient',
    'InflowCoefficient',
    'ConstantCoefficient',
    'VectorConstCoefficient',

    # vis.py
    'GlVis_Visualizer',
    'Matplotlib_Visualizer',
    'Mesh_Report',
    'GlVis_2D',
    'GlVis_3D',

    # solver.py
    'GMRES_solver',
    'Initial_Guess',
    'Solve_Phi',
    'Solve_Psi',

    # assemble.py
    'BoundaryConditions',
    'FESpace',
    'UniformRefinement',
    'Bilinear_Form',
    'Linear_Form',
    ]
