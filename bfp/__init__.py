"""
Problem Package for the CSDA Transport Equation Solver using PyMFEM.

This package contains modules for mesh generation, system assembly,
solvers, and visualization.
"""

__version__ = "0.0.1"

from .mesh import create_2D_mesh, create_3D_mesh
from .utils import read_data, gauss_legendre_dirs, compute_S_derivative
from .coeff import (
    TotalXSCoefficient,
    ScatteringXSCoefficient,
    StoppingPowerCoefficient,
    StoppingPowerDerivativeCoefficient,
    InflowCoefficient,
    QCoefficient,
    EnergyDependentCoefficient,
    VelocityCoefficient
)
from .vis import glvis_visualize, matplotlib_visualize

__all__ = [
    # mesh.py
    'create_2D_mesh',
    'create_3D_mesh',

    # utils.py
    'read_data',
    'gauss_legendre_dirs',
    'compute_S_derivative',

    # coeff.py
    'TotalXSCoefficient',
    'ScatteringXSCoefficient',
    'StoppingPowerCoefficient',
    'StoppingPowerDerivativeCoefficient',
    'InflowCoefficient',
    'QCoefficient',
    'EnergyDependentCoefficient',
    'VelocityCoefficient',

    # vis.py
    'glvis_visualize',
    'matplotlib_visualize',
]
