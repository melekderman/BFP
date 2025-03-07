"""
Problem Package for the CSDA Transport Equation Solver using PyMFEM.

This package contains modules for mesh generation, system assembly,
solvers, and visualization.
"""
__version__ = "0.0.1"
import mfem.ser as mfem
import numpy as np
import matplotlib.pyplot as plt
import os
import glvis

from .mesh import create_2D_mesh, create_3D_mesh