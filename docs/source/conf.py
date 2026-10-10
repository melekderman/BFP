"""Sphinx configuration for the source-based API reference."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

project = "BFP"
copyright = "2026, Melek Derman"
author = "Melek Derman"
release = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
]

# Import BFP's real source without installing the numerical backend.
autodoc_mock_imports = ["mfem", "numpy", "matplotlib", "h5py", "glvis"]
autodoc_member_order = "bysource"
napoleon_numpy_docstring = True
napoleon_google_docstring = True

html_theme = "classic"
html_title = "BFP documentation"
html_baseurl = "https://melekderman.github.io/BFP/"
