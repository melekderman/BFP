# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os, sys

sys.path.insert(0, os.path.abspath('../../bfp'))

# On Read the Docs, need to mock any python packages that would require c
from unittest.mock import MagicMock

MOCK_MODULES = [
    "matplotlib.pyplot",
    "mfem",
    "pyglvis",
]

# -- Project information -----------------------------------------------------

project = 'BFP'
project_copyright_info = '2025, Melek Derman'
author = 'Melek Derman'
release = 'v0.0.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = ["sphinx.ext.autodoc", 
              "sphinx.ext.napoleon", 
              "sphinx.ext.mathjax", 
              "sphinx.ext.intersphinx",
              "sphinx.ext.githubpages",
              "sphinx.ext.autosummary",
              "sphinx_toolbox.github",
              "sphinx_toolbox.sidebar_links",
              "sphinx.ext.autosectionlabel",
              'myst_parser', 
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3.11', None)
    # 'PyMFEM': ('https://github.com/mfem/PyMFEM', None),
    # 'pyglvis': ('https://github.com/GLVis/pyglvis', None),
}

autosummary_generate = True
napoleon_numpy_docstring = True
napoleon_google_docstring = True

github_username = "melekderman"
github_repository = "BFP"
github_url = "https://github.com/{github_username}/{github_repository}"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "classic"
#html_logo = "images/home/logo.svg"

# html_permalinks = ['https://cement-psaap.github.io/', 'https://github.com/CEMeNT-PSAAP/MCDC']

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
