# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath('../../bfp'))

project = 'BFP'
copyright = '2025, Melek Derman'
author = 'Melek Derman'
release = 'v0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser', 
              'sphinx.ext.autodoc', 
              'sphinx.ext.napoleon', 
              'sphinx.ext.mathjax', 
              'sphinx.ext.intersphinx',
              'sphinx.ext.githubpages' ]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3.11', None),
    # 'PyMFEM': ('https://github.com/mfem/PyMFEM', None),
    # 'pyglvis': ('https://github.com/GLVis/pyglvis', None),
}

napoleon_numpy_docstring = False
napoleon_google_docstring = True

autodoc_default_options = {
    'members': True
}
autoclass_content = 'class'

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
