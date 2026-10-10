Installation
============

The distribution name is ``bfp-transport``; the Python import remains ``bfp``.
The project is not yet published on PyPI. Install from a repository checkout
using Python 3.9, 3.10 or 3.11::

    git clone https://github.com/melekderman/BFP.git
    cd BFP
    python3.11 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip
    export PIP_ONLY_BINARY=mfem
    python -m pip install .

For an editable development install, replace the last command with::

    python -m pip install -e .

Runtime dependencies are PyMFEM ``mfem==4.7.0.1``, ``numpy>=1.20,<2``,
``glvis>=1.1.1``, ``requests``, ``matplotlib`` and ``h5py``. The binary-only
setting prevents an unexpected MFEM source build. Prebuilt wheels support
Linux x86-64 for the supported Python versions; on another platform without
a compatible wheel, install PyMFEM 4.7.0.1 separately first. Package installation
does not clone or compile MFEM. BFP currently imports the serial backend;
installing an MPI-enabled PyMFEM does not make these solvers parallel.

Check an installation and run the existing tests from the checkout::

    python -c "import bfp, mfem.ser, glvis; print(bfp.__version__)"
    bfp --help
    python -m pip install pytest
    python -m pytest tests/unit/

Building this documentation
---------------------------

Only Sphinx and its dependencies are needed; the documentation builder mocks
external numerical packages while importing BFP's actual Python source.
Use a separate Python 3.11 environment from the repository root::

    python3.11 -m venv /tmp/bfp-docs-env
    /tmp/bfp-docs-env/bin/python -m pip install -r docs/requirements.txt
    /tmp/bfp-docs-env/bin/python -m sphinx -W --keep-going -b html docs/source docs/_build/html

Open ``docs/_build/html/index.html`` to browse locally. The existing Makefile
also uses ``_build``; generated output is ignored by Git. The historical files
in ``docs/build/`` are retained but are not the publication source.

GitHub Actions builds on ``dev``/``main`` pushes and pull requests, treating
warnings as errors. Only successful ``main`` pushes upload and deploy a Pages
artifact, with write permissions restricted to the deployment job. A maintainer
must select **Settings → Pages → Build and deployment → Source → GitHub Actions**
to enable this workflow. The ``github-pages`` environment must permit ``main``
deployments. Dev and PR builds never publish the website.
