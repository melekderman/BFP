# Boltzman Fokker-Plank Transport Solver

#![mcdc_logo v1](https://user-images.githubusercontent.com/26186244/173467190-74d9b09a-ef7d-4f0e-8bdf-4a076de7c43c.svg)

#[![Build](https://github.com/CEMeNT-PSAAP/MCDC/actions/workflows/mpi_numba_reg.yml/badge.svg)](https://github.com/CEMeNT-PSAAP/MCDC/actions/workflows/mpi_numba_reg.yml)
#[![DOI](https://joss.theoj.org/papers/10.21105/joss.06415/status.svg)](https://doi.org/10.21105/joss.06415)
[![ReadTheDocs](https://readthedocs.org/projects/bfp/badge/?version=latest&style=flat)](https://bfp.readthedocs.org/en/latest/ )
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This package offers a solver for Boltzman Fokker-Plank charged-particle transport problems, seamlessly integrating with [PyMFEM](https://github.com/mfem/PyMFEM) and [pyglvis](https://github.com/GLVis/pyglvis). The following instructions provide guidance on installation, running the solver, and testing.

To clone this repository, use the following command:
```
git clone https://github.com/melekderman/BFP.git
```

## Installation

It is recommended to work within isolated environments when installing this package.

### Environment Setup

#### Virtual Environment

Create and activate a virtual environment using Python 3.11:
```
    cd BFP
    python3.11 -m venv .venv
    source .venv/bin/activate
```
#### Conda Environment

If you prefer using Conda, create and activate an environment as follows:
```
    conda create --name BFP-env python==3.11
    conda activate BFP-env
```
### Known issues with pip version

Some versions of pip may produce errors. The version known to work reliably is pip 23.2.1. Please ensure your pip version is compatible before proceeding with the installation.

## Testing

Details on how to run tests and verify the installation will be added here.

```
python test.py -serial
```
```
python test.py -parallel
```

## Acknowledgements

This work was supported by the Center for Exascale Monte-Carlo Neutron Transport (CEMeNT) a PSAAP-III project funded by the Department of Energy, grant number: DE-NA003967.

## License

BFP is licensed under a BSD-3 clause license. We believe in open source software.
