# Boltzman Fokker-Plank Transport Solver

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

## Links
- PyMFEM ([github](https://github.com/mfem/pymfem), [pypi](https://pypi.org/project/mfem/))
- pyglvis ([github](https://github.com/GLVis/pyglvis), [pypi]())
