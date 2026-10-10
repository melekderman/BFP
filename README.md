# Boltzmann Fokker–Planck Charged-Particle Transport (BFP)

[![DOI](https://zenodo.org/badge/938897515.svg)](https://doi.org/10.5281/zenodo.15107100) [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE) [![Unit Tests](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml) [![Coverage](badges/coverage.svg)](badges/coverage.svg) [![Documentation](https://img.shields.io/badge/docs-source-blue)](docs/source/)

## Overview

BFP is a research software tool for developing and evaluating computational methods for strongly forward-peaked charged-particle transport. The project aims to investigate different angular scattering models, numerical discretizations, and solution strategies, with a long-term emphasis on scalable and GPU-accelerated algorithms.

The current implementation uses a coupled spatial–energy discontinuous Galerkin (DG) discretization with discrete-ordinates angular quadrature. It supports CSDA, Fokker–Planck–CSDA, and GFP2–CSDA transport in slab geometry using serial [PyMFEM](https://github.com/mfem/PyMFEM). Parallel and GPU-accelerated solvers remain future development goals.

| Solver | Implemented model |
| --- | --- |
| `bfp.Solve_Psi` | CSDA streaming and removal, solved independently for each angle; no angular scattering coupling. |
| `bfp.Solve_FP` | FP–CSDA with Morel angular differencing and elastic transport scattering. |
| `bfp.Solve_GFP2` | GFP2–CSDA with a rational angular operator and an auxiliary DG variable. |

## Installation

The distribution is named `bfp-transport`, while the Python import is `bfp`. **It has not yet been published on PyPI.** The supported environment is Python 3.9–3.11 with `mfem==4.7.0.1`; other runtime dependencies are declared in [pyproject.toml](pyproject.toml). Prebuilt PyMFEM wheels are available for supported Python versions on Linux x86-64.

```bash
git clone https://github.com/melekderman/BFP.git
cd BFP
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
export PIP_ONLY_BINARY=mfem
python -m pip install .
```

For development, use `python -m pip install -e .` instead. The `PIP_ONLY_BINARY=mfem` setting prevents an unexpected PyMFEM source build; on platforms without a compatible wheel, install a compatible PyMFEM 4.7.0.1 separately.

## Usage

### Command line

The CLI runs five existing CSDA examples; FP/GFP2 are currently available through the Python API rather than CLI models.

```bash
bfp --help
bfp -p 4 --nx 4 --nE 4 --N_ang 4 --iter_ 1000 --tol 1e-12 --p_level 0
```

The problem choices are **1** (infinite medium), **2** (spatial attenuation with a constant source), **3** (linear-in-position manufactured solution), **4** (linear-in-energy manufactured solution), and **5** (bilinear manufactured solution). Run `python main.py` with the same options from the repository checkout. The CLI saves and visualizes results by default; it reports solver nonconvergence but does not currently return a nonzero exit code for that condition.

### Python API

This FP–CSDA example has constant angular flux `ψ = 1` and scalar flux `φ = 2`:

```python
import bfp
import mfem.ser as mfem

mesh = bfp.create_2D_mesh(4, 4, 0.0, 1.0, 0.0, 1.0)
fes = bfp.FESpace(1, mesh)
mu, weights = bfp.gauss_legendre_dirs(4)
dir_bdr1, dir_bdr2 = bfp.BoundaryConditions(mesh)
one = mfem.ConstantCoefficient(1.0)

psi, info = bfp.Solve_FP(
    mu, weights, mesh, fes,
    absorption_coeff=one, transport_coeff=one, S_coeff=one,
    q_coeff=one, inflow=one, dir_bdr1=dir_bdr1, dir_bdr2=dir_bdr2,
    visualize=False, save=False, return_info=True, p_level=0,
)
phi = bfp.Solve_Phi(fes, psi, save=False)
print(info["converged"], info["relative_true_residual"])
```

For `Solve_GFP2`, replace `transport_coeff` with `alpha_coeff` and `beta_coeff` (MFEM coefficients satisfying `alpha > 0` and `beta >= 0`). `bfp.gfp2_coefficients` computes these parameters from the first two elastic angular relaxation cross sections. Sources and inflows may be shared across angles or supplied individually.

All three solvers return `(mu, weight, GridFunction)` tuples for each angle. With `return_info=True`, `Solve_Psi` provides per-angle convergence diagnostics, while FP/GFP2 return global diagnostics including the true coupled residual. `Solve_Phi` integrates the angular flux. Use `visualize=False` and `save=False` to disable solver output; notebook visualization is available through [GLVis](https://github.com/GLVis/pyglvis). See the [example notebooks](examples/) for the original CSDA problems.

## Current limitations

- The current transport solvers use the serial `mfem.ser` backend and one spatial coordinate coupled to energy; MPI, GPU, and 3D transport solvers are not yet implemented.
- Large or strongly coupled systems may need better preconditioning. For GFP2 with both stopping power and absorption equal to zero, an even angular order avoids a zero diagonal in the current smoother.

## Tests and documentation

Run the unit tests from the repository root:

```bash
python -m pip install pytest
python -m pytest tests/unit/
```

The [CI workflow](.github/workflows/unit_tests.yml) checks Python 3.9, 3.10, and 3.11. The current API documentation sources are in [docs/source/](docs/source/); [hosted documentation](https://melekderman.github.io/BFP/) and checked-in generated HTML may describe earlier versions. The [Zenodo DOI](https://doi.org/10.5281/zenodo.15107100) likewise refers to an archived release rather than necessarily to the current development code.

## Contributing and license

Contributions, bug reports, and feature requests are welcome through [GitHub issues](https://github.com/melekderman/BFP/issues) and [pull requests](https://github.com/melekderman/BFP/pulls). BFP is licensed under the [BSD 3-Clause License](LICENSE).

## Acknowledgements

The initial year of BFP development was supported by the **Center for Exascale Monte-Carlo Neutron Transport (CEMeNT)**, a U.S. Department of Energy (DOE) PSAAP-III project (grant DE-NA003967). Subsequent development has been supported by the **Center for Advancing the Radiation Resilience of Electronics (CARRE)**, a DOE PSAAP-IV project (grant DE-NA0004268).
## Contact

Melek Derman — dermanm@oregonstate.edu
