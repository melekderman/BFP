# Boltzmann Fokker-Planck Charged Particle Transport Solver (BFP)

[![DOI](https://zenodo.org/badge/938897515.svg)](https://doi.org/10.5281/zenodo.15107100)[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![Unit Tests](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml)[![ReadTheDocs](https://readthedocs.org/projects/bfp/badge/?version=latest&style=flat)](https://bfp.readthedocs.org/en/latest/ )

## 👀 Overview

BFP provides a solver for Boltzman Fokker-Plank charged-particle transport problems, seamlessly integrating with: 
[PyMFEM](https://github.com/mfem/PyMFEM): Python wrapper for MFEM library, supporting FEM simulations.
[PyGLVis](https://github.com/GLVis/pyglvis). Interactive visualization tool for finite element methods, designed especially for Jupyter notebooks.

This package simplifies setup and visualization of BFP simulations, designed specifically for Charged Particle Transport Problems.

### Step 1: Clone the Repository

```bash
git clone https://github.com/melekderman/BFP.git
cd BFP
```

### Step 2: Set Up Your Python Environment

We recommend using an environment.

#### Conda Environment

```bash
conda create -n BFP-env python==3.11
conda activate BFP-env
```
### Step 3: Install BFP

Install BFP along with dependencies using:

```bash
pip install .
```

---

## 📈 Visualization

In a Jupyter notebook:

```python
GlVis_2D(mesh, solution)
```
or 
```python
GlVis_2D(mesh, solution)
```

---

## 🔬 Running Tests

After installation, run the provided tests to verify your setup.

#### Testing in Serial Mode

```bash
cd tests/unit/
python test_coeff.py
python test_mesh.py
python test_prob.py
```

#### Testing in Parallel Mode *(upcoming feature)*

```bash
python test_coeff.py --parallel
python test_mesh.py --parallel
```

---

## ⚛️ Examples

```bash
Run a selected problem input function with customizable parameters.

Available problems:
1: Infinite Medium: ψ = Q/σₜ
2: Exponential Attenuation: ψ(x) = ψₗ * exp(-σₜ * x / μ)
3: MMS - Linear in x: ψ = a + b * x
4: MMS - Linear in E: ψ = a + b * E
5: Mixed: ψ = a + b * xE

The following parameters can be provided:
nx (int): Number of cells in the x-direction (default: 10).
nE (int): Number of energy cells (default: 10).
N_ang (int): Number of angles for the SN method (default: 2).
iter_ (int): Maximum number of solver iterations (default: 1000).
tol (float): Solver tolerance (default: 1e-12).
p_level (int): Print level (1 for verbose, 0 for silent; default: 1).
``` 

```bash
Example usage:
python main.py -p 4 --nx 15 --nE 12
```

---

## 📚 Documentation
For more information about the modules, please visit the [documentation](https://melekderman.github.io/BFP/).
(The documentation page will be updated soon to include all functions.)

---

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue: https://github.com/melekderman/BFP/issues or submit a pull request: https://github.com/melekderman/BFP/pulls.

---

## 📜 License

BFP is released under the BSD-3 Clause License. See LICENSE file for details.

---

## 💬 Acknowledgements

This project is supported by the **Center for Exascale Monte-Carlo Neutron Transport (CEMeNT)**, a PSAAP-III project funded by the Department of Energy (DOE), grant number: DE-NA003967.

---

📮 **Contact**

For any questions or further details, please contact:

📧 Melek Derman – dermanm@oregonstate.edu