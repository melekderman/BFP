# Boltzmann Fokker-Planck Charged Particle Transport Solver (BFP)

[![ReadTheDocs](https://readthedocs.org/projects/bfp/badge/?version=latest&style=flat)](https://bfp.readthedocs.org/en/latest/ )
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Unit Tests](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/melekderman/BFP/actions/workflows/unit_tests.yml)

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

We recommend using a virtual or Conda environment.

#### Virtual Environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

#### Conda Environment

```bash
conda create --name BFP-env python==3.11
conda activate BFP-env
```
### Step 3: Install BFP

Install BFP along with dependencies using:

```bash
pip install .
```
## 📈 Usage & Visualization

BFP integrates with **PyGLVis**, which is a Jupyter-compatible visualization toolkit. To effectively visualize results:

### In a Jupyter notebook, initialize PyGLVis:

```python
import glvis
glvis(mesh, solution)
```
Refer to the PyGLVis Documentation: [PyGLVis](https://github.com/GLVis/pyglvis) for more examples.

---

## 🔬 Running Tests

After installation, run the provided tests to verify your setup.

#### Serial Testing

```bash
cd tests/unit/
python test_coeff.py
python test_mesh.py
```

#### Parallel Testing *(upcoming feature)*

```bash
python test_coeff.py --parallel
python test_mesh.py --parallel
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