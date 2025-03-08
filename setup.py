import os
import sys
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.install import install

# This is an option for the future, BFT does not work in parallel for now.
# Parallel modes will be added soon...
  
PARALLEL = False
if "--parallel" in sys.argv:
    PARALLEL = True
    sys.argv.remove("--parallel")

BANNER = r"""
 _______        _______      ___________  
|   _  "\      /"     "|    ("     _   ") 
(. |_)  :)    (: ______)     )__/  \\__/  
|:     \/      \/    |          \\_ /     
(|  _  \\      // ___)          |.  |     
|: |_)  :)    (:  (             \:  |     
(_______/      \__/              \__|     
                                          
          BFP Installation
"""

class CustomInstallCommand(install):
    def run(self):
        print(BANNER)
        cwd = os.getcwd()

        dependencies = [
            "numpy>=1.20.0,<2.0.0",
            "numba",
            "numba-scipy",
            "scipy",
            "swig>=4.2.1",
            "cmake",
            "anywidget>=0.9.9",
            "ipywidgets>=8.0.0",
            "requests"
        ]

        print("Installing dependencies...")
        for dep in dependencies:
            subprocess.check_call([sys.executable, "-m", "pip", "install", dep])

        # PyMFEM
        pymfem_dir = os.path.join(cwd, "PyMFEM")
        if os.path.exists(pymfem_dir):
            print("Removing existing PyMFEM directory...")
            shutil.rmtree(pymfem_dir)

        print("Cloning PyMFEM repository...")
        subprocess.check_call(["git", "clone", "https://github.com/melekderman/PyMFEM.git", pymfem_dir])

        os.chdir(pymfem_dir)
        try:
            if PARALLEL:
                print("Installing PyMFEM with parallel support...")
                subprocess.check_call([sys.executable, "setup.py", "install", "--with-parallel"])
            else:
                print("Installing PyMFEM (standard mode)...")
                subprocess.check_call([sys.executable, "setup.py", "install"])
        except subprocess.CalledProcessError as e:
            print("PyMFEM installation failed:", e)
            sys.exit(1)
        finally:
            os.chdir(cwd)
            print("Cleaning up PyMFEM directory...")
            shutil.rmtree(pymfem_dir)

        # PyGLVis
        pyglvis_dir = os.path.join(cwd, "pyglvis")
        if os.path.exists(pyglvis_dir):
            print("Removing existing pyglvis directory...")
            shutil.rmtree(pyglvis_dir)

        print("Cloning pyglvis repository...")
        subprocess.check_call(["git", "clone", "https://github.com/melekderman/pyglvis.git", pyglvis_dir])

        os.chdir(pyglvis_dir)
        try:
            print("Installing pyglvis...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", "."])
        except subprocess.CalledProcessError as e:
            print("pyglvis installation failed:", e)
            sys.exit(1)
        finally:
            os.chdir(cwd)
            print("Cleaning up pyglvis directory...")
            shutil.rmtree(pyglvis_dir)

        # Finally, install the main package
        install.run(self)

setup(
    name="BFP",
    version="0.0.1",
    packages=find_packages(),
    include_package_data=True,
    description="BFP: A Boltzmann Fokker-Planck Charged Particle Transport Solver integrating with PyMFEM",
    long_description=open("README.md", encoding="utf-8").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    install_requires=[
    "numpy>=1.20.0,<2.0.0",
    "numba",
    "numba-scipy",
    "scipy",
    "swig>=4.2.1",
    "cmake",
    "anywidget>=0.9.9",
    "ipywidgets>=8.0.0",
    "requests"
],
    author="Melek Derman",
    author_email="dermanm@oregonstate.edu",
    url="https://github.com/melekderman/BFP",
    cmdclass={'install': CustomInstallCommand},
)