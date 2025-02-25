import os
import sys
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

# Check for the "--auto" flag on the command line.
AUTO = False
if "--auto" in sys.argv:
    AUTO = True
    sys.argv.remove("--auto")

BANNER = r"""
 _______        _______      ___________  
|   _  "\      /"     "|    ("     _   ") 
(. |_)  :)    (: ______)     )__/  \\__/  
|:     \/      \/    |          \\_ /     
(|  _  \\      // ___)          |.  |     
|: |_)  :)    (:  (             \:  |     
(_______/      \__/              \__|     
                                          
          BFT Installation
"""

class CustomInstallCommand(install):
    def initialize_options(self):
        install.initialize_options(self)
        # Use the global AUTO flag.
        self.auto = AUTO

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        print(BANNER)
        cwd = os.getcwd()

        if self.auto:
            print("Automatic clone and installation mode selected.")

            # First, install additional dependencies from requirements.txt if it exists.
            req_file = os.path.join(cwd, "requirements.txt")
            if os.path.exists(req_file):
                print("Installing additional dependencies from requirements.txt...")
                try:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", req_file])
                except subprocess.CalledProcessError as e:
                    print("Failed to install requirements:", e)
                    sys.exit(1)
            else:
                print("No requirements.txt found. Skipping additional dependency installation.")

            # --- Clone and install PyMFEM (always in parallel mode) ---
            pymfem_dir = os.path.join(cwd, "PyMFEM")
            if not os.path.exists(pymfem_dir):
                try:
                    print("Cloning PyMFEM repository...")
                    subprocess.check_call(["git", "clone", "https://github.com/melekderman/PyMFEM.git"])
                except subprocess.CalledProcessError as e:
                    print("Failed to clone PyMFEM repository:", e)
                    sys.exit(1)
            else:
                print("PyMFEM repository already exists. Skipping clone.")

            # --- Ensure mpi4py is installed (required by PyMFEM) ---
            try:
                import mpi4py
            except ImportError:
                print("mpi4py is not installed. Installing mpi4py...")
                try:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", "mpi4py"])
                except subprocess.CalledProcessError as e:
                    print("Failed to install mpi4py:", e)
                    sys.exit(1)

            try:
                print("Installing PyMFEM (parallel mode)...")
                os.chdir(pymfem_dir)
                subprocess.check_call([sys.executable, "setup.py", "install", "--with-parallel"])
                os.chdir(cwd)
            except subprocess.CalledProcessError as e:
                print("PyMFEM installation failed:", e)
                sys.exit(1)

            # --- Clone and install pyglvis ---
            pyglvis_dir = os.path.join(cwd, "pyglvis")
            if not os.path.exists(pyglvis_dir):
                try:
                    print("Cloning pyglvis repository...")
                    subprocess.check_call(["git", "clone", "https://github.com/melekderman/pyglvis.git"])
                except subprocess.CalledProcessError as e:
                    print("Failed to clone pyglvis repository:", e)
                    sys.exit(1)
            else:
                print("pyglvis repository already exists. Skipping clone.")

            try:
                print("Installing pyglvis...")
                os.chdir(pyglvis_dir)
                subprocess.check_call([sys.executable, "-m", "pip", "install", "."])
                os.chdir(cwd)
            except subprocess.CalledProcessError as e:
                print("pyglvis installation failed:", e)
                sys.exit(1)
        else:
            print("Automatic clone and installation not selected. Installing only the main module.")

        # Finally, install the main package (e.g. with matplotlib as a dependency)
        install.run(self)

setup(
    name="BFP",
    version="0.0.1",
    packages=find_packages(),
    include_package_data=True,
    description="BFP: A Boltzman Fokker-Plank Transport Solver integrating with PyMFEM",
    long_description=open("README.md", encoding="utf-8").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    install_requires=["requirements.txt"],
    author="Melek Derman",
    author_email="dermanm@oregonstate.edu",
    url="https://github.com/melekderman/BFP",
    cmdclass={'install': CustomInstallCommand},
)
