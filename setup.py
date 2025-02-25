import os
import sys
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

# ASCII banner to be printed during installation
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
    """
    Custom installation command that:
    1. Prints an ASCII banner "BFT" at the start.
    2. Clones the PyMFEM repository and installs it in both serial and parallel versions.
    3. Clones and installs pyglvis.
    4. Installs additional dependencies listed in requirements.txt.
    """
    def run(self):
        # Print the ASCII banner
        print(BANNER)
        
        cwd = os.getcwd()
        
        # 1. Clone and install PyMFEM
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
        
        try:
            # Install serial version
            print("Installing PyMFEM (serial version)...")
            os.chdir(pymfem_dir)
            subprocess.check_call([sys.executable, "setup.py", "install"])
            
            # Install parallel version
            print("Installing PyMFEM (parallel version)...")
            subprocess.check_call([sys.executable, "setup.py", "install", "--with-parallel"])
            os.chdir(cwd)
        except subprocess.CalledProcessError as e:
            print("PyMFEM installation failed:", e)
            sys.exit(1)
        
        # 2. Clone and install pyglvis
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
        
        # 3. Install additional dependencies from requirements.txt
        req_file = os.path.join(cwd, "requirements.txt")
        if os.path.exists(req_file):
            try:
                print("Installing additional dependencies from requirements.txt...")
                with open(req_file, "r", encoding="utf-8") as f:
                    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
                for dep in requirements:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", dep])
            except subprocess.CalledProcessError as e:
                print("Failed to install a dependency from requirements.txt:", e)
                sys.exit(1)
        else:
            print("No requirements.txt file found. Skipping additional dependencies.")
        
        # Continue with the standard installation process.
        install.run(self)

setup(
    name="BFP",
    version="0.0.1",
    packages=find_packages(),
    include_package_data=True,
    description="BFP: A repository integrating with PyMFEM (both serial and parallel builds)",
    long_description=open("README.md", encoding="utf-8").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    cmdclass={'install': CustomInstallCommand},
)
