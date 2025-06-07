from setuptools import setup, find_packages # find_packages might not be strictly needed if packages are simple
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
import os

# Define the C++ extension module
# The name must match what you import in Python, e.g., from microaggregation.microaggregation_cpp import ...
ext_modules = [
    Pybind11Extension(
        "microaggregation.microaggregation_cpp", # This will create microaggregation_cpp.pyd inside the microaggregation package
        sources=[
            "src/microaggregation_cpp.cpp",
            "src/utility.cpp",
            "src/myTiming.cpp"
        ],
        include_dirs=[
            "src",  # For your local headers like myTiming.h, utility.h
            pybind11.get_include(),
            # Add other necessary include paths here if any (e.g., for NumPy if not handled by pybind11)
            # os.path.join(pybind11.get_include(True), '..', 'include'), # Sometimes needed for pybind11 internals
            # np.get_include() # If you were using NumPy C API directly in C++
        ],
        language='c++',
        cxx_std=14, # or 17 if your code uses C++17 features
        # Optional: Add compiler/linker arguments if needed
        # extra_compile_args=['/EHsc', '/bigobj'] if platform.system() == "Windows" else [],
        # extra_link_args=[],
    ),
]

# Read the long description from README for PyPI
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    # name, version, author, etc., are now primarily managed by pyproject.toml
    # However, some can be kept for compatibility or if setup.py is run directly by older tools.
    # For a pyproject.toml-driven build, these might be ignored or cause warnings if duplicated.
    # Let's keep it minimal here, relying on pyproject.toml.
    long_description=long_description, # Often good to provide for setup.py context
    long_description_content_type="text/markdown",
    packages=find_packages(where="."), # Explicitly tell it where to find packages
                                       # e.g. ['microaggregation'] if your microaggregation Python module is in the root
                                       # If microaggregation/microaggregation.py, then find_packages() is fine.
                                       # Your structure is microaggregation/microaggregation/__init__.py
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False, # C-extensions generally make the package not zip-safe
    # install_requires, python_requires, classifiers are better in pyproject.toml
)
