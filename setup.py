# setup.py
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
import os

# Define the C++ extension module
ext_modules = [
    Pybind11Extension(
        "microaggregation.microaggregation_cpp",
        sources=[
            "src/microaggregation_cpp.cpp",
            "src/utility.cpp",
            "src/myTiming.cpp"
        ],
        include_dirs=[
            "src",
            pybind11.get_include(),
        ],
        language='c++',
        cxx_std=14,
    ),
]

# Read the long description from README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="microaggregation",  # CRUCIAL
    version="0.1.8",          # CRUCIAL - Update this to your current version
    author="Reza Mortazavi",
    author_email="ir1979@gmail.com",
    description="Fast C++ implementation of microaggregation algorithms for data privacy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ir1979/microaggregation", # Update if different
    license="MIT", # Or use license_files if preferred with newer setuptools
    packages=find_packages(), # Should find ['microaggregation']
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    install_requires=[
        "numpy>=1.19.0",
    ],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "License :: OSI Approved :: MIT License", # Still useful for older setuptools
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: C++",
    ],
)
