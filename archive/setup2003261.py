# setup.py

"""
ESTIF-Gravity v2.0 - Package Installation
Strong-field modifications to General Relativity with testable predictions
"""

from setuptools import setup, find_packages
import os

# Read the README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Core dependencies (minimal set)
core_requirements = [
    "numpy>=1.24.0,<2.0.0",
    "scipy>=1.10.0,<2.0.0",
    "matplotlib>=3.7.0,<4.0.0",
    "astropy>=5.2.0,<7.0.0",
]

setup(
    name="estif-gravity",
    version="2.0.0",
    author="Peter Angelov",
    author_email="tervion@gmail.com",
    description="Strong-field modifications to General Relativity with testable predictions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tervion/estif-publication",
    project_urls={
        "Bug Tracker": "https://github.com/tervion/estif-publication/issues",
        "Documentation": "https://github.com/tervion/estif-publication/tree/main/docs",
        "Source Code": "https://github.com/tervion/estif-publication",
        "Previous Version (ESTIF-FD)": "https://zenodo.org/records/17261725",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
    python_requires=">=3.8",
    install_requires=core_requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
        ],
        "analysis": [
            "corner>=2.2.0",
            "emcee>=3.1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "estif-gravity=estif_ec_gr_run_simulation:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["data/*.txt"],
    },
    zip_safe=False,
    keywords=[
        "general-relativity",
        "modified-gravity",
        "gravitational-waves",
        "LISA",
        "strong-field",
        "black-holes",
        "event-horizon-telescope",
        "testable-predictions",
    ],
)

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2

