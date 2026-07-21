# setup.py
"""
ESTIF v6.1 — Emergent Spacetime from Inward Flow
Geometric derivation of gravity, dark matter, and the MOND acceleration constant.

Key results (v6.1):
  - a₀ = H₀cx₀/√3 derived from geometry, zero free parameters, 1.72% from MOND empirical
  - 87 SPARC galaxies: RMS 15.6%, within observed baryonic Tully-Fisher scatter
  - B = L/3 multiplier derived from 3D spatial isotropy
  - DESI DR2 constraint: Ω_tilt(z) fails at chi²/N=10.8 — cosmology sector under revision
  - Gravity sector (EHT, Planck Λ, LISA): all pass, 0 free parameters
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

core_requirements = [
    "numpy>=1.24.0,<2.0.0",
    "scipy>=1.10.0,<2.0.0",
    "matplotlib>=3.7.0,<4.0.0",
    "astropy>=5.2.0,<7.0.0",
]

setup(
    name="estif-gravity",
    version="6.1.0",
    author="Peter Angelov",
    author_email="tervion@gmail.com",
    description=(
        "ESTIF: Geometric derivation of MOND acceleration, dark matter identity, "
        "and strong-field gravity from 4D hypersurface tilt"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tervion/estif-publication",
    project_urls={
        "Bug Tracker":   "https://github.com/tervion/estif-publication/issues",
        "Documentation": "https://github.com/tervion/estif-publication/tree/main/docs",
        "Source Code":   "https://github.com/tervion/estif-publication",
        "Zenodo":        "https://zenodo.org/records/17418087",
        "Previous Version (ESTIF-FD v1.0)": "https://zenodo.org/records/17261725",
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
            "estif-run=estif_ec_gr_run_simulation:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["data/*.txt", "data/*.tsv"],
    },
    zip_safe=False,
    keywords=[
        "general-relativity",
        "modified-gravity",
        "MOND",
        "dark-matter",
        "dark-energy",
        "SPARC",
        "baryonic-Tully-Fisher",
        "gravitational-waves",
        "LISA",
        "event-horizon-telescope",
        "4D-geometry",
        "hypersurface",
        "emergent-spacetime",
        "cosmological-constant",
        "DESI",
    ],
)

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2
