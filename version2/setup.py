"""Setup configuration for BSDWormer."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
requirements = []
if requirements_file.exists():
    requirements = requirements_file.read_text().strip().split('\n')

# Read dev requirements
dev_requirements_file = Path(__file__).parent / "requirements-dev.txt"
dev_requirements = []
if dev_requirements_file.exists():
    dev_requirements = dev_requirements_file.read_text().strip().split('\n')

setup(
    name="bsdwormer",
    version="2.0.0",
    author="Franklin G. Horowitz",
    author_email="frank@horow.net",
    description="Poisson wavelet multiscale edge detection for geophysical potential fields",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RichardScottOZ/BSDWormer",
    packages=find_packages(exclude=["tests", "examples", "docs"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: GIS",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": dev_requirements,
    },
    entry_points={
        "console_scripts": [
            "bsdwormer=bsdwormer.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="geophysics potential-fields edge-detection wavelet fourier",
    project_urls={
        "Bug Reports": "https://github.com/RichardScottOZ/BSDWormer/issues",
        "Source": "https://github.com/RichardScottOZ/BSDWormer",
    },
)
