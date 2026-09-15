# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="eqdsk2vmec",
    version="0.1.0",
    packages=find_packages(where='source'),
    package_dir={'': 'source'},
    install_requires=[
        "numpy",
        "scipy>=1.8",
        "matplotlib>=3.5",
        "f90nml>=1.4",
    ],
    author="Arunav Kumar",
    license="MIT",
    license_files=["LICENSE"],
    description="A Python package to convert EQDSK files to VMEC input format.",
    long_description=Path(__file__).with_name("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)


