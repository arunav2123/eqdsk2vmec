# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
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
    description="A Python package to convert EQDSK files to VMEC input format.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)


