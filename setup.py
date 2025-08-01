from setuptools import setup, find_packages

setup(
    name="eqdsk2vmec",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "freeqdsk",
    ],
    author="Arunav Kumar",
    author_email="arunavk@mit.edu",
    description="A Python package to convert EQDSK files to VMEC input format.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)


