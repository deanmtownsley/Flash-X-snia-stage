"""Build and installation script for myproject."""

# standard libraries
import re
from setuptools import setup, find_packages

# get long description from README.rst
with open("README.md", mode="r") as readme:
    long_description = readme.read()

# get package metadata by parsing __meta__ module
with open("flashparser/__meta__.py", mode="r") as source:
    content = source.read().strip()
    metadata = {
        key: re.search(key + r'\s*=\s*[\'"]([^\'"]*)[\'"]', content).group(1)
        for key in [
            "__pkgname__",
            "__version__",
            "__description__",
        ]
    }

# core dependancies
DEPENDENCIES = ["numpy"]

setup(
    name=metadata["__pkgname__"],
    version=metadata["__version__"],
    description=metadata["__description__"],
    packages=find_packages(where="./"),
    package_dir={"": "./"},
    include_package_data=True,
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python :: 3.8",
    ],
    install_requires=DEPENDENCIES,
)
