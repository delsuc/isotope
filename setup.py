#!/usr/bin/env python 
# encoding: utf-8

"""
setup.py

"""
import setuptools

version = '1.1'        
# get description
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='isotope',
    version=version,
    author="M-A. Delsuc",
    author_email="madelsuc@unistra.fr",
    description="This module allows the precise computation of isotopic distributions, average mass, and monoisotopic masses of arbitrary chemical formula",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/delsuc/isotope",
    packages=setuptools.find_packages(),
    include_package_data = True,
    package_data={ "elements": ["elements.asc"]},
    license="CeCILL-2.1",
    provides=["isotope"],
    requires=["matplotlib", "numpy", 'scipy'],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
)
