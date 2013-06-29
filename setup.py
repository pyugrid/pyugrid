#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup # to support "develop" mode
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy # for the includes for the Cython code


setup(
    name = "pyugrid",
    version = "0.1.1",
    author = "Dharhas Pothina, Alex Crosby, Chris Barker",
    author_email = "Chris.Barker@noaa.gov",
    description = ("A package for working with unstructured grids, and the data on them"),
    license = "BSD",
    keywords = "unstructured numpy models",
    url = "https://github.com/pyugrid/pyugrid",
#    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
        ],
    packages = ["pyugrid", "tests"],
    scripts = [],
    )
