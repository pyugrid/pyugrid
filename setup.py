#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup # to support "develop" mode
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
from setuptools.command.test import test as TestCommand

#import numpy # for the includes for the Cython code

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name="pyugrid",
    version="0.1.2",
    author="Chris Barker, Chris Calloway, Rich Signell",
    author_email="Chris.Barker@noaa.gov",
    description=("A package for working with triangular unstructured grids, and the data on them"),
    license="BSD",
    keywords="unstructured numpy models",
    url="https://github.com/pyugrid/pyugrid",
    long_description=open('README.rst').read(),
    install_requires=[
        'numpy',
        'netCDF4',
        ],
    tests_require=[
        'pytest>=2.3.2',
        ],
    cmdclass={'test': PyTest},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
        ],
    packages=["pyugrid", "test"],
    scripts=[],
    )
