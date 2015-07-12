from __future__ import (absolute_import, division, print_function)

import os
from setuptools import setup
from setuptools.command.test import test as TestCommand


rootpath = os.path.abspath(os.path.dirname(__file__))
long_description = open(os.path.join(rootpath, 'README.rst')).read()


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        import sys
        sys.exit(errno)


def version():
    """Get the version number."""
    with open(os.path.join(rootpath, "VERSION.txt")) as v:
        _version = v.read()
    return _version.strip()


__version__ = version()


def set_version(filename, new_string=__version__):
    with open(filename) as f:
        lines = f.readlines()
    for k, line in enumerate(lines):
        if (line.startswith('__version__')):
            lines[k] = "__version__ = '{}'\n".format(__version__)
            break
    with open(filename, 'w') as f:
        f.writelines(lines)

set_version(os.path.join(rootpath, 'pyugrid', '__init__.py'))

setup(
    name="pyugrid",
    version=__version__,
    author="Chris Barker, Chris Calloway, Rich Signell",
    author_email="Chris.Barker@noaa.gov",
    description=("A package for working with triangular unstructured grids, "
                 "and the data on them"),
    license="BSD",
    keywords="unstructured numpy models",
    url="https://github.com/pyugrid/pyugrid",
    long_description=long_description,
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
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        ],
    packages=["pyugrid", "test"],
    scripts=[],
    )
