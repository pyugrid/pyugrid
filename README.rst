.. image:: https://badge.fury.io/py/pyugrid.png
   :target: http://badge.fury.io/py/pyugrid
.. image:: https://travis-ci.org/pyugrid/pyugrid.svg?branch=master
   :target: https://travis-ci.org/pyugrid/pyugrid


pyugrid
=======

A Python API to utilize data written using the netCDF unstructured grid conventions
[(UGRID)](https://github.com/ugrid-conventions/ugrid-conventions).

Background
----------

For many years, folks in the met-ocean community have been able to exchange data,
model results, etc using the [CF Conventions](http://cfconventions.org/).

However, the Convention does not specify standard ways to work with results from unstructured grid models.
The UGRID effort is an attempt to remedy that.

This package is a small Python implementation of the data model specified by the UGRID project,
as well as code that reads and writes UGRID-compliant netCDF files.

Status
------

The package currently only covers triangular mesh grids,
and has limited functionality for manipulating and visualizing the data.
It does however, provide the ability to read and write netCDF files,
and provides a basic structure an API for adding capability.

Development is managed on GitHub:

https://github.com/pyugrid/pyugrid/
