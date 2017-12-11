.. image:: https://badge.fury.io/py/pyugrid.png
   :target: http://badge.fury.io/py/pyugrid
.. image:: https://travis-ci.org/pyugrid/pyugrid.svg?branch=master
   :target: https://travis-ci.org/pyugrid/pyugrid


pyugrid
=======

A Python API to utilize data written using the netCDF unstructured grid conventions:
[(UGRID)](https://github.com/ugrid-conventions/ugrid-conventions).


Background
----------

For many years, folks in the met-ocean community have been able to exchange data,
model results, etc using the [CF Conventions](http://cfconventions.org/).

However, the Convention does not specify standard ways to work with results
from unstructured grid models.  The UGRID effort is an attempt to remedy that.

This package is a Python implementation of the data model specified by the UGRID project.

It provides code that reads and writes UGRID-compliant netCDF files, a start on
code to read/write other formats, and code to work with unstructured grids.


Status
------

**NOTE:** This may be the last release of pyugrid -- continued development will be taking place as part of the "gridded" project:  https://github.com/NOAA-ORR-ERD/gridded

The package currently covers triangular mesh grids, quad grids, and mixed traingle/quad grids.

It provides some limited functionality for manipulating and visualizing the data.

It provides functionality for accessing and interpolating fields on the Grid.

It also provides the ability to read and write netCDF files, and provides a basic
structure an API for adding capability.

Development is managed on GitHub:

https://github.com/pyugrid/pyugrid/


