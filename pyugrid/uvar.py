#!/usr/bin/env python

"""
UVar object, used to hold variables that are associated with a ugrid
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from collections import OrderedDict

try:
    from .util import asarraylike, isarraylike
except ValueError:
    from util import asarraylike, isarraylike


class UVar(object):
    """
    A class to hold a variable associated with the UGrid. Data can be on the 
    nodes, edges, etc. -- "UGrid Variable"

    It holds an array of the data, as well as the attributes associated
    with that data  -- this is mapped to a netcdf variable with 
    attributes(attributes get stored in the netcdf file)
    """

    def __init__(self, name, location='none', data=None, attributes=None):
        """
        create a UVar object
        :param name: the name of the data (depth, u_velocity, etc.)
        :type name: string

        :param location: the type of grid element the data is associated with:
                         'node', 'edge', or 'face' the data is assigned to

        :param data: the data
        :type data: 1-d numpy array or array-like object ().
                    IF you have a list or tuple, it should be converted or something compatible (list, etc.)        
        """
        self.name = name

        if location not in ['node', 'edge', 'face', 'boundary', 'none']:
            raise ValueError(
                "location must be one of: 'node', 'edge', 'face', 'boundary', or 'none'")

        self.location = location

        if data is None:
            # could be any data type
            self._data = np.zeros((0,), dtype=np.float64)
        else:
            self._data = asarraylike(data)

        if attributes is None and data is not None and hasattr(data, '__dict__'):
            self.update(data.__dict__)
        else:
            self.update(attributes)
        self._cache = OrderedDict()

    def update(self, attr):
        """

        :param attr: Dict containing attributes to be added to the object
        """
        for key, val in attr.items():
            setattr(self, key, val)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = asarraylike(data)

    @data.deleter
    def data(self):
        self._data = self._data = np.zeros((0,), dtype=np.float64)

    @property
    def shape(self):
        return self.data.shape

    @property
    def max(self):
        return np.max(self._data)

    @property
    def min(self):
        return np.min(self._data)

    @property
    def dtype(self):
        return self.data.dtype

    @property
    def ndim(self):
        return self.data.ndim

    def __getitem__(self, item):
        """
        Transfers responsibility to the data's __getitem__ if not cached
        """
        rv = None
        if str(item) in self._cache:
            rv = self._cache[str(item)]
        else:
            rv = self._data.__getitem__(item)
            self._cache[str(item)] = rv
            if len(self._cache) > 3:
                self._cache.popitem(last=False)
        return rv

    def __str__(self):
        return "UVar object: {0:s}, on the {1:s}s, and {2:d} data points\nAttributes: {3}".format(self.name, self.location, len(self.data), self.attributes)

    def __len__(self):
        return len(self.data)


class UMVar(object):
    """
    A class to group multiple UVars (or other data sources) and retrieve common information. All the variables
    grouped in this class must have the same shape, location, and unique names.
    TODO: Add attribues that all grouped variables have in common to the UMVar?
    """

    def __init__(self, name, location='none', data=None, attributes=None):
        """
        :param name: the name of the data (depth, u_velocity, etc.)
        :type name: string

        :param location: the type of grid element the data is associated with:
                         'node', 'edge', or 'face' the data is assigned to

        :param data: the data
        :type data: list-like of data sources that satisfy the conditions of util.asarraylike. All data sources
        must have the same shape.
        Examples: netCDF Dataset, numpy array
        """
        self.name = name

        if location not in ['node', 'edge', 'face', 'boundary', 'none']:
            raise ValueError(
                "location must be one of: 'node', 'edge', 'face', 'boundary', or 'none'")

        self.location = location

        if len(data) == 1:
            raise ValueError(
                "UMVar need at least 2 data sources of the same size and shape")

        if not all([isarraylike(d) for d in data]):
            raise ValueError("Data must satisfy isarraylike or be a UVar")

        self.shape = data[0].shape
        if not all([d.shape == self.shape for d in data]):
            raise ValueError(
                "All data sources must be the same size and shape")

        for d in data:
            setattr(self, d.name, d)

        self.variables = [d.name for d in data]
        self._cache = OrderedDict()

    def add_var(self, var):
        if var.shape != self.shape:
            raise ValueError(
                'Variable {0} has incorrect shape {1}'.format(var.name, var.shape))
        if var.name in self.variables:
            raise ValueError(
                'Variable {0} already exists in UMVar'.format(var.name))
        self.variables.append(var.name)
        setattr(self, var.name, var)

    def __getitem__(self, item):
        if str(item) in self._cache:
            return self._cache[str(item)]
        else:
            rv = np.ma.column_stack(
                [self.__getattribute__(var).__getitem__(item) for var in self.variables])
            self._cache[str(item)] = rv
            if len(self._cache) > 3:
                self._cache.popitem(last=False)
            return rv


if __name__ == "__main__":
    import netCDF4 as ncdf
    df = ncdf.Dataset('../test/files/21_tri_mesh.nc')
    u = UVar('EW_water_velocity', 'node', df['u'])
    v = UVar('NS_water_velocity', 'node', df['v'])
    vels = UMVar('velocity', 'node', [u, v])
    vels.add_var(u)
    pass
