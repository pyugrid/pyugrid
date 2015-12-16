#!/usr/bin/env python

"""
UVar object, used to hold variables that are associated with a ugrid
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from util import asarraylike

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
            raise ValueError("location must be one of: 'node', 'edge', 'face', 'boundary', or 'none'")
 
        self.location = location 

        if data is None:
            self._data = np.zeros((0,), dtype=np.float64) # could be any data type
        else:
            self._data = asarraylike(data)

        if attributes is None and data is not None and hasattr(data,'__dict__'):
            self.update(data.__dict__)
        else:
            self.update(attributes)

    def update(self, attr):
        """

        :param attr: Dict containing attributes to be added to the object
        """
        for key,val in attr.items():
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
    def dimensions(self):
        return zip(self.data.dimensions, self.data.shape)

    @property
    def shape(self):
        return self.data.shape

    def __getitem__(self, item):
        """
        Transfers responsibility to the data's __getitem__
        """
        return self.data.__getitem__(item)

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
        self.name = name

        if location not in ['node', 'edge', 'face', 'boundary', 'none']:
            raise ValueError("location must be one of: 'node', 'edge', 'face', 'boundary', or 'none'")

        self.location = location

        if len(data) == 1:
            raise ValueError("UMVar need at least 2 data sources of the same size and shape")

        shape = data[0].shape
        if not all([d.shape == shape for d in data]):
            raise ValueError("All data sources must be the same size and shape")

        for d in data:
            setattr(self, d.name, d)

        self.variables = [d.name for d in data]

    def dimensions(self):
        self.__getattribute__(self.variables[0]).shape

    def __getitem__(self, item):
        return np.column_stack((self.__getattribute__(var).__getitem__(item) for var in self.variables))

if __name__ == "__main__":
    import netCDF4 as ncdf
    df = ncdf.Dataset('../test/data/21_tri_mesh.nc')
    u = UVar('EW_water_velocity', 'node', df['u'])
    v = UVar('NS_water_velocity', 'node', df['v'])
    vels = UMVar('velocity', 'node', [u,v])
    pass