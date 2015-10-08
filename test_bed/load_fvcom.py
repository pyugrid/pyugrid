#!usr/bin/env python

"""
test to see how to write a custom ugrid file reader
"""

import numpy as np
import netCDF4
import pyugrid

def load_fvcom_gnome(filename):
    """
    load an fvcom/gnome netcdf file

    :param filename: the path to the file to load

    :returns: an UGRID object with the data from the file
    """
    nc = netCDF4.Dataset("small_trigrid_example.nc")

    # check if it's the file type we are looking for
    if nc.getncattr('grid_type').lower() != 'triangular':
        raise ValueError('This does not appear to be a valid triangular grid file:\nIt does not have the "grid_type"="Triangular" global attribute')

    # create an empty UGrid:
    ug = pyugrid.UGrid()
    # load the nodes
    lon = nc.variables['lon']
    lat = nc.variables['lat']

    print lon.shape
    num_nodes = lon.shape[0]
    ug.nodes = np.zeros((num_nodes, 2), dtype=lon.dtype)
    ug.nodes[:,0] = lon[:]
    ug.nodes[:,1] = lat[:]

    # load the faces
    ug.faces = nc.variables['nv'][:].T

    # load the connectivity array
    ug.face_face_connectivity = nc.variables['nbe'][:].T

    # load the center points of the faces
    ug.face_coordinates = np.zeros((len(ug.faces), 2), dtype=lon.dtype)
    ug.face_coordinates[:,0] = nc.variables['lonc'][:]
    ug.face_coordinates[:,1] = nc.variables['latc'][:]


    bounds = nc.variables['bnd']
    ug.boundaries = bounds[:,:2] # ignoring the second two fields -- what are they???

    return ug

def load_from_varnames(filename, names_mapping):
    """
    load a UGrid from a netcdf file where the roles are defined by the names of the variables

    :param filename: the names of the file to load (or opendap url)

    :param names_mapping: dict that maps the variable names to the UGRid components

    """

    nc = netCDF4.Dataset("small_trigrid_example.nc")

    # check if it's the file type we are looking for
    try:
        name, value = names_mapping['attribute_check']
        if nc.getncattr(name).lower() != value:
            raise ValueError('This does not appear to be a valid triangular grid file:\nIt does not have the "grid_type"="Triangular" global attribute')
    except KeyError:
        pass

    # create an empty UGrid:
    ug = pyugrid.UGrid()
    # load the nodes
    lon = nc.variables['lon']
    lat = nc.variables['lat']

    print lon.shape
    num_nodes = lon.shape[0]
    ug.nodes = np.zeros((num_nodes, 2), dtype=lon.dtype)
    ug.nodes[:,0] = lon[:]
    ug.nodes[:,1] = lat[:]

    # load the faces
    ug.faces = nc.variables['nv'][:].T

    # load the connectivity array
    ug.face_face_connectivity = nc.variables['nbe'][:].T

    # load the center points of the faces
    ug.face_coordinates = np.zeros((len(ug.faces), 2), dtype=lon.dtype)
    ug.face_coordinates[:,0] = nc.variables['lonc'][:]
    ug.face_coordinates[:,1] = nc.variables['latc'][:]


    bounds = nc.variables['bnd']
    ug.boundaries = bounds[:,:2] # ignoring the second two fields -- what are they???

    return ug

if __name__ == "__main__":

    names_mapping = {
                     'attribute_check': ('grid_type','triangular'),
                     'nodes': ('lon', 'lat'),

                    }


    ug = load_from_varnames("small_trigrid_example.nc", names_mapping)

    print ug




