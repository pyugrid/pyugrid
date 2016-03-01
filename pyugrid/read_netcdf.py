#!/usr/bin/env python

"""
code to read the netcdf unstructured grid standard:

https://github.com/ugrid-conventions/ugrid-conventions/

This code is called by the UGrid class to load into a UGRID object.

    NOTE: passing the UGrid object in to avoid circular references,
    while keeping the netcdf reading code in its own file.

"""

from __future__ import (absolute_import, division, print_function)


__all__ = ['find_mesh_names',
           'is_valid_mesh']


def find_mesh_names(nc):
    """
    Find all the meshes in a `netcCDF4.Dataset`.

    :param nc: the `netCDF4.Dataset` object to look for mesh names in

    NOTE: checks for 2D topology_dimension.
    """
    mesh_names = []
    for varname in nc.variables.keys():
        if is_valid_mesh(nc, varname):
                    mesh_names.append(varname)
    return mesh_names


def is_valid_mesh(nc, varname):
    """
    Determine if the given variable name is a valid mesh definition.

    :param nc: a `netCDF4.Dataset` to check

    :param varname: name of the candidate mesh variable

    """
    try:
        mesh_var = nc.variables[varname]
    except KeyError:
        return False
    try:
        if (mesh_var.cf_role.strip() == 'mesh_topology' and
           int(mesh_var.topology_dimension) == 2):
            return True
    except AttributeError:
        return False

# Defining properties of various connectivity arrays
# so that the same code can load all of them.
grid_defs = [{'grid_attr': 'faces',  # Name in UGrid object.
              'role': 'face_node_connectivity',  # Name in mesh variable.
              'num_ind': 3,  # Number of idx (3 for faces, 2 for segments).
              },
             {'grid_attr': 'face_face_connectivity',  # Name in UGrid object.
              'role': 'face_face_connectivity',  # Name in mesh variable.
              'num_ind': 3,  # Number of idx (3 for faces, 2 for segments).
              },
             {'grid_attr': 'boundaries',  # Name in UGrid object.
              'role': 'boundary_node_connectivity',  # Name in mesh var.
              'num_ind': 2,  # Number of idx (3 for faces, 2 for segments).
              },
             {'grid_attr': 'edges',  # Name in UGrid object.
              'role': 'edge_node_connectivity',  # Name in mesh variable.
              'num_ind': 2,  # Number of idx (3 for faces, 2 for segments).
              },
             ]

# Definitions for various coordinate arrays.
coord_defs = [{'grid_attr': 'nodes',  # Attribute name in UGrid object.
               'role': 'node_coordinates',  # Attribute name in mesh variable.
               'required': True,  # Is this required?
               },
              {'grid_attr': 'face_coordinates',  # Name in UGrid object.
               'role': 'face_coordinates',  # Name in mesh variable.
               'required': False,  # Is this required?
               },
              {'grid_attr': 'edge_coordinates',  # Name in UGrid object.
               'role': 'edge_coordinates',  # Name in mesh variable.
               'required': False,  # is this required?
               },
              {'grid_attr': 'boundary_coordinates',  # Name in UGrid object.
               'role': 'boundary_coordinates',  # Name in mesh variable.
               'required': False,  # Is this required?
               },
              ]
