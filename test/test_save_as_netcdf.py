#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""

import netCDF4

from pyugrid.ugrid import UGrid, DataSet
from pyugrid.test_examples import two_triangles

# code to check netcdf files for stuff:
def nc_has_variable(ds, var_name):
    """
    checks that a netcdf file has the given variable defined
    
    :param ds: a netCDF4 Dataset object, or a netcdf file name

    """
    if not isinstance(ds, netCDF4.Dataset):
        ds = netCDF4.Dataset(ds)

    if ds.variables.has_key(var_name):
        return True
    else:
        print var_name, " is not a variable in the dataset"
        return False

def nc_has_dimension(ds, dim_name):
    """
    checks that a netcdf file has the given dimension defined
    
    :param ds: a netCDF4 Dataset object, or a netcdf file name

    """
    if not isinstance(ds, netCDF4.Dataset):
        ds = netCDF4.Dataset(ds)

    if ds.dimensions.has_key(dim_name):
        return True
    else:
        print dim_name, " is not a dimension in the dataset"
        return False


def nc_var_has_attr(ds, var_name, att_name):
    """
    checks that the variable, var_name, has the attribute, att_name
    """
    if not isinstance(ds, netCDF4.Dataset):
        ds = netCDF4.Dataset(ds)

    try:
        getattr(ds.variables[var_name], att_name)
        return True
    except AttributeError:
        print att_name, "is not in the var:", var_name
        return False    

def nc_var_has_attr_vals(ds, var_name, att_dict):
    """
    checks that the variable, var_name, as teh attribtes (and values) in the att_dict
    """
    if not isinstance(ds, netCDF4.Dataset):
        ds = netCDF4.Dataset(ds)

    for key, val in att_dict.items():
        try:
            if val != getattr(ds.variables[var_name], key):
                print "attribute:", key
                print "expected val:", val
                print "val in file:", repr( getattr(ds.variables[var_name], key) )
                return False
        except AttributeError:
            print key, "is not an attribute of var:", var_name
            return False
    return True


def test_simple_write():

    fname = 'two_triangles.nc'
    grid = two_triangles()

    grid.save_as_netcdf(fname)

    ## could be lots of tests here...
    ds = netCDF4.Dataset(fname)
    assert nc_has_variable(ds, 'mesh')

    assert nc_var_has_attr_vals(ds, 'mesh', {'cf_role':'mesh_topology',
                                            'topology_dimension' : 2,
                                            'long_name': u'Topology data of 2D unstructured mesh'
                                            })
    ds.close()

def test_set_mesh_name():
    fname = 'two_triangles2.nc'
    grid  =  two_triangles()

    grid.save_as_netcdf(fname, mesh_name = "mesh_2")

    ds = netCDF4.Dataset(fname)
    assert nc_has_variable(ds, 'mesh_2')

    assert nc_var_has_attr_vals(ds, 'mesh_2', {'cf_role':'mesh_topology',
                                               'topology_dimension' : 2,
                                               'long_name': u'Topology data of 2D unstructured mesh'
                                               })

    assert nc_var_has_attr_vals(ds, 'mesh_2', {'cf_role':'mesh_topology',
                                               'topology_dimension' : 2,
                                               'long_name': u'Topology data of 2D unstructured mesh',
                                               'node_coordinates': 'mesh_2_node_lon mesh_2_node_lat',
                                               })

    assert nc_has_variable(ds, 'mesh_2_node_lon')
    assert nc_has_variable(ds, 'mesh_2_node_lat')
    assert nc_has_variable(ds, 'mesh_2_face_nodes')
    assert nc_has_variable(ds, 'mesh_2_edge_nodes')

    assert nc_has_dimension(ds, "mesh_2_num_nodes")
    assert nc_has_dimension(ds, "mesh_2_num_edges")
    assert nc_has_dimension(ds, "mesh_2_num_faces")
    assert nc_has_dimension(ds, "mesh_2_num_vertices")

    assert not nc_var_has_attr(ds, 'mesh_2', "face_edge_connectivity")

    ds.close()

def test_write_with_depths():
    '''
    tests writting a netcdf file with depth data
    '''

    fname = 'two_triangles_depth.nc'

    grid = two_triangles()

    # create a dataset object for the depths:
    depths = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    depths.attributes['units'] = 'm'
    depths.attributes["standard_name"] = "sea_floor_depth_below_geoid"
    depths.attributes["positive"] = "down"

    grid.add_data(depths)

    grid.save_as_netcdf(fname)

    ds = netCDF4.Dataset(fname)

    assert nc_has_variable(ds, 'mesh')
    assert nc_has_variable(ds, 'depth')

    assert nc_var_has_attr_vals(ds, 'depth', {"coordinates" : "mesh_node_lon mesh_node_lat",
                                              "location" : "node"})


def test_write_with_velocities():
    '''
    tests writting a netcdf file with velocities on the faces
    '''

    fname = 'two_triangles_vel.nc'

    grid = two_triangles()

    # create a dataset object for u velocity:
    u_vel = DataSet('u', location='face', data=[1.0, 2.0])
    u_vel.attributes['units'] = 'm/s'
    u_vel.attributes["standard_name"] = "eastward_sea_water_velocity"

    grid.add_data(u_vel)

    # create a dataset object for v velocity:
    v_vel = DataSet('v', location='face', data=[3.2, 4.3])
    v_vel.attributes['units'] = 'm/s'
    v_vel.attributes["standard_name"] = "northward_sea_water_velocity"

    grid.add_data(v_vel)

    grid.save_as_netcdf(fname)

    ds = netCDF4.Dataset(fname)

    assert nc_has_variable(ds, 'mesh')
    assert nc_has_variable(ds, 'u')
    assert nc_has_variable(ds, 'v')

    assert nc_var_has_attr_vals(ds, 'u', {
                                          #"coordinates" : "mesh_node_lon mesh_node_lat",
                                          "location" : "face",
                                          })


    assert False



   

 