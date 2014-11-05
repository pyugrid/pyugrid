#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""

import numpy as np
import netCDF4

from pyugrid.ugrid import UGrid, DataSet
from pyugrid.test_examples import two_triangles, twenty_one_triangles

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

    fname = 'temp.nc'
    grid = two_triangles()

    grid.save_as_netcdf(fname)

    ## could be lots of tests here...
    with netCDF4.Dataset(fname) as ds:

        assert nc_has_variable(ds, 'mesh')

        assert nc_var_has_attr_vals(ds, 'mesh', {'cf_role':'mesh_topology',
                                                'topology_dimension' : 2,
                                                'long_name': u'Topology data of 2D unstructured mesh'
                                                })

def test_set_mesh_name():
    fname = 'temp.nc'
    
    grid  =  two_triangles()
    grid.mesh_name = "mesh_2"
    
    grid.save_as_netcdf(fname)

    with netCDF4.Dataset(fname) as ds:
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

        assert nc_has_dimension(ds, "mesh_2_num_node")
        assert nc_has_dimension(ds, "mesh_2_num_edge")
        assert nc_has_dimension(ds, "mesh_2_num_face")
        assert nc_has_dimension(ds, "mesh_2_num_vertices")

        assert not nc_var_has_attr(ds, 'mesh_2', "face_edge_connectivity")


def test_write_with_depths():
    '''
    tests writing a netcdf file with depth data
    '''

    fname = 'temp.nc'

    grid = two_triangles()
    grid.mesh_name='mesh1'

    # create a dataset object for the depths:
    depths = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    depths.attributes['units'] = 'm'
    depths.attributes["standard_name"] = "sea_floor_depth_below_geoid"
    depths.attributes["positive"] = "down"

    grid.add_data(depths)

    grid.save_as_netcdf(fname)

    with netCDF4.Dataset(fname) as ds:

        assert nc_has_variable(ds, 'mesh1')
        assert nc_has_variable(ds, 'depth')

        assert nc_var_has_attr_vals(ds, 'depth', {"coordinates" : "mesh1_node_lon mesh1_node_lat",
                                                  "location" : "node",
                                                  "mesh": "mesh1"})


def test_write_with_velocities():
    '''
    tests writing a netcdf file with velocities on the faces
    '''

    fname = 'temp.nc'

    grid = two_triangles()
    grid.mesh_name = 'mesh2'

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

    # add coordinates for face data
    grid.build_face_coordinates()

    grid.save_as_netcdf(fname)

    with netCDF4.Dataset(fname) as ds:
        assert nc_has_variable(ds, 'mesh2')
        assert nc_has_variable(ds, 'u')
        assert nc_has_variable(ds, 'v')

        assert nc_var_has_attr_vals(ds, 'u', {
                                              "coordinates" : "mesh2_face_lon mesh2_face_lat",
                                              "location" : "face",
                                              "mesh": "mesh2",
                                              })

def test_write_with_edge_data():
    '''
    tests writing a netcdf file with data on the edges (fluxes, maybe?)
    '''
    fname = 'temp.nc'

    grid = two_triangles()
    grid.mesh_name = 'mesh2'

    # create a dataset object for fluxes:
    flux = DataSet('flux', location='edge', data=[0.0, 0.0, 4.1, 0.0, 5.1, ])
    flux.attributes['units'] = 'm^3/s'
    flux.attributes["long_name"] = "volume flux between cells"
    flux.attributes["standard_name"] = "ocean_volume_transport_across_line"

    grid.add_data(flux)
    #add coordinates for edges
    grid.build_edge_coordinates()

    grid.save_as_netcdf(fname)

    with netCDF4.Dataset(fname) as ds:

        assert nc_has_variable(ds, 'mesh2')
        assert nc_has_variable(ds, 'flux')

        assert nc_var_has_attr_vals(ds, 'flux', {
                                              "coordinates" : "mesh2_edge_lon mesh2_edge_lat",
                                              "location" : "edge",
                                              'units' : 'm^3/s',
                                              "mesh": "mesh2",
                                              })
        assert np.array_equal( ds.variables['mesh2_edge_lon'], grid.edge_coordinates[:,0] )
        assert np.array_equal( ds.variables['mesh2_edge_lat'], grid.edge_coordinates[:,1] )

def test_write_with_bound_data():
    '''
    tests writing a netcdf file with data on the boundaries
    suitable for boundary conditions, for example --  (fluxes, maybe?)
    '''
    fname = 'temp.nc'

    grid = two_triangles() # using default mesh name
    # add the boundary definitions:
    grid.boundaries = [(0,1),
                       (0,2),
                       (1,3),
                       (2,3),
                      ]


    # create a dataset object for boundary conditions:
    bnds = DataSet('bnd_cond', location='boundary', data=[0, 1, 0, 0])
    bnds.attributes["long_name"] = "model boundary conditions"
    bnds.attributes["flag_values"] = "0 1"
    bnds.attributes["flag_meanings"] = "no_flow_boundary  open_boundary"

    grid.add_data(bnds)

    grid.save_as_netcdf(fname)

    with netCDF4.Dataset(fname) as ds:

        assert nc_has_variable(ds, 'mesh')
        assert nc_has_variable(ds, 'bnd_cond')

        assert nc_var_has_attr_vals(ds, 'mesh', {
                                              "boundary_node_connectivity" : "mesh_boundary_nodes",
                                              })

        assert nc_var_has_attr_vals(ds, 'bnd_cond', {
                                              "location" : "boundary",
                                              "flag_values" : "0 1",
                                              "flag_meanings" : "no_flow_boundary  open_boundary",
                                              "mesh": "mesh",
                                              })
        ## there should be no coordinates attribute or variable for the boundaries
        ##  as there is no boundaries_coordinates defined
        assert not nc_has_variable(ds, 'mesh_boundary_lon')
        assert not nc_has_variable(ds, 'mesh_boundary_lat')
        assert not nc_var_has_attr(ds, 'bnd_cond', 'coordinates')

def test_write_everything():
    """ An example with all features enabled, and a less trivial grid """

    # use a small, but interesting grid
    fname = 'full_example.nc'

    grid = twenty_one_triangles() # using default mesh name
    grid.build_face_face_connectivity()
    grid.build_edges()

    grid.build_edge_coordinates()
    grid.build_face_coordinates()
    grid.build_boundary_coordinates()

    # depth on the nodes
    depths = DataSet('depth', location='node', data=np.linspace(1,10,20))
    depths.attributes['units'] = 'm'
    depths.attributes["standard_name"] = "sea_floor_depth_below_geoid"
    depths.attributes["positive"] = "down"

    grid.add_data(depths)

    # velocities on the faces:
    u_vel = DataSet('u', location='face', data=np.sin(np.linspace(3,12,21)))
    u_vel.attributes['units'] = 'm/s'
    u_vel.attributes["standard_name"] = "eastward_sea_water_velocity"

    grid.add_data(u_vel)

    # create a dataset object for v velocity:
    v_vel = DataSet('v', location='face', data=np.sin(np.linspace(12,15,21)))
    v_vel.attributes['units'] = 'm/s'
    v_vel.attributes["standard_name"] = "northward_sea_water_velocity"

    grid.add_data(v_vel)

    # fluxes on the edges:
    flux = DataSet('flux', location='edge', data=np.linspace(1000,2000,41))
    flux.attributes['units'] = 'm^3/s'
    flux.attributes["long_name"] = "volume flux between cells"
    flux.attributes["standard_name"] = "ocean_volume_transport_across_line"

    grid.add_data(flux)

    # Some boundary conditions:

    bounds = np.zeros( (19,), dtype=np.uint8 )
    bounds[7] = 1
    bnds = DataSet('bnd_cond', location='boundary', data=bounds)
    bnds.attributes["long_name"] = "model boundary conditions"
    bnds.attributes["flag_values"] = "0 1"
    bnds.attributes["flag_meanings"] = "no_flow_boundary  open_boundary"

    grid.add_data(bnds)

    grid.save_as_netcdf(fname)

    ## now the tests:
    with netCDF4.Dataset(fname) as ds:

        assert nc_has_variable(ds, 'mesh')
        assert nc_has_variable(ds, 'depth')

        assert nc_var_has_attr_vals(ds, 'depth', {"coordinates" : "mesh_node_lon mesh_node_lat",
                                                  "location" : "node"})

        assert nc_has_variable(ds, 'u')
        assert nc_has_variable(ds, 'v')

        assert nc_var_has_attr_vals(ds, 'u', {
                                              "coordinates" : "mesh_face_lon mesh_face_lat",
                                              "location" : "face",
                                              "mesh": "mesh"
                                              })

        assert nc_var_has_attr_vals(ds, 'v', {
                                              "coordinates" : "mesh_face_lon mesh_face_lat",
                                              "location" : "face",
                                              "mesh": "mesh",
                                              })

        assert nc_has_variable(ds, 'flux')

        assert nc_var_has_attr_vals(ds, 'flux', {
                                              "coordinates" : "mesh_edge_lon mesh_edge_lat",
                                              "location" : "edge",
                                              'units' : 'm^3/s',
                                              "mesh": "mesh",
                                              })
        assert nc_has_variable(ds, 'mesh')
        assert nc_has_variable(ds, 'bnd_cond')

        assert nc_var_has_attr_vals(ds, 'mesh', {
                                              "boundary_node_connectivity" : "mesh_boundary_nodes",
                                              })

        assert nc_var_has_attr_vals(ds, 'bnd_cond', {
                                              "location" : "boundary",
                                              "flag_values" : "0 1",
                                              "flag_meanings" : "no_flow_boundary  open_boundary",
                                              "mesh": "mesh",
                                              })
    # and make sure pyugrid can reload it!
    grid = UGrid.from_ncfile(fname,load_data=True)
    # and that some things are the same:
    # note:  more testing might be good here...
    #        maybe some grid comparison functions? 

    assert grid.mesh_name == 'mesh'

    print "grid data:", grid.data
    assert len(grid.nodes) == 20


    depth = grid.data['depth']
    assert depth.attributes['units'] == 'm'

    u = grid.data['u']
    assert u.attributes['units'] == 'm/s'


if __name__ == "__main__":
    # run the tests:

    test_simple_write()
    test_set_mesh_name()
    test_write_with_depths()
    test_write_with_velocities()
    test_write_with_edge_data()
   

 