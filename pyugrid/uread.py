import os, sys, netCDF4
import numpy as np
from pyugrid import ugrid.ugrid as ug

def open_cf_todict( filename ):
    nc = netCDF4.Dataset(filename, 'r')
    ncvars = nc.variables
    meshes = {}
    for varname in ncvars.iterkeys():
        meshname = ncvars[varname].getncattr('mesh')
        if (meshname != None) or (meshname not in set(meshes.getkeys())):
            meshatts = ncvars[meshname].getncatts()
            
            ## Make sure that this mesh style is supported in this codebase
            assert meshatts['cf_role'] == 'mesh_topology'
            if meshatts['topology_dimension'] != '2':
                raise error("Unfortuntely, only meshes that are unstructured in 2 dimensions are supported")
     
            ## Grab node coordinates from mesh meta-variable, and pull out the coord values
            node_coordinates = meshatts.get('node_coordinates', None)
            if node_coordinates == None:
                raise error("Unstructured meshes must include node coordinates, specified with the 'node_coordinates' attribute")
            node_coordinates = node_coordinates.split(" ")
            nodes = [None, None]
            for coord in node_coordinates:
                units = ncvars[coord].units
                if 'north' in units:
                    nodes[0] = ncvars[coord][:]
                elif ('east' in units) or ('west' in units):
                    nodes[1] = ncvars[coord][:]
                else:
                    raise error("Node coordinates don't contain 'units' attribute!")  
            
            ## Grab Face and Edge node connectivity arrays
            face_node_conn_name = meshatts.get('face_node_connectivity', None)
            edge_node_conn_name = meshatts.get('edge_node_connectivity', None)
            faces = []
            edges = []
	    if face_node_conn_name != None:
                faces = ncvars[face_node_conn_name][:,:]
            if edge_node_conn_name != None:
                edges = ncvars[edge_node_conn_name][:,:]
                          
            ## Add to dictionary of meshes
            meshes[meshname] = ug(nodes, faces, edges)
    return meshes # Return dictionary of ugrid objects
    