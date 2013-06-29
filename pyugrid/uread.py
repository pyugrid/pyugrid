import netCDF4
import numpy as np
from ugrid import ugrid as ug

def open_cf_todict( filename ):
    nc = netCDF4.Dataset(filename, 'r')
    ncvars = nc.variables
    meshes = {}
    for varname in ncvars.iterkeys():
        try:
            meshname = ncvars[varname].getncattr('mesh')
        except AttributeError:
            meshname = None
        if (meshname != None) and (meshname not in set(meshes.viewkeys())):
            meshatt_names = ncvars[meshname].ncattrs()
            
            ## Make sure that this mesh style is supported in this codebase
            meshatts = {}
            for attname in meshatt_names:
                meshatts[attname] = ncvars[meshname].getncattr(attname)
            assert meshatts['cf_role'] == 'mesh_topology'
            #if meshatts['topology_dimension'] != '2':
            #    raise ValueError("Unfortuntely, only meshes that are unstructured in 2 dimensions are supported")
            
            ## Grab node coordinates from mesh meta-variable, and pull out the coord values
            node_coordinates = meshatts.get('node_coordinates', None)
            if node_coordinates == None:
                raise AttributeError("Unstructured meshes must include node coordinates, specified with the 'node_coordinates' attribute")
            node_coordinates = node_coordinates.split(" ")
            nodes = [None, None]
            for coord in node_coordinates:
                units = ncvars[coord].units
                if 'north' in units:
                    nodes[0] = ncvars[coord][:]
                elif ('east' in units) or ('west' in units):
                    nodes[1] = ncvars[coord][:]
                else:
                    raise AttributeError("Node coordinates don't contain 'units' attribute!")  
            
            ## Grab Face and Edge node connectivity arrays
            face_node_conn_name = meshatts.get('face_node_connectivity', None)
            edge_node_conn_name = meshatts.get('edge_node_connectivity', None)
            faces = []
            edges = []
	    if face_node_conn_name != None:
                faces = ncvars[face_node_conn_name][:,:]
                index_base = np.min(np.min(faces))
                if index_base  >= 1:
                    faces = faces - index_base
            try:
                if edge_node_conn_name != None:
                    edges = ncvars[edge_node_conn_name][:,:]
                    index_base = np.min(np.min(edges))
                    if index_base >= 1:
                        edges = edges - index_base
            except:
                pass #TODO: Generate edge node topology if none exists, perhaps optional
  
            ## Add to dictionary of meshes
            meshes[meshname] = ug(nodes, faces, edges)
    return meshes # Return dictionary of ugrid objects
    
