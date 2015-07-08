# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ##Test out UGRID-0.9 compliant unstructured grid model datasets with PYUGRID

# <codecell>

from __future__ import (absolute_import, division, print_function)

import matplotlib.tri as tri
import datetime as dt

# <codecell>

import cartopy.crs as ccrs
import iris
iris.FUTURE.netcdf_promote = True
import pyugrid

# <codecell>

#ADCIRC
#url =  'http://comt.sura.org/thredds/dodsC/data/comt_1_archive/inundation_tropical/UND_ADCIRC/Hurricane_Ike_3D_final_run_with_waves'

#FVCOM
#url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'

#SELFE
url = 'http://comt.sura.org/thredds/dodsC/data/comt_1_archive/inundation_tropical/VIMS_SELFE/Hurricane_Ike_2D_final_run_with_waves'

# <codecell>

cube = iris.load_cube(url,'sea_surface_height_above_geoid')

# <codecell>

print(cube)

# <codecell>

# Desired time for snapshot
# ....right now (or some number of hours from now) ...
start = dt.datetime.utcnow() + dt.timedelta(hours=6)
# ... or specific time (UTC)
#start = dt.datetime(2013,3,2,15,0,0)

# <codecell>

ug = pyugrid.UGrid.from_ncfile(url)

# What's in there?
#print "There are %i nodes"%ug.nodes.shape[0]
#print "There are %i edges"%ug.edges.shape[0]
#print "There are %i faces"%ug.faces.shape[0]

# <codecell>

cube.mesh = ug
cube.mesh_dimension = 1  # (0:time,1:node)

# <codecell>

lon = cube.mesh.nodes[:,0]
lat = cube.mesh.nodes[:,1]
nv = cube.mesh.faces

# <codecell>

triang = tri.Triangulation(lon,lat,triangles=nv)

# <codecell>

# skip trying to find the closest time index to requested time, because it's messy
ind = -1 # just take the last time index for now
zcube = cube[ind]

# <codecell>

figure(figsize=(12,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-90, -60, 5, 50])
ax.coastlines()
levs=arange(-1,5,.2)
tricontourf(triang, zcube.data, levels=levs)
colorbar()
tricontour(triang, zcube.data, colors='k',levels=levs)
tvar = cube.coord('time')
tstr = tvar.units.num2date(tvar.points[ind])
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = False
gl.ylabels_right = False
title('%s: Elevation (m): %s' % (zcube.attributes['title'],tstr));

# <codecell>


