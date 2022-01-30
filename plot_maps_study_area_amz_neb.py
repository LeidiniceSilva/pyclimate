# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/29/2021"
__description__ = "This script plot study area AMZ_NEB"

import os
import sys
import conda
import numpy as np
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def map_RegCMtopo(ax, lat, lon, topo, latc , lonc ,
        lat_start, lat_end, lon_start, lon_end, fontsize=10):
	"""
	use: map = map_RegCMdomain(ax, latc, lonc, lat_start, lat_end,
							   lon_start, lon_end) # to create a basemap object
	"""
	m = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='i', area_thresh=10000., projection='mill', lon_0=lonc, lat_0=latc, lat_ts=0)	
	llevels = (1, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500)
	x, y = m(lon,lat)

	path = '/home/nice/Documents/github_projects/shp'
	m.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)
	m.readshapefile('{0}/lim_semiarido/LIM_Semiarido_OFICIAL_POLIGONAL'.format(path), 'world', drawbounds=True, color='red', linewidth=1.)
	
	im = m.contourf(x, y, topo, llevels, cmap=plt.cm.RdYlGn_r, extend='max')
	m.drawmapboundary(fill_color='deepskyblue')
	m.drawparallels(np.arange(lat_start, lat_end,  5.), labels=[1,0,0,0], fontsize=fontsize, linewidth=1., color='black')
	m.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[0,0,0,1], fontsize=fontsize, linewidth=1., color='black')                  
	cbar = fig.colorbar(im, drawedges=True, pad=0.05, orientation='horizontal', aspect=40)

	plt.text(7000000, 150000, u'\u25B2 \nN', fontsize=10, fontweight='bold')
	plt.text(2100000, 1500000, u'AMZ', fontsize=10, fontweight='bold')
	plt.text(4600000, 1100000, u'NEB', fontsize=10, fontweight='bold')

	x1,i1 = m(-72,-12)
	x2,i2 = m(-72,0)
	x3,i3 = m(-55,0)
	x4,i4 = m(-55,-12)

	poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.5)
	plt.gca().add_patch(poly1)

	y1,j1 = m(-47,-18)
	y2,j2 = m(-47,-2)
	y3,j3 = m(-35,-2)
	y4,j4 = m(-35,-18)

	poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.5)
	plt.gca().add_patch(poly2)
	
	return m


# Specify directories 
dirnc = '/home/nice/Downloads'
domname = 'reg_amz_neb'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_historical_STS.2005110100.nc'), mode='r')
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
topo = RCMf.variables['topo'][:,:]
mask = RCMf.variables['mask'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()

lat_start  = -20
lat_end    = 10
lon_start  = -85
lon_end    = -15

# Plot study area
fig = plt.figure()
ax = plt.subplot(1,1,1)
m = map_RegCMtopo(ax, lat, lon, topo, latc, lonc, lat_start, lat_end, lon_start, lon_end)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_study_area_amz_neb.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

