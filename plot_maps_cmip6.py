import matplotlib.pyplot as plt
# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology maps from CMIP6 models end OBS basedata"

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_mdl(param, exp):
	
	path = '/home/nice/Documentos/data_file/cmip_data/cmip6/historical'
	arq  = '{0}/{1}_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412.nc'.format(path, param)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] * 86400
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	exp  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, exp


def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	map.drawmeridians(np.arange(-85.,-5.,10.), size=4, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.,15.,5.), size=4, labels=[1,0,0,0], linewidth=0.4)
	map.drawcoastlines(linewidth=1, color='k')
	map.drawcountries(linewidth=1, color='k')
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
def plot_maps_clim(mdl1, mdl2, mdl3, mdl4, mdl5, mdl6, mdl7, mdl8, mdl9):
		
	fig = plt.figure(figsize=(8,4))

	levs = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 12]

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(331)
	plt.title('A) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl1[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(332)
	plt.title('B) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl2[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(333)
	plt.title('C) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl3[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(334)
	plt.title('D) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl4[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(335)
	plt.title('E) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl5[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	
	# Plot firt map cmip6 models 
	ax = fig.add_subplot(336)
	plt.title('F) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl6[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(337)
	plt.title('G) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl7[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	
	# Plot firt map cmip6 models 
	ax = fig.add_subplot(338)
	plt.title('H) BCC-CSM2-MR (mm/d)', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl8[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	
	# Plot firt map cmip6 models 
	ax = fig.add_subplot(339)
	plt.title('I) CRU (mm/d)', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	plt.text(-37, -18, u'Reference period: 1985-2004', fontsize=4, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl9[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 
			
	return plt_maps_clim
	
	
# Import regcm exp and cru databases 	   
lat, lon, mdl1 = import_mdl('pr', 'BCC-ESM1')
lat, lon, mdl2 = import_mdl('pr', 'BCC-ESM1')
lat, lon, mdl3 = import_mdl('pr', 'CanESM5')
lat, lon, mdl4 = import_mdl('pr', 'CESM2')
lat, lon, mdl5 = import_mdl('pr', 'MRI-ESM2-0')
lat, lon, mdl6 = import_mdl('pr', 'BCC-CSM2-MR')
lat, lon, mdl7 = import_mdl('pr', 'BCC-CSM2-MR')
lat, lon, mdl8 = import_mdl('pr', 'BCC-CSM2-MR')
lat, lon, mdl9 = import_mdl('pr', 'BCC-CSM2-MR')

# Plot maps with the function
plt_map = plot_maps_clim(mdl1, mdl2, mdl3, mdl4, mdl5, mdl6, mdl7, mdl8, mdl9)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_mean_pr_cmip6_obs_1850-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()

