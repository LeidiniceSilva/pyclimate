# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2019"
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


def import_mdl(param, model):
	
	path = '/home/nice/Documentos/data_file/cmip_data/cmip6/historical'
	arq  = '{0}/{1}_Amon_{2}_historical_r1i1p1f1_gn_185001-201412.nc'.format(path, param, model)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] * 86400
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	mdl  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mdl


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
	
	map = Basemap(projection='cyl', llcrnrlon=np.nanmin(lon),urcrnrlon=np.nanmax(lon), llcrnrlat=lat.min(),urcrnrlat=lat.max(), resolution='c')
	map.drawmeridians(np.arange(-180.,180.,60.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-90.,90.,30.), size=6, labels=[1,0,0,0], linewidth=0.4)
	map.drawcoastlines(linewidth=1, color='k')
	map.drawcountries(linewidth=1, color='k')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	

def plot_maps_clim(mdl_mean, mdl_bias):
		
	fig = plt.figure(figsize=(8,4))

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(221)
	plt.title('A) Ensemble Mean Rainfall (mm/d)', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	levs = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mdl_mean[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(222)
	plt.title('B) CRU Mean Rainfall (mm/d)', fontsize=8, fontweight='bold')
	levs = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mdl_mean[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(223)
	plt.title('C) Mean Bias Error', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	levs = [-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6]
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mdl_mean[:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(224)
	plt.title('D) Root Mean Square Error', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	levs = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20]
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mdl_bias[:,:], levels=levs, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	plt.tight_layout()
		
	return plot_maps_clim
		

# Import regcm exp and cru databases 	   
lat, lon, mdl1 = import_mdl('pr', 'BCC-ESM1')

# Compute and plot bias from regcm exp and cru database
mdl_mean = mdl1
mdl_bias = mdl1

plt_map = plot_maps_clim(mdl_mean, mdl_bias)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_bias_pr_cmip6_obs_1850-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
