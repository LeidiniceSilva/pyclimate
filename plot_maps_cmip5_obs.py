# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot maps from CMIP5 models end OBS basedata"

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap


def import_cmip5_mean(model):
	
	param = 'tas' # pr or tas
	area  = 'amz_neb' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:] 
	sim   = np.nanmean(var[:][:,:,:], axis=0)
			
	return lat, lon, sim


def import_obs_mean(database):
	
	param = 'tmp' # pre or tmp
	area  = 'amz_neb' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, 
	database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]	
	obs   = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, obs


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


def plot_maps_clim(mdl1, mdl2, mdl3, mdl4):

	fig = plt.figure()

	levs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16]
	levs = [12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(141)
	plt.title(u'A) CanESM2', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl1[:,:], levels=levs, latlon=True, cmap=cm.Reds)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(142)
	plt.title(u'B) HadGEM2-ES', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl2[:,:], levels=levs, latlon=True, cmap=cm.Reds)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(143)
	plt.title(u'C) MPI-ESM-MR', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl3[:,:], levels=levs, latlon=True, cmap=cm.Reds)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(144)
	plt.title(u'D) NorESM1-M', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl4[:,:], levels=levs, latlon=True, cmap=cm.Reds)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 
	
	return plt_maps_clim


# Import regcm exp and cru databases 	   
lat, lon, mdl1 = import_cmip5_mean('CanESM2')
lat, lon, mdl2 = import_cmip5_mean('HadGEM2-ES')
lat, lon, mdl3 = import_cmip5_mean('MPI-ESM-MR')
lat, lon, mdl4 = import_cmip5_mean('NorESM1-M')

# Plot maps with the function
plt_map = plot_maps_clim(mdl1, mdl2, mdl3, mdl4)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_tmp_cmip5_obs_1975-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
