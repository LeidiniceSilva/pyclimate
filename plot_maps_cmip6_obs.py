# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot maps from CMIP6 models end OBS basedata"

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


def import_cmip6_mean(param, model, exp):

	path = '/home/nice/Documentos/data_file/cmip_data/cmip6/historical'
	arq  = '{0}/{1}_Amon_{2}_{3}.nc'.format(path, param, model, exp)	

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

	levs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16]

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(331)
	plt.title(u'A) BCC-CSM2-MR', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl1[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(332)
	plt.title(u'B) CanESM5', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl2[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(333)
	plt.title(u'C) GFDL-ESM4', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl3[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(334)
	plt.title(u'D) GISS-E2-1-H', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl4[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(335)
	plt.title(u'E) HadGEM3', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl5[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(336)
	plt.title(u'F) IPSL-CM6A-LR', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl6[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 

	# Plot firt map cmip6 models 
	ax = fig.add_subplot(337)
	plt.title(u'G) LASG-FGOALS-f3-L', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mdl7[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	ax = fig.add_subplot(338)
	plt.title(u'H) MIROC6', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl8[:,:], levels=levs, latlon=True, cmap=cm.YlGn)

	ax = fig.add_subplot(339)
	plt.title(u'I) MPI-ESM1-2-HR', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mdl9[:,:], levels=levs, latlon=True, cmap=cm.YlGn)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=4) 

	return plt_maps_clim


# Import regcm exp and cru databases 	   
lat, lon, mdl1 = import_cmip6_mean('pr', 'BCC-CSM2-MR', 'historical_r1i1p1f1_gn_185001-201412_remap')
lat, lon, mdl2 = import_cmip6_mean('pr', 'CanESM5', 'historical_r1i1p1f1_gn_185001-201412_remap')
lat, lon, mdl3 = import_cmip6_mean('pr', 'GFDL-ESM4', 'esm-hist_r1i1p1f1_gr1_195001-201412_remap')
lat, lon, mdl4 = import_cmip6_mean('pr', 'GISS-E2-1-H', 'historical_r1i1p1f1_gn_195101-201412_remap')
lat, lon, mdl5 = import_cmip6_mean('pr', 'HadGEM3-GC31-LL', 'historical_r4i1p1f3_gn_195001-201412_remap')
lat, lon, mdl6 = import_cmip6_mean('pr', 'IPSL-CM6A-LR', 'historical_r1i1p1f1_gr_185001-201412_remap')
lat, lon, mdl7 = import_cmip6_mean('pr', 'LASG-FGOALS-f3-L', 'amip_r1i1p1f1_gr_197901-201412_remap')
lat, lon, mdl8 = import_cmip6_mean('pr', 'MIROC6', 'historical_r1i1p1f1_gn_195001-201412_remap')
lat, lon, mdl9 = import_cmip6_mean('pr', 'MPI-ESM1-2-HR', 'highresSST-present_r1i1p1f1_gn_197901-201412_remap')

print(mdl1.shape)
print(mdl2.shape)
print(mdl3.shape)
print(mdl4.shape)
print(mdl5.shape)
print(mdl6.shape)
print(mdl7.shape)
print(mdl8.shape)
print(mdl9.shape)
exit()

# Plot maps with the function
plt_map = plot_maps_clim(mdl1, mdl2, mdl3, mdl4, mdl5, mdl6, mdl7, mdl8, mdl9)

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/papers/cmip6/results'
name_out = 'pyplt_maps_pr_cmip6_obs_1979-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



