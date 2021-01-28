# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot climatology maps from Reg and Had models end obs database"

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

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from comp_statist_indices import compute_added_value


def import_obs(var):
	
	path = '/home/nice/Downloads'
	arq  = '{0}/cru_ts4.04.1901.2019.{1}.dat_lonlat.nc'.format(path, var)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	clim_8110  = np.nanmean(var[:][960:1319,:,:], axis=0)
	clim_0119  = np.nanmean(var[:][:,:,:], axis=0)
	clim_anomaly = clim_8110 - clim_0119
	
	value_8110 = np.nanmean(clim_8110, axis=1)
	value_0119 = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

	#~ year_anomaly = []
	#~ for i in range(0, 1428):
		#~ anomaly = value_0119[i] - value_8110
		#~ print(anomaly)
		#~ exit()
		#~ year_anomaly.append(np.squeeze(anomaly))

	return lat, lon, clim_anomaly
	

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
	map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4)
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
def plot_maps_mean(clim_anomaly_pre, lim_anomaly_tas):
		
	fig = plt.figure()

	levs1 = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
	levs2 = [0, 0.5, 1, 1.5, 2, 3, 4, 5]

	ax = fig.add_subplot(2, 2, 1)
	plt.title(u'A) Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, clim_anomaly_pre, levels=levs1, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(2, 2, 2)
	plt.title(u'B) Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, clim_anomaly_pre, levels=levs1, latlon=True, cmap=cm.BrBG)

	ax = fig.add_subplot(2, 2, 3)
	plt.title(u'C) Reg - Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, lim_anomaly_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(2, 2, 4)
	plt.title(u'C) Reg - Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, lim_anomaly_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd) 

	fig.tight_layout()
	
	return plot_maps_mean


# Import regcm exp and cru databases 	   
lat, lon, clim_anomaly_pre = import_obs('pre')
lat, lon, clim_anomaly_tas = import_obs('tmp')

# Plot maps with the function
plt_map = plot_maps_mean(clim_anomaly_pre, clim_anomaly_tas)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.25, hspace=0.10)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_anomaly_cru.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



