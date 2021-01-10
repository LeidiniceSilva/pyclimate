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

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_rcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_reg_had_{2}_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_Amon_HadGEM2-ES_{2}_r1i1p1_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, gcm
	

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
	
	
def plot_maps_diff(diff_pre_reg_r26_hist, diff_pre_had_r26_hist, diff_pre_reg_r85_hist, diff_pre_had_r85_hist, diff_tas_reg_r26_hist, diff_tas_had_r26_hist, diff_tas_reg_r85_hist, diff_tas_had_r85_hist):
		
	fig = plt.figure()

	levs1 = [-3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3]
	levs2 = [0.1, 0.5, 1, 1.5, 2, 2.5, 3]
	levs3 = [0.5, 1, 2, 3, 4, 5, 6, 7, 8]

	ax = fig.add_subplot(4, 2, 1)
	plt.title(u'A) Reg RCP26-Hist (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_pre_reg_r26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 2, 2)
	plt.title(u'B) Had RCP26-Hist (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_pre_had_r26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 2, 3)
	plt.title(u'C) Reg RCP85-Hist (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_pre_reg_r85_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 2, 4)
	plt.title(u'D) Had RCP85-Hist  (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_pre_had_r85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 2, 5)
	plt.title(u'E) Reg RCP26-Hist (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_tas_reg_r26_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 2, 6)
	plt.title(u'F) Had RCP26-Hist (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_tas_had_r26_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 2, 7)
	plt.title(u'G) Reg RCP85-Hist (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_tas_reg_r85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 2, 8)
	plt.title(u'H) Had RCP85-Hist (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, diff_tas_had_r85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	fig.tight_layout()
	
	return plot_maps_diff


# Import regcm and hadgem exp
lat, lon, pre_reg_hist = import_rcm('pr', 'hist', '1986-2005')
lat, lon, pre_had_hist = import_gcm('pr', 'hist', '1986-2005')
lat, lon, pre_reg_rcp26 = import_rcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_had_rcp26 = import_gcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_reg_rcp85 = import_rcm('pr', 'rcp85', '2080-2099')
lat, lon, pre_had_rcp85 = import_gcm('pr', 'rcp85', '2080-2099')

lat, lon, tas_reg_hist = import_rcm('tas', 'hist', '1986-2005')
lat, lon, tas_had_hist = import_gcm('tas', 'hist', '1986-2005')
lat, lon, tas_reg_rcp26 = import_rcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_had_rcp26 = import_gcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_reg_rcp85 = import_rcm('tas', 'rcp85', '2080-2099')
lat, lon, tas_had_rcp85 = import_gcm('tas', 'rcp85', '2080-2099')

# Compute change from rcp and historical period
diff_pre_reg_rcp26_hist = pre_reg_rcp26 - pre_reg_hist
diff_pre_had_rcp26_hist = pre_had_rcp26 - pre_had_hist
diff_pre_reg_rcp85_hist = pre_reg_rcp85 - pre_reg_hist
diff_pre_had_rcp85_hist = pre_had_rcp85 - pre_had_hist

diff_tas_reg_rcp26_hist = np.nanmean(tas_reg_rcp26, axis=0) - np.nanmean(tas_reg_hist, axis=0)
diff_tas_had_rcp26_hist = tas_had_rcp26 - tas_had_hist
diff_tas_reg_rcp85_hist = np.nanmean(tas_reg_rcp85, axis=0) - np.nanmean(tas_reg_hist, axis=0)
diff_tas_had_rcp85_hist = tas_had_rcp85 - tas_had_hist

# Plot bias maps 
plt_map = plot_maps_diff(diff_pre_reg_rcp26_hist, diff_pre_had_rcp26_hist, diff_pre_reg_rcp85_hist, diff_pre_had_rcp85_hist, diff_tas_reg_rcp26_hist, diff_tas_had_rcp26_hist, diff_tas_reg_rcp85_hist, diff_tas_had_rcp85_hist)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



