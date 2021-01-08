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

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, gcm

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
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
	
	
def plot_maps_mean(rcm_pre, gcm_pre, diff_rcm_gcm_pre, cru_pre, udel_pre, chirps_pre, era5_pre, rcm_tas, gcm_tas, diff_rcm_gcm_tas, cru_tas, udel_tas, era5_tas):
		
	fig = plt.figure()

	levs1 = [1, 2, 4, 6, 8, 10, 15]
	levs2 = [20, 22, 24, 26, 28, 30]
	levs3 = [-4, -3, -2, -1, 1, 2, 3, 4]

	ax = fig.add_subplot(4, 4, 1)
	plt.title(u'A) Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	
	ax = fig.add_subplot(4, 4, 2)
	plt.title(u'B) Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, gcm_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	ax = fig.add_subplot(4, 4, 3)
	plt.title(u'C) Reg - Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, diff_rcm_gcm_pre, levels=levs3, latlon=True, cmap=cm.BrBG) 
	cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 5)
	plt.title(u'E) CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, cru_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	
	ax = fig.add_subplot(4, 4, 6)
	plt.title(u'F) UDEL (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, udel_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	
	ax = fig.add_subplot(4, 4, 7)
	plt.title(u'G) CHIRPS (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, chirps_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	
	ax = fig.add_subplot(4, 4, 8)
	plt.title(u'H) ERA5 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, era5_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 9)
	plt.title(u'I) Reg (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_tas[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 10)
	plt.title(u'J) Had (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, gcm_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 11)
	plt.title(u'L) Reg - Had (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, diff_rcm_gcm_tas[0,:,:], levels=levs3, latlon=True, cmap=cm.bwr)
	cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(4, 4, 13)
	plt.title(u'N) CRU (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, cru_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 14)
	plt.title(u'O) UDEL (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, udel_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 15)
	plt.title(u'P) ERA5 (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, era5_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	fig.tight_layout()
	
	return plot_maps_mean


# Import regcm exp and cru databases 	   
lat, lon, rcm_pre = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, gcm_pre = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, cru_pre = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, udel_pre = import_obs('pre', 'amz_neb', 'udel_v301', '1986-2005')
lat, lon, chirps_pre = import_obs('precip', 'amz_neb', 'chirps-v2.0', '1986-2005')
lat, lon, era5_pre = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')

lat, lon, rcm_tas = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, gcm_tas = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, cru_tas = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, udel_tas = import_obs('temp', 'amz_neb', 'udel_v301', '1986-2005')
lat, lon, era5_tas = import_obs('t2m', 'amz_neb', 'era5', '1986-2005')

# Compute regcm - hagdem outut
diff_rcm_gcm_pre = rcm_pre - gcm_pre
diff_rcm_gcm_tas = rcm_tas - gcm_tas

# Plot maps with the function
plt_map = plot_maps_mean(rcm_pre, gcm_pre, diff_rcm_gcm_pre, cru_pre, udel_pre, chirps_pre, era5_pre, rcm_tas, gcm_tas, diff_rcm_gcm_tas, cru_tas, udel_tas, era5_tas)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.25, hspace=0.10)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



