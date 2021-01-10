# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot bias maps from Reg and Had models output"

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
	rcm  = np.nanmean(var[:][:-1,:,:], axis=0)

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
	obs  = np.nanmean(var[:][:-1,:,:], axis=0)
	
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
	
	
def plot_maps_av(av_rcm_gcm_cru_pre, av_rcm_gcm_udel_pre, av_rcm_gcm_chirps_pre, av_rcm_gcm_era5_pre, av_rcm_gcm_cru_tas, av_rcm_gcm_udel_tas, av_rcm_gcm_era5_tas):
		
	fig = plt.figure()

	levs = [-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1]

	ax = fig.add_subplot(4, 2, 1)
	plt.title(u'A) AV Reg Had CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_cru_pre, levels=levs, latlon=True, cmap=cm.PiYG)
	
	ax = fig.add_subplot(4, 2, 2)
	plt.title(u'B) AV Reg Had UDEL (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_udel_pre, levels=levs, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 2, 3)
	plt.title(u'C) AV Reg Had CHIRPS (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_chirps_pre, levels=levs, latlon=True, cmap=cm.PiYG)
	
	ax = fig.add_subplot(4, 2, 4)
	plt.title(u'D) AV Reg Had ERA5 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_era5_pre, levels=levs, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 2, 5)
	plt.title(u'E) AV Reg Had CRU (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_cru_tas, levels=levs, latlon=True, cmap=cm.PiYG)
	
	ax = fig.add_subplot(4, 2, 6)
	plt.title(u'F) AV Reg Had UDEL (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_udel_tas, levels=levs, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
		
	ax = fig.add_subplot(4, 2, 7)
	plt.title(u'G) AV Reg Had ERA5 (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', labelpad=15, fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_av = map.contourf(xx, yy, av_rcm_gcm_era5_tas, levels=levs, latlon=True, cmap=cm.PiYG)

	fig.tight_layout()

	return plt_maps_av


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

# Compute added value from regcm output
av_rcm_gcm_cru_pre = compute_added_value(rcm_pre, gcm_pre, cru_pre)
av_rcm_gcm_udel_pre = compute_added_value(rcm_pre, gcm_pre, udel_pre)
av_rcm_gcm_chirps_pre = compute_added_value(rcm_pre, gcm_pre, chirps_pre)
av_rcm_gcm_era5_pre = compute_added_value(rcm_pre, gcm_pre, era5_pre)

av_rcm_gcm_cru_tas = compute_added_value(np.nanmean(rcm_tas, axis=0), gcm_tas, cru_tas)
av_rcm_gcm_udel_tas = compute_added_value(np.nanmean(rcm_tas, axis=0), gcm_tas, udel_tas)
av_rcm_gcm_era5_tas = compute_added_value(np.nanmean(rcm_tas, axis=0), gcm_tas, era5_tas)

# Plot maps with the function
plt_map = plot_maps_av(av_rcm_gcm_cru_pre, av_rcm_gcm_udel_pre, av_rcm_gcm_chirps_pre, av_rcm_gcm_era5_pre, av_rcm_gcm_cru_tas, av_rcm_gcm_udel_tas, av_rcm_gcm_era5_tas)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_av_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



