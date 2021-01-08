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

	
def import_obs(var, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_amz_neb_{2}_obs_mon_{3}_lonlat.nc'.format(path, var, dataset, dt)	
			
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
	
	
def plot_maps_bias(p_reg_cru, p_reg_udel, p_reg_chirps, p_reg_era5, p_had_cru, p_had_udel, p_had_chirps, p_had_era5, t_reg_cru, t_reg_udel, t_reg_era5, t_had_cru, t_had_udel, t_had_era5):
		
	fig = plt.figure()

	levs = [-4, -3, -2, -1, 1, 2, 3, 4]
	
	ax = fig.add_subplot(4, 4, 1)
	plt.title(u'A) Reg - CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_reg_cru, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 2)
	plt.title(u'B) Reg - UDEL (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_reg_udel, levels=levs, latlon=True, cmap=cm.BrBG)

	ax = fig.add_subplot(4, 4, 3)
	plt.title(u'C) Reg - CHIRPS (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_reg_chirps, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 4)
	plt.title(u'D) Reg - ERA5 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_reg_era5, levels=levs, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 5)
	plt.title(u'E) Had - CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_had_cru, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 6)
	plt.title(u'F) Had - UDEL (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_had_udel, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 7)
	plt.title(u'G) Had - CHIRPS (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_had_chirps, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 8)
	plt.title(u'H) Had - ERA5 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, p_had_era5, levels=levs, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 9)
	plt.title(u'I) Reg - CRU (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_reg_cru[0,:,:], levels=levs, latlon=True, cmap=cm.bwr)
	
	ax = fig.add_subplot(4, 4, 10)
	plt.title(u'J) Reg - UDEL (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_reg_udel[0,:,:], levels=levs, latlon=True, cmap=cm.bwr)
	
	ax = fig.add_subplot(4, 4, 11)
	plt.title(u'L) Reg - ERA5 (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_reg_era5[0,:,:], levels=levs, latlon=True, cmap=cm.bwr) 
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(4, 4, 13)
	plt.title(u'M) Had - CRU (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_had_cru, levels=levs, latlon=True, cmap=cm.bwr) 
	
	ax = fig.add_subplot(4, 4, 14)
	plt.title(u'N) Had - UDEL (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_had_udel, levels=levs, latlon=True, cmap=cm.bwr) 
	
	ax = fig.add_subplot(4, 4, 15)
	plt.title(u'O) Had - ERA5 (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, t_had_era5, levels=levs, latlon=True, cmap=cm.bwr) 
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	fig.tight_layout()
	
	return plt_maps_bias


# Import regcm exp and cru databases 	   
lat, lon, reg_pre = import_rcm('pr', 'hist', '1986-2005')
lat, lon, had_pre = import_gcm('pr', 'hist', '1986-2005')
lat, lon, cru_pre = import_obs('pre', 'cru_ts4.04', '1986-2005')
lat, lon, udel_pre = import_obs('pre', 'udel_v301', '1986-2005')
lat, lon, chirps_pre = import_obs('precip', 'chirps-v2.0', '1986-2005')
lat, lon, era5_pre = import_obs('mtpr', 'era5', '1986-2005')

lat, lon, reg_tas = import_rcm('tas', 'hist', '1986-2005')
lat, lon, had_tas = import_gcm('tas', 'hist', '1986-2005')
lat, lon, cru_tas = import_obs('tmp', 'cru_ts4.04', '1986-2005')
lat, lon, udel_tas = import_obs('temp', 'udel_v301', '1986-2005')
lat, lon, era5_tas = import_obs('t2m', 'era5', '1986-2005')

# Compute and plot bias from regcm exp and cru database
p_reg_cru = reg_pre - cru_pre
p_reg_udel = reg_pre - udel_pre
p_reg_chirps = reg_pre - chirps_pre
p_reg_era5 = reg_pre - era5_pre
p_had_cru = had_pre - cru_pre
p_had_udel = had_pre - udel_pre
p_had_chirps = had_pre - chirps_pre
p_had_era5 = had_pre - era5_pre

t_reg_cru = reg_tas - cru_tas
t_reg_udel = reg_tas - udel_tas
t_reg_era5 = reg_tas - era5_tas
t_had_cru = had_tas- cru_tas
t_had_udel = had_tas- udel_tas
t_had_era5 = had_tas - era5_tas

# Plot maps with the function
plt_map = plot_maps_bias(p_reg_cru, p_reg_udel, p_reg_chirps, p_reg_era5, p_had_cru, p_had_udel, p_had_chirps, p_had_era5, t_reg_cru, t_reg_udel, t_reg_era5, t_had_cru, t_had_udel, t_had_era5)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.25, hspace=0.10)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



