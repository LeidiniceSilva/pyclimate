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


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, obs
	
	
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
	
	
def plot_maps_mean(obs_pre, rcm_pre, gcm_pre, diff_pre_rcm_gcm, pre_rcm_obs, pre_gcm_obs, av_pre_rcm_gcm_obs, obs_tas, rcm_tas, gcm_tas, diff_tas_rcm_gcm, tas_rcm_obs, tas_gcm_obs, av_tas_rcm_gcm_obs):
		
	fig = plt.figure(figsize=(6,3))

	levs1 = [1, 2, 4, 6, 8, 10, 15]
	levs2 = [20, 22, 24, 26, 28, 30]
	levs3 = [-4, -2, -1, -.5, 0.5, 1, 2, 4]
	levs4 = [-1, -0.7, -0.4, 0, 0.4, 0.7, 1]

	ax = fig.add_subplot(4, 4, 1)
	plt.title(u'A) CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	ax = fig.add_subplot(4, 4, 2)
	plt.title(u'B) Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	ax = fig.add_subplot(4, 4, 3)
	plt.title(u'C) Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, gcm_pre, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 
		
	ax = fig.add_subplot(4, 4, 5)
	plt.title(u'D) Reg - Had (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, diff_pre_rcm_gcm, levels=levs3, latlon=True, cmap=cm.BrBG)  
	
	ax = fig.add_subplot(4, 4, 6)
	plt.title(u'E) Reg - CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, pre_rcm_obs, levels=levs3, latlon=True, cmap=cm.BrBG)

	ax = fig.add_subplot(4, 4, 7)
	plt.title(u'F) Had - CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, pre_gcm_obs, levels=levs3, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 8)
	plt.title(u'G) AV', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, av_pre_rcm_gcm_obs, levels=levs4, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 

	ax = fig.add_subplot(4, 4, 9)
	plt.title(u'H) CRU (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 10)
	plt.title(u'I) Reg (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_tas[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 11)
	plt.title(u'J) Had (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, gcm_tas, levels=levs2, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(4, 4, 13)
	plt.title(u'K) Reg - Had (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, diff_tas_rcm_gcm, levels=levs3, latlon=True, cmap=cm.bwr)
	
	ax = fig.add_subplot(4, 4, 14)
	plt.title(u'L) Reg - CRU (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, tas_rcm_obs, levels=levs3, latlon=True, cmap=cm.bwr) 

	ax = fig.add_subplot(4, 4, 15)
	plt.title(u'M) Had - CRU (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, tas_gcm_obs, levels=levs3, latlon=True, cmap=cm.bwr) 
	cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6)  
		
	ax = fig.add_subplot(4, 4, 16)
	plt.title(u'N) AV', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, av_tas_rcm_gcm_obs, levels=levs4, latlon=True, cmap=cm.PiYG) 
	cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6)  

	fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
	
	return plot_maps_mean


# Import regcm exp and cru databases 	
lat, lon, obs_pre = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')   
lat, lon, rcm_pre = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, gcm_pre = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, obs_tas = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, rcm_tas = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, gcm_tas = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute regcm - hagdem outut
diff_pre_rcm_gcm = rcm_pre - gcm_pre
diff_tas_rcm_gcm = np.nanmean(rcm_tas, axis=0) - gcm_tas

# Compute and plot bias from regcm exp and cru database
pre_rcm_obs = rcm_pre - obs_pre
pre_gcm_obs = gcm_pre - obs_pre
tas_rcm_obs = np.nanmean(rcm_tas, axis=0) - obs_tas
tas_gcm_obs = gcm_tas - obs_tas

# Compute added value from regcm output
av_pre_rcm_gcm_obs = compute_added_value(gcm_pre, rcm_pre, obs_pre)
av_tas_rcm_gcm_obs = compute_added_value(gcm_tas, np.nanmean(rcm_tas, axis=0), obs_tas)

# Plot maps with the function
plt_map = plot_maps_mean(obs_pre, rcm_pre, gcm_pre, diff_pre_rcm_gcm, pre_rcm_obs, pre_gcm_obs, av_pre_rcm_gcm_obs, obs_tas, rcm_tas, gcm_tas, diff_tas_rcm_gcm, tas_rcm_obs, tas_gcm_obs, av_tas_rcm_gcm_obs)
#~ plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.50, hspace=0.5)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



