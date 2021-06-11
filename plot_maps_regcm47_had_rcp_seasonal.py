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
    
    
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)

	return lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(season_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(season_gcm[2:80:4], axis=0)

	return lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm
	

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
	lats = np.arange(-20.,15.,0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
def plot_maps_diff(pre_djf_rcm_rcp26_hist,pre_mam_rcm_rcp26_hist,pre_jja_rcm_rcp26_hist,pre_son_rcm_rcp26_hist,pre_djf_gcm_rcp26_hist,pre_mam_gcm_rcp26_hist,pre_jja_gcm_rcp26_hist,pre_son_gcm_rcp26_hist,pre_djf_rcm_rcp85_hist,pre_mam_rcm_rcp85_hist,pre_jja_rcm_rcp85_hist,pre_son_rcm_rcp85_hist,pre_djf_gcm_rcp85_hist,pre_mam_gcm_rcp85_hist,pre_jja_gcm_rcp85_hist,pre_son_gcm_rcp85_hist,tas_djf_rcm_rcp26_hist,tas_mam_rcm_rcp26_hist,tas_jja_rcm_rcp26_hist,tas_son_rcm_rcp26_hist,tas_djf_gcm_rcp26_hist,tas_mam_gcm_rcp26_hist,tas_jja_gcm_rcp26_hist,tas_son_gcm_rcp26_hist,tas_djf_rcm_rcp85_hist,tas_mam_rcm_rcp85_hist,tas_jja_rcm_rcp85_hist,tas_son_rcm_rcp85_hist,tas_djf_gcm_rcp85_hist,tas_mam_gcm_rcp85_hist,tas_jja_gcm_rcp85_hist,tas_son_gcm_rcp85_hist):
		
	fig = plt.figure(figsize=(6,3))

	#~ levs1 = [-4, -2, -1, -0.5, 0.5, 1, 2, 4]

	#~ ax = fig.add_subplot(4, 4, 1)
	#~ plt.title(u'A) Reg RCP26-Hist \n DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 2)
	#~ plt.title(u'B) Reg RCP26-Hist \n MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 3)
	#~ plt.title(u'C) Reg RCP26-Hist \n JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 4)
	#~ plt.title(u'D) Reg RCP26-Hist  \n SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	#~ cbar.ax.tick_params(labelsize=6) 
	
	#~ ax = fig.add_subplot(4, 4, 5)
	#~ plt.title(u'E) Had RCP26-Hist \n DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 6)
	#~ plt.title(u'F) Had RCP26-Hist \n MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 7)
	#~ plt.title(u'G) Had RCP26-Hist \n JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 8)
	#~ plt.title(u'H) Had RCP26-Hist \n SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	#~ cbar.ax.tick_params(labelsize=6) 

	#~ ax = fig.add_subplot(4, 4, 9)
	#~ plt.title(u'I) Reg RCP26-Hist \n DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 10)
	#~ plt.title(u'J) Reg RCP85-Hist \n MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 11)
	#~ plt.title(u'K) Reg RCP85-Hist \n JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 12)
	#~ plt.title(u'L) Reg RCP85-Hist  \n SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	#~ cbar.ax.tick_params(labelsize=6) 
	
	#~ ax = fig.add_subplot(4, 4, 13)
	#~ plt.title(u'M) Had RCP85-Hist \n DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 14)
	#~ plt.title(u'N) Had RCP85-Hist \n MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	
	#~ ax = fig.add_subplot(4, 4, 15)
	#~ plt.title(u'O) Had RCP85-Hist \n JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
	
	#~ ax = fig.add_subplot(4, 4, 16)
	#~ plt.title(u'P) Had RCP85-Hist \n SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	#~ cbar.ax.tick_params(labelsize=6) 

	levs1 = [0.5, 1, 1.5, 2, 3, 4, 5, 7, 9]

	ax = fig.add_subplot(4, 4, 1)
	plt.title(u'A) Reg RCP26-Hist \n DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 2)
	plt.title(u'B) Reg RCP26-Hist \n MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 3)
	plt.title(u'C) Reg RCP26-Hist \n JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 4)
	plt.title(u'D) Reg RCP26-Hist \n SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 5)
	plt.title(u'E) Had RCP26-Hist \n DJF (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 6)
	plt.title(u'F) Had RCP26-Hist \n MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 7)
	plt.title(u'G) Had RCP26-Hist \n JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 8)
	plt.title(u'H) Had RCP26-Hist \n SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 

	ax = fig.add_subplot(4, 4, 9)
	plt.title(u'I) Reg RCP26-Hist \n DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 10)
	plt.title(u'J) Reg RCP85-Hist \n MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 11)
	plt.title(u'K) Reg RCP85-Hist \n JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 12)
	plt.title(u'L) Reg RCP85-Hist  \n SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 13)
	plt.title(u'M) Had RCP85-Hist \n DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 14)
	plt.title(u'N) Had RCP85-Hist \n MAM (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_mam_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
	ax = fig.add_subplot(4, 4, 15)
	plt.title(u'O) Had RCP85-Hist \n JJA (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	
	ax = fig.add_subplot(4, 4, 16)
	plt.title(u'P) Had RCP85-Hist \n SON (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6) 

	fig.tight_layout()
		
	return plot_maps_diff


# Import regcm and hadgem exp
lat, lon, pre_djf_rcm_hist, pre_mam_rcm_hist, pre_jja_rcm_hist, pre_son_rcm_hist = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_djf_gcm_hist, pre_mam_gcm_hist, pre_jja_gcm_hist, pre_son_gcm_hist = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_djf_rcm_rcp26, pre_mam_rcm_rcp26, pre_jja_rcm_rcp26, pre_son_rcm_rcp26 = import_rcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_djf_gcm_rcp26, pre_mam_gcm_rcp26, pre_jja_gcm_rcp26, pre_son_gcm_rcp26 = import_gcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_djf_rcm_rcp85, pre_mam_rcm_rcp85, pre_jja_rcm_rcp85, pre_son_rcm_rcp85 = import_rcm('pr', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, pre_djf_gcm_rcp85, pre_mam_gcm_rcp85, pre_jja_gcm_rcp85, pre_son_gcm_rcp85 = import_gcm('pr', 'amz_neb', 'rcp85', '2080-2099')

lat, lon, tas_djf_rcm_hist, tas_mam_rcm_hist, tas_jja_rcm_hist, tas_son_rcm_hist = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_djf_gcm_hist, tas_mam_gcm_hist, tas_jja_gcm_hist, tas_son_gcm_hist = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_djf_rcm_rcp26, tas_mam_rcm_rcp26, tas_jja_rcm_rcp26, tas_son_rcm_rcp26 = import_rcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_djf_gcm_rcp26, tas_mam_gcm_rcp26, tas_jja_gcm_rcp26, tas_son_gcm_rcp26 = import_gcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_djf_rcm_rcp85, tas_mam_rcm_rcp85, tas_jja_rcm_rcp85, tas_son_rcm_rcp85 = import_rcm('tas', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, tas_djf_gcm_rcp85, tas_mam_gcm_rcp85, tas_jja_gcm_rcp85, tas_son_gcm_rcp85 = import_gcm('tas', 'amz_neb', 'rcp85', '2080-2099')

# Compute change from rcp and historical period
pre_djf_rcm_rcp26_hist = pre_djf_rcm_rcp26 - pre_djf_rcm_hist
pre_mam_rcm_rcp26_hist = pre_mam_rcm_rcp26 - pre_mam_rcm_hist
pre_jja_rcm_rcp26_hist = pre_jja_rcm_rcp26 - pre_jja_rcm_hist
pre_son_rcm_rcp26_hist = pre_son_rcm_rcp26 - pre_son_rcm_hist
pre_djf_rcm_rcp85_hist = pre_djf_rcm_rcp85 - pre_djf_rcm_hist
pre_mam_rcm_rcp85_hist = pre_mam_rcm_rcp85 - pre_mam_rcm_hist
pre_jja_rcm_rcp85_hist = pre_jja_rcm_rcp85 - pre_jja_rcm_hist
pre_son_rcm_rcp85_hist = pre_son_rcm_rcp85 - pre_son_rcm_hist
pre_djf_gcm_rcp26_hist = pre_djf_gcm_rcp26 - pre_djf_gcm_hist
pre_mam_gcm_rcp26_hist = pre_mam_gcm_rcp26 - pre_mam_gcm_hist
pre_jja_gcm_rcp26_hist = pre_jja_gcm_rcp26 - pre_jja_gcm_hist
pre_son_gcm_rcp26_hist = pre_son_gcm_rcp26 - pre_son_gcm_hist
pre_djf_gcm_rcp85_hist = pre_djf_gcm_rcp85 - pre_djf_gcm_hist
pre_mam_gcm_rcp85_hist = pre_mam_gcm_rcp85 - pre_mam_gcm_hist
pre_jja_gcm_rcp85_hist = pre_jja_gcm_rcp85 - pre_jja_gcm_hist
pre_son_gcm_rcp85_hist = pre_son_gcm_rcp85 - pre_son_gcm_hist

tas_djf_rcm_rcp26_hist = np.nanmean(tas_djf_rcm_rcp26, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp26_hist = np.nanmean(tas_mam_rcm_rcp26, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp26_hist = np.nanmean(tas_jja_rcm_rcp26, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp26_hist = np.nanmean(tas_son_rcm_rcp26, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_djf_rcm_rcp85_hist = np.nanmean(tas_djf_rcm_rcp85, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp85_hist = np.nanmean(tas_mam_rcm_rcp85, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp85_hist = np.nanmean(tas_jja_rcm_rcp85, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp85_hist = np.nanmean(tas_son_rcm_rcp85, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_djf_gcm_rcp26_hist = tas_djf_gcm_rcp26 - tas_djf_gcm_hist
tas_mam_gcm_rcp26_hist = tas_mam_gcm_rcp26 - tas_mam_gcm_hist
tas_jja_gcm_rcp26_hist = tas_jja_gcm_rcp26 - tas_jja_gcm_hist
tas_son_gcm_rcp26_hist = tas_son_gcm_rcp26 - tas_son_gcm_hist
tas_djf_gcm_rcp85_hist = tas_djf_gcm_rcp85 - tas_djf_gcm_hist
tas_mam_gcm_rcp85_hist = tas_mam_gcm_rcp85 - tas_mam_gcm_hist
tas_jja_gcm_rcp85_hist = tas_jja_gcm_rcp85 - tas_jja_gcm_hist
tas_son_gcm_rcp85_hist = tas_son_gcm_rcp85 - tas_son_gcm_hist

# Plot bias maps 
plt_map = plot_maps_diff(pre_djf_rcm_rcp26_hist,pre_mam_rcm_rcp26_hist,pre_jja_rcm_rcp26_hist,pre_son_rcm_rcp26_hist,
pre_djf_gcm_rcp26_hist,pre_mam_gcm_rcp26_hist,pre_jja_gcm_rcp26_hist,pre_son_gcm_rcp26_hist,
pre_djf_rcm_rcp85_hist,pre_mam_rcm_rcp85_hist,pre_jja_rcm_rcp85_hist,pre_son_rcm_rcp85_hist,
pre_djf_gcm_rcp85_hist,pre_mam_gcm_rcp85_hist,pre_jja_gcm_rcp85_hist,pre_son_gcm_rcp85_hist,
tas_djf_rcm_rcp26_hist,tas_mam_rcm_rcp26_hist,tas_jja_rcm_rcp26_hist,tas_son_rcm_rcp26_hist,
tas_djf_gcm_rcp26_hist,tas_mam_gcm_rcp26_hist,tas_jja_gcm_rcp26_hist,tas_son_gcm_rcp26_hist,
tas_djf_rcm_rcp85_hist,tas_mam_rcm_rcp85_hist,tas_jja_rcm_rcp85_hist,tas_son_rcm_rcp85_hist,
tas_djf_gcm_rcp85_hist,tas_mam_gcm_rcp85_hist,tas_jja_gcm_rcp85_hist,tas_son_gcm_rcp85_hist)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.95)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp-hist_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()









