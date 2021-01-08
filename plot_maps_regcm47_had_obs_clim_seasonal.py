# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal climatology maps from Reg and Had models end obs database"

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


def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]
	#~ djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	#~ mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	#~ jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	#~ son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)

	djf_rcm = np.nanmean(np.nanmean(season_rcm[3:80:4], axis=0), axis=0)
	mam_rcm = np.nanmean(np.nanmean(season_rcm[0:80:4], axis=0), axis=0)
	jja_rcm = np.nanmean(np.nanmean(season_rcm[1:80:4], axis=0), axis=0)
	son_rcm = np.nanmean(np.nanmean(season_rcm[2:80:4], axis=0), axis=0)

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

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_obs = value[2:240:3,:,:]
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)
	mam_obs = np.nanmean(season_obs[0:80:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:80:4], axis=0)
	son_obs = np.nanmean(season_obs[2:80:4], axis=0)

	return lat, lon, djf_obs, mam_obs, jja_obs, son_obs
	

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
	
	
#~ def plot_maps_seasonal(djf_rcm,mam_rcm,jja_rcm,son_rcm,djf_gcm,mam_gcm,jja_gcm,son_gcm,djf_cru,mam_cru,jja_cru,son_cru,djf_udel,mam_udel,jja_udel,son_udel,djf_chirps,mam_chirps,jja_chirps,son_chirps,djf_era5,mam_era5,jja_era5,son_era5):
def plot_maps_seasonal(djf_rcm,mam_rcm,jja_rcm,son_rcm,djf_gcm,mam_gcm,jja_gcm,son_gcm,djf_cru,mam_cru,jja_cru,son_cru,djf_udel,mam_udel,jja_udel,son_udel,djf_era5,mam_era5,jja_era5,son_era5):

	fig = plt.figure()
	
	#~ levs1 = [1, 2, 4, 6, 8, 10, 15]
	#~ ax = fig.add_subplot(6, 4, 1)
	#~ plt.title(u'A) Reg DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 2)
	#~ plt.title(u'B) Reg MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 3)
	#~ plt.title(u'C) Reg JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 4)
	#~ plt.title(u'D) Reg SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6) 

	#~ ax = fig.add_subplot(6, 4, 5)
	#~ plt.title(u'E) Had DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 6)
	#~ plt.title(u'F) Had MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 7)
	#~ plt.title(u'G) Had JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 8)
	#~ plt.title(u'H) Had SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6) 

	#~ ax = fig.add_subplot(6, 4, 9)
	#~ plt.title(u'I) CRU DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 10)
	#~ plt.title(u'J) CRU MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 11)
	#~ plt.title(u'K) CRU JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 12)
	#~ plt.title(u'L) CRU SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6)

	#~ ax = fig.add_subplot(6, 4, 13)
	#~ plt.title(u'M) UDEL DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_udel, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 14)
	#~ plt.title(u'N) UDEL MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_udel, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 15)
	#~ plt.title(u'O) UDEL JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_udel, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 16)
	#~ plt.title(u'P) UDEL SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_udel, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6)	

	#~ ax = fig.add_subplot(6, 4, 17)
	#~ plt.title(u'Q) CHIRPS DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_chirps, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 18)
	#~ plt.title(u'R) CHIRPS MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_chirps, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 19)
	#~ plt.title(u'S) CHIRPS JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_chirps, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 20)
	#~ plt.title(u'T) CHIRPS SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_chirps, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6)	

	#~ ax = fig.add_subplot(6, 4, 21)
	#~ plt.title(u'U) ERA5 DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, djf_era5, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 22)
	#~ plt.title(u'V) ERA5 MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, mam_era5, levels=levs1, latlon=True, cmap=cm.YlGnBu)

	#~ ax = fig.add_subplot(6, 4, 23)
	#~ plt.title(u'W) ERA5 JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, jja_era5, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

	#~ ax = fig.add_subplot(6, 4, 24)
	#~ plt.title(u'X) ERA5 SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	#~ plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	#~ map, xx, yy = basemap(lat, lon)
	#~ plot_maps_mean = map.contourf(xx, yy, son_era5, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=6)	
	
	levs1 = [22, 24, 26, 28, 30, 32]
	ax = fig.add_subplot(5, 4, 1)
	plt.title(u'A) Reg DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, djf_rcm, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 2)
	plt.title(u'B) Reg MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, mam_rcm, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 3)
	plt.title(u'C) Reg JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, jja_rcm, levels=levs1, latlon=True, cmap=cm.YlOrRd) 

	ax = fig.add_subplot(5, 4, 4)
	plt.title(u'D) Reg SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, son_rcm, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	ax = fig.add_subplot(5, 4, 5)
	plt.title(u'E) Had DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, djf_gcm, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 6)
	plt.title(u'F) Had MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, mam_gcm, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 7)
	plt.title(u'G) Had JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, jja_gcm, levels=levs1, latlon=True, cmap=cm.YlOrRd) 

	ax = fig.add_subplot(5, 4, 8)
	plt.title(u'H) Had SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, son_gcm, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	ax = fig.add_subplot(5, 4, 9)
	plt.title(u'I) CRU DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, djf_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 10)
	plt.title(u'J) CRU MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, mam_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 11)
	plt.title(u'K) CRU JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, jja_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd) 

	ax = fig.add_subplot(5, 4, 12)
	plt.title(u'L) CRU SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, son_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)

	ax = fig.add_subplot(5, 4, 13)
	plt.title(u'M) UDEL DJF (°C)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, djf_udel, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 14)
	plt.title(u'N) UDEL MAM (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, mam_udel, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 15)
	plt.title(u'O) UDEL JJA (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, jja_udel, levels=levs1, latlon=True, cmap=cm.YlOrRd) 

	ax = fig.add_subplot(5, 4, 16)
	plt.title(u'P) UDEL SON (°C)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, son_udel, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)	

	ax = fig.add_subplot(5, 4, 17)
	plt.title(u'Q) ERA5 DJF (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, djf_era5, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 18)
	plt.title(u'R) ERA5 MAM (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, mam_era5, levels=levs1, latlon=True, cmap=cm.YlOrRd)

	ax = fig.add_subplot(5, 4, 19)
	plt.title(u'S) ERA5 JJA (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, jja_era5, levels=levs1, latlon=True, cmap=cm.YlOrRd) 

	ax = fig.add_subplot(5, 4, 20)
	plt.title(u'T) ERA5 SON (°C)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, son_era5, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)	

	fig.tight_layout()
	
	return plot_maps_mean


# Import regcm exp and cru databases 	   
#~ lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
#~ lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
#~ lat, lon, djf_cru, mam_cru, jja_cru, son_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')
#~ lat, lon, djf_udel, mam_udel, jja_udel, son_udel = import_obs('pre', 'amz_neb', 'udel_v301', '1986-2005')
#~ lat, lon, djf_chirps, mam_chirps, jja_chirps, son_chirps = import_obs('precip', 'amz_neb', 'chirps-v2.0', '1986-2005')
#~ lat, lon, djf_era5, mam_era5, jja_era5, son_era5 = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')

lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, djf_cru, mam_cru, jja_cru, son_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, djf_udel, mam_udel, jja_udel, son_udel = import_obs('temp', 'amz_neb', 'udel_v301', '1986-2005')
lat, lon, djf_era5, mam_era5, jja_era5, son_era5 = import_obs('t2m', 'amz_neb', 'era5', '1986-2005')

# Plot maps with the function
#~ plt_map = plot_maps_seasonal(djf_rcm,mam_rcm,jja_rcm,son_rcm,djf_gcm,mam_gcm,jja_gcm,son_gcm,djf_cru,mam_cru,jja_cru,son_cru,djf_udel,mam_udel,jja_udel,son_udel,djf_chirps,mam_chirps,jja_chirps,son_chirps,djf_era5,mam_era5,jja_era5,son_era5)
plt_map = plot_maps_seasonal(djf_rcm,mam_rcm,jja_rcm,son_rcm,djf_gcm,mam_gcm,jja_gcm,son_gcm,djf_cru,mam_cru,jja_cru,son_cru,djf_udel,mam_udel,jja_udel,son_udel,djf_era5,mam_era5,jja_era5,son_era5)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_seasonal_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






