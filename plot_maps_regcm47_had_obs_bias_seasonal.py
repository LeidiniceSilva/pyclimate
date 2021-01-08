# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal bias maps from Reg and Had models output"

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
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)

	#~ djf_rcm = np.nanmean(np.nanmean(season_rcm[3:80:4], axis=0), axis=0)
	#~ mam_rcm = np.nanmean(np.nanmean(season_rcm[0:80:4], axis=0), axis=0)
	#~ jja_rcm = np.nanmean(np.nanmean(season_rcm[1:80:4], axis=0), axis=0)
	#~ son_rcm = np.nanmean(np.nanmean(season_rcm[2:80:4], axis=0), axis=0)

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
	
	
def plot_maps_bias_seasonal(djf_b1,djf_b2,djf_b3,djf_b4,mam_b1,mam_b2,mam_b3,mam_b4,jja_b1,jja_b2,jja_b3,jja_b4,son_b1,son_b2,son_b3,son_b4):
		
	fig = plt.figure()

	levs = [-4, -3, -2, -1, 1, 2, 3, 4]
	
	ax = fig.add_subplot(4, 4, 1)
	plt.title(u'A) Had - CRU DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, djf_b1, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 2)
	plt.title(u'B) Had - CRU MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mam_b1, levels=levs, latlon=True, cmap=cm.BrBG)

	ax = fig.add_subplot(4, 4, 3)
	plt.title(u'C) Had - CRU JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, jja_b1, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 4)
	plt.title(u'D) Had - CRU SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, son_b1, levels=levs, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 5)
	plt.title(u'E) Had - UDEL DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, djf_b2, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 6)
	plt.title(u'F) Had - UDEL MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mam_b2, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 7)
	plt.title(u'G) Had - UDEL JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, jja_b2, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 8)
	plt.title(u'H) Had - UDEL SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, son_b2, levels=levs, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(4, 4, 9)
	plt.title(u'I) Had - CHIRPS DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, djf_b3, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 10)
	plt.title(u'J) Had - CHIRPS MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mam_b3, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 11)
	plt.title(u'L) Had - CHIRPS JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, jja_b3, levels=levs, latlon=True, cmap=cm.BrBG)
	
	ax = fig.add_subplot(4, 4, 12)
	plt.title(u'M) Had - CHIRPS SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, son_b3, levels=levs, latlon=True, cmap=cm.BrBG) 
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(4, 4, 13)
	plt.title(u'N) Had - ERA5 DJF (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, djf_b4, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 14)
	plt.title(u'O) Had - ERA5 MAM (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, mam_b4, levels=levs, latlon=True, cmap=cm.BrBG) 
	
	ax = fig.add_subplot(4, 4, 15)
	plt.title(u'P) Had - ERA5 JJA (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, jja_b4, levels=levs, latlon=True, cmap=cm.BrBG) 

	ax = fig.add_subplot(4, 4, 16)
	plt.title(u'Q) Had - ERA5 SON (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)	
	plt_maps_bias = map.contourf(xx, yy, son_b4, levels=levs, latlon=True, cmap=cm.BrBG) 
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	fig.tight_layout()
	
	return plt_maps_bias


# Import regcm exp and cru databases 	   
lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, djf_cru, mam_cru, jja_cru, son_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, djf_udel, mam_udel, jja_udel, son_udel = import_obs('pre', 'amz_neb', 'udel_v301', '1986-2005')
lat, lon, djf_chirps, mam_chirps, jja_chirps, son_chirps = import_obs('precip', 'amz_neb', 'chirps-v2.0', '1986-2005')
lat, lon, djf_era5, mam_era5, jja_era5, son_era5 = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')

# Compute and plot bias from regcm exp and cru database
#~ djf_b1 = djf_rcm - djf_cru
#~ djf_b2 = djf_rcm - djf_udel
#~ djf_b3 = djf_rcm - djf_chirps
#~ djf_b4 = djf_rcm - djf_era5

#~ mam_b1 = mam_rcm - mam_cru
#~ mam_b2 = mam_rcm - mam_udel
#~ mam_b3 = mam_rcm - mam_chirps
#~ mam_b4 = mam_rcm - mam_era5

#~ jja_b1 = jja_rcm - jja_cru
#~ jja_b2 = jja_rcm - jja_udel
#~ jja_b3 = jja_rcm - jja_chirps
#~ jja_b4 = jja_rcm - jja_era5

#~ son_b1 = son_rcm - son_cru
#~ son_b2 = son_rcm - son_udel
#~ son_b3 = son_rcm - son_chirps
#~ son_b4 = son_rcm - son_era5

djf_b1 = djf_rcm - djf_cru
djf_b2 = djf_rcm - djf_udel
djf_b3 = djf_rcm - djf_chirps
djf_b4 = djf_rcm - djf_era5

mam_b1 = mam_rcm - mam_cru
mam_b2 = mam_rcm - mam_udel
mam_b3 = mam_rcm - mam_chirps
mam_b4 = mam_rcm - mam_era5

jja_b1 = jja_rcm - jja_cru
jja_b2 = jja_rcm - jja_udel
jja_b3 = jja_rcm - jja_chirps
jja_b4 = jja_rcm - jja_era5

son_b1 = son_rcm - son_cru
son_b2 = son_rcm - son_udel
son_b3 = son_rcm - son_chirps
son_b4 = son_rcm - son_era5

# Plot maps with the function
plt_map = plot_maps_bias_seasonal(djf_b1,djf_b2,djf_b3,djf_b4,mam_b1,mam_b2,mam_b3,mam_b4,jja_b1,jja_b2,jja_b3,jja_b4,son_b1,son_b2,son_b3,son_b4)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.50)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_seasonal_pre_gcm_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



