# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot seasonal climatology bias maps from regcm47 and hadgem models and obs database"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib as mpl 
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from comp_statist_indices import compute_av


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_obs = value[2:240:3,:,:]	
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:80:4], axis=0)

	return lat, lon, djf_obs, jja_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)

	return lat, lon, djf_rcm, jja_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp2'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)

	return lat, lon, djf_gcm, jja_gcm
	

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
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	

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
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	

# Import models and obs database 
lat, lon, pre_djf_cru,  pre_jja_cru  = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	
lat, lon, pre_djf_gpcp, pre_jja_gpcp = import_obs('precip', 'amz_neb', 'gpcp_v2.2', '1986-2005')	   
lat, lon, pre_djf_era5, pre_jja_era5 = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, pre_djf_rcm, pre_jja_rcm = import_rcm('pr', 'amz_neb', 'historical', '1986-2005')
lat, lon, pre_djf_gcm, pre_jja_gcm = import_gcm('pr', 'amz_neb', 'historical', '1986-2005')

lat, lon, tas_djf_cru, tas_jja_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_djf_era5, tas_jja_era5 = import_obs('t2m', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, tas_djf_rcm, tas_jja_rcm = import_rcm('tas', 'amz_neb', 'historical', '1986-2005')
lat, lon, tas_djf_gcm, tas_jja_gcm = import_gcm('tas', 'amz_neb', 'historical', '1986-2005')
	
# Compute bias from models and obs database 
pre_djf_rcm_cru = pre_djf_rcm - pre_djf_cru
pre_djf_rcm_gpcp = pre_djf_rcm - pre_djf_gpcp
pre_djf_rcm_era5 = pre_djf_rcm - pre_djf_era5
pre_djf_gcm_cru = pre_djf_gcm - pre_djf_cru
pre_djf_gcm_gpcp = pre_djf_gcm - pre_djf_gpcp
pre_djf_gcm_era5 = pre_djf_gcm - pre_djf_era5

tas_djf_rcm_cru = np.nanmean(tas_djf_rcm, axis=0) - tas_djf_cru
tas_djf_rcm_era5 = np.nanmean(tas_djf_rcm, axis=0) - tas_djf_era5
tas_djf_gcm_cru = tas_djf_gcm - tas_djf_cru
tas_djf_gcm_era5 = tas_djf_gcm - tas_djf_era5

# Plot models and obs database 
#~ fig = plt.figure(figsize=(8,4))
#~ levs1 = [-6, -4, -2, 2, 4, 6]

#~ ax1 = fig.add_subplot(2, 3, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_rcm_cru, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax2 = fig.add_subplot(2, 3, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_rcm_gpcp, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax3 = fig.add_subplot(2, 3, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_rcm_era5, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax4 = fig.add_subplot(2, 3, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gcm_cru, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax4 = fig.add_subplot(2, 3, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gcm_gpcp, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax4 = fig.add_subplot(2, 3, 6)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gcm_era5, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_bias_reg_exp2_pre.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()


# Plot models and obs database 
fig = plt.figure(figsize=(6,4))
levs1 = [-5, -3, -1, 1, 3, 5]

ax1 = fig.add_subplot(2, 2, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_rcm_cru, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax2 = fig.add_subplot(2, 2, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_rcm_era5, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

ax3 = fig.add_subplot(2, 2, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_gcm_cru, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax4 = fig.add_subplot(2, 2, 4)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_gcm_era5, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_exp2_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



