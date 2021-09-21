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

# Plot models and obs database 
#~ fig = plt.figure(figsize=(8,4))
#~ levs1 = [0, 2, 4, 6, 8, 10, 12, 14]

#~ ax1 = plt.subplot2grid((2,6), (0,0), colspan=2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_cru, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gpcp, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_era5, levels=levs1, latlon=True, cmap=cm.Blues, extend='max')
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_rcm, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gcm, levels=levs1, latlon=True, cmap=cm.Blues, extend='max')
#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ fig_coord = [0.35,0.08,0.3,0.03]
#~ cbar_ax = fig.add_axes(fig_coord)
#~ cb1 = plt.colorbar(plt_maps_bias, cax=cbar_ax, orientation='horizontal', ticks=levs1)

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_clim_reg_exp2_pre.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()

fig = plt.figure(figsize=(6,4))
levs2 = [18, 20, 22, 24, 26, 28, 30, 32]

ax1 = fig.add_subplot(2, 2, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_cru, levels=levs2, latlon=True, cmap=cm.OrRd)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax2 = fig.add_subplot(2, 2, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_era5, levels=levs2, latlon=True, cmap=cm.OrRd, extend='max')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

ax3 = fig.add_subplot(2, 2, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.OrRd)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax4 = fig.add_subplot(2, 2, 4)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_gcm, levels=levs2, latlon=True, cmap=cm.OrRd, extend='max')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



