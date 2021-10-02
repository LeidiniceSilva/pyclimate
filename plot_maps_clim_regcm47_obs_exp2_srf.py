# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot seasonal climatology maps from regcm47 and hadgem models and obs database"

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
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_obs = value[2:240:3,:,:]	
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)

	return lat, lon, djf_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)

	return lat, lon, djf_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)

	return lat, lon, djf_gcm
	

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
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)

	return map, xx, yy
	

# Import models and obs database 
lat, lon, pre_djf_cru  = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	
lat, lon, pre_djf_gpcp = import_obs('precip', 'amz_neb', 'gpcp_v2.2', '1986-2005')	   
lat, lon, pre_djf_era5 = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, pre_djf_rcm = import_rcm('pr', 'amz_neb', 'historical', '1986-2005')
lat, lon, pre_djf_gcm = import_gcm('pr', 'amz_neb', 'historical', '1986-2005')

lat, lon, tas_djf_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_djf_era5 = import_obs('t2m', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, tas_djf_rcm = import_rcm('tas', 'amz_neb', 'historical', '1986-2005')
lat, lon, tas_djf_gcm = import_gcm('tas', 'amz_neb', 'historical', '1986-2005')

lat, lon, uas_djf_era5 = import_obs('u10', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, ua_djf_rcm = import_rcm('uas', 'amz_neb', 'historical', '1986-2005')
lat, lon, uas_djf_gcm = import_gcm('uas', 'amz_neb', 'historical', '1986-2005')

lat, lon, vas_djf_era5 = import_obs('v10', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, va_djf_rcm = import_rcm('vas', 'amz_neb', 'historical', '1986-2005')
lat, lon, vas_djf_gcm = import_gcm('vas', 'amz_neb', 'historical', '1986-2005')

lat, lon, psl_djf_era5 = import_obs('msl', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, psl_djf_rcm = import_rcm('psl', 'amz_neb', 'historical', '1986-2005')
lat, lon, psl_djf_gcm = import_gcm('psl', 'amz_neb', 'historical', '1986-2005')

uas_djf_rcm = ua_djf_rcm[0,:,:]
vas_djf_rcm = va_djf_rcm[0,:,:]

# Plot models and obs database 
#~ fig = plt.figure(figsize=(8,4))
#~ levs1 = [0, 2, 4, 6, 8, 10, 12, 14]

#~ ax1 = fig.add_subplot(2, 3, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_cru, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax2 = fig.add_subplot(2, 3, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gpcp, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax2 = fig.add_subplot(2, 3, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_era5, levels=levs1, latlon=True, cmap=cm.Blues, extend='max')
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax4 = fig.add_subplot(2, 3, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_rcm, levels=levs1, latlon=True, cmap=cm.Blues)
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax5 = fig.add_subplot(2, 3, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, pre_djf_gcm, levels=levs1, latlon=True, cmap=cm.Blues, extend='max')
#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_clim_reg_exp2_pre.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()

#~ fig = plt.figure(figsize=(5,4))
#~ levs2 = [18, 20, 22, 24, 26, 28, 30, 32]

#~ ax1 = fig.add_subplot(2, 2, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, tas_djf_cru, levels=levs2, latlon=True, cmap=cm.OrRd)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax2 = fig.add_subplot(2, 2, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, tas_djf_era5, levels=levs2, latlon=True, cmap=cm.OrRd, extend='max')
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax3 = fig.add_subplot(2, 2, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, tas_djf_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.OrRd)
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax4 = fig.add_subplot(2, 2, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, tas_djf_gcm, levels=levs2, latlon=True, cmap=cm.OrRd, extend='max')
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_clim_reg_exp2_tas.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()

fig = plt.figure(figsize=(5,4))
levs2 = [1006, 1008, 1010, 1012, 1014, 1016, 1018, 1020]

ax1 = fig.add_subplot(3, 1, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, psl_djf_era5, levels=levs2, latlon=True, cmap=cm.bwr, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], uas_djf_era5[::10,::10], vas_djf_era5[::10,::10]) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
a=map.contour(xx, yy, psl_djf_era5, np.arange(1006, 1020, 2), linewidths=0.5, colors='black')
plt.clabel(a, fmt='%d', fontsize=5, colors='black')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')

ax2 = fig.add_subplot(3, 1, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, psl_djf_rcm, levels=levs2, latlon=True, cmap=cm.bwr, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], uas_djf_rcm[::10,::10], vas_djf_rcm[::10,::10]) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
a=map.contour(xx, yy, psl_djf_rcm, np.arange(1006, 1020, 2), linewidths=0.5, colors='black')
plt.clabel(a, fmt='%d', fontsize=5, colors='black')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')

ax3 = fig.add_subplot(3, 1, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, psl_djf_gcm, levels=levs2, latlon=True, cmap=cm.bwr, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], uas_djf_gcm[::10,::10], vas_djf_gcm[::10,::10]) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.contour(xx, yy, psl_djf_gcm, np.arange(1006, 1020, 2), linewidths=0.5, colors='black')
a=map.contour(xx, yy, psl_djf_gcm, np.arange(1006, 1020, 2), linewidths=0.5, colors='black')
plt.clabel(a, fmt='%d', fontsize=5, colors='black')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_psl.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


