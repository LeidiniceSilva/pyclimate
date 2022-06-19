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
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	dec = np.nanmean(value[11:240:11], axis=0)
	jan = np.nanmean(value[0:240:11], axis=0)
	feb = np.nanmean(value[1:240:11], axis=0)
	djf = (dec + jan + feb)/3

	return lat, lon, djf
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	dec = np.nanmean(value[11:240:11], axis=0)
	jan = np.nanmean(value[0:240:11], axis=0)
	feb = np.nanmean(value[1:240:11], axis=0)
	djf = (dec + jan + feb)/3

	return lat, lon, djf


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	dec = np.nanmean(value[11:240:11], axis=0)
	jan = np.nanmean(value[0:240:11], axis=0)
	feb = np.nanmean(value[1:240:11], axis=0)
	djf = (dec + jan + feb)/3

	return lat, lon, djf
	

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
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)

	x1,i1 = map(-72,-12)
	x2,i2 = map(-72,0)
	x3,i3 = map(-55,0)
	x4,i4 = map(-55,-12)

	poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
	plt.gca().add_patch(poly1)

	y1,j1 = map(-47,-18)
	y2,j2 = map(-47,-2)
	y3,j3 = map(-35,-2)
	y4,j4 = map(-35,-18)

	poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.)
	plt.gca().add_patch(poly2)
	
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
fig = plt.figure()
levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]

ax1 = fig.add_subplot(3, 2, 1)
map, xx, yy = basemap(lat, lon)
pltfig=map.contourf(xx, yy, pre_djf_rcm_cru, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.92, 0.24, 0.02, 0.5]))
cbar.ax.tick_params(labelsize=8)

ax2 = fig.add_subplot(3, 2, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, pre_djf_gcm_cru, levels=levs1, latlon=True, cmap=cm.BrBG)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax3 = fig.add_subplot(3, 2, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, pre_djf_rcm_gpcp, levels=levs1, latlon=True, cmap=cm.BrBG)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax4 = fig.add_subplot(3, 2, 4)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, pre_djf_gcm_gpcp, levels=levs1, latlon=True, cmap=cm.BrBG)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax4 = fig.add_subplot(3, 2, 5)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, pre_djf_rcm_era5, levels=levs1, latlon=True, cmap=cm.BrBG)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax4 = fig.add_subplot(3, 2, 6)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, pre_djf_gcm_era5, levels=levs1, latlon=True, cmap=cm.BrBG)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_exp2_pre.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

fig = plt.figure(figsize=(8,4))
levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]

ax1 = fig.add_subplot(2, 2, 1)
map, xx, yy = basemap(lat, lon)
pltfig=map.contourf(xx, yy, tas_djf_rcm_cru, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.92, 0.24, 0.02, 0.5]))
cbar.ax.tick_params(labelsize=8)

ax2 = fig.add_subplot(2, 2, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_gcm_cru, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax3 = fig.add_subplot(2, 2, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, tas_djf_rcm_era5, levels=levs1, latlon=True, cmap=cm.bwr)
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
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_exp2_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


