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
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_djf_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	obs_850hPa = np.nanmean(var[:][:,30,:,:], axis=0)
	obs_200hPa = np.nanmean(var[:][:,16,:,:], axis=0)
	
	return lat, lon, obs_850hPa, obs_200hPa
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	rcm_850hPa = np.nanmean(var[:][:,12,:,:], axis=0)
	rcm_200hPa = np.nanmean(var[:][:,4,:,:], axis=0)

	return lat, lon, rcm_850hPa, rcm_200hPa
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp2'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	gcm_850hPa = np.nanmean(var[:][:,2,:,:], axis=0)
	gcm_200hPa = np.nanmean(var[:][:,8,:,:], axis=0)

	return lat, lon, gcm_850hPa, gcm_200hPa
	

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
lat, lon, hus_djf_rea_850hPa, hus_djf_rea_200hPa = import_obs('q', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, hus_djf_rcm_850hPa, hus_djf_rcm_200hPa = import_rcm('hus', 'amz_neb', 'historical', '1986-2005')
lat, lon, hus_djf_gcm_850hPa, hus_djf_gcm_200hPa = import_gcm('hus', 'amz_neb', 'historical', '1986-2005')

lat, lon, ua_djf_rea_850hPa, ua_djf_rea_200hPa = import_obs('u', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, ua_djf_rcm_850hPa, ua_djf_rcm_200hPa = import_rcm('ua', 'amz_neb', 'historical', '1986-2005')
lat, lon, ua_djf_gcm_850hPa, ua_djf_gcm_200hPa = import_gcm('ua', 'amz_neb', 'historical', '1986-2005')

lat, lon, va_djf_rea_850hPa, va_djf_rea_200hPa = import_obs('v', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, va_djf_rcm_850hPa, va_djf_rcm_200hPa = import_rcm('va', 'amz_neb', 'historical', '1986-2005')
lat, lon, va_djf_gcm_850hPa, va_djf_gcm_200hPa = import_gcm('va', 'amz_neb', 'historical', '1986-2005')

wind_djf_rea = ua_djf_rea_850hPa - ua_djf_rea_200hPa
wind_djf_rcm = ua_djf_rcm_850hPa - ua_djf_rcm_200hPa
wind_djf_gcm = ua_djf_gcm_850hPa - ua_djf_gcm_200hPa

# Plot models and obs database 
fig = plt.figure(figsize=(8,4))
levs1 = [2, 4, 6, 8, 10, 12, 14, 16]

ax1 = fig.add_subplot(1, 3, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, hus_djf_rea_850hPa*1000, levels=levs1, latlon=True, cmap=cm.RdYlGn)
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rea_850hPa[::10,::10], va_djf_rea_850hPa[::10,::10]) 
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

y1,j1 = map(-35,-4)
y2,j2 = map(-35,7)
y3,j3 = map(-19,7)
y4,j4 = map(-19,-4)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.0)
plt.gca().add_patch(poly2)
	
ax2 = fig.add_subplot(1, 3, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, hus_djf_rcm_850hPa*1000, levels=levs1, latlon=True, cmap=cm.RdYlGn)
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rcm_850hPa[::10,::10], va_djf_rcm_850hPa[::10,::10]) 
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

y1,j1 = map(-35,-4)
y2,j2 = map(-35,7)
y3,j3 = map(-19,7)
y4,j4 = map(-19,-4)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.0)
plt.gca().add_patch(poly2)

ax2 = fig.add_subplot(1, 3, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, hus_djf_gcm_850hPa*1000, levels=levs1, latlon=True, cmap=cm.RdYlGn, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_gcm_850hPa[::10,::10], va_djf_gcm_850hPa[::10,::10]) 
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

y1,j1 = map(-35,-4)
y2,j2 = map(-35,7)
y3,j3 = map(-19,7)
y4,j4 = map(-19,-4)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.0)
plt.gca().add_patch(poly2)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_hus.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

# Plot models and obs database 
fig = plt.figure(figsize=(8,4))
levs1 = [-20, -15, -10, -5, 5, 10, 15, 20]

ax1 = fig.add_subplot(1, 3, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, wind_djf_rea, levels=levs1, latlon=True, cmap=cm.PuOr)
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rea_200hPa[::10,::10], va_djf_rea_200hPa[::10,::10]) 
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax2 = fig.add_subplot(1, 3, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, wind_djf_rcm, levels=levs1, latlon=True, cmap=cm.PuOr)
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rcm_200hPa[::10,::10], va_djf_rcm_200hPa[::10,::10]) 
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

ax2 = fig.add_subplot(1, 3, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, wind_djf_gcm, levels=levs1, latlon=True, cmap=cm.PuOr, extend='both')
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_gcm_200hPa[::10,::10], va_djf_gcm_200hPa[::10,::10]) 
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_wind.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


