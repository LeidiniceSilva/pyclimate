# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/09/2021"
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
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_djf_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs_850hPa = np.nanmean(var[:][:,16:30,:,:], axis=0)
	obs_200hPa = np.nanmean(var[:][:,16,:,:], axis=0)

	return lat, lon, obs_850hPa, obs_200hPa
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:,:]
	obs_850hPa = np.nanmean(var[:][:,12,:,:], axis=0)
	obs_200hPa = np.nanmean(var[:][:,4,:,:], axis=0)

	return lat, lon, obs_850hPa, obs_200hPa


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp2'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:,:]
	gcm_850hPa = np.nanmean(var[:][:,2,:,:], axis=0)
	gcm_200hPa = np.nanmean(var[:][:,8,:,:], axis=0)

	return lat, lon, gcm_850hPa, gcm_200hPa
	

def comp_vimf(hus1, ua1, va1):
	
	p1 = (ua1 * hus1) / 9.8
	p2 = (va1 * hus1) / 9.8
	p3 = (p1 ** 2 + p2 ** 2) ** 0.5
	
	return p3
	
	
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
lat, lon, hus_djf_era5_850hPa, hus_djf_era5_250hPa = import_obs('q', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, hus_djf_rcm_850hPa, hus_djf_rcm_250hPa = import_rcm('hus', 'amz_neb', 'historical', '1986-2005')
lat, lon, hus_djf_gcm_850hPa, hus_djf_gcm_250hPa = import_gcm('hus', 'amz_neb', 'historical', '1986-2005')

lat, lon, ua_djf_era5_850hPa, ua_djf_era5_250hPa = import_obs('u', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, ua_djf_rcm_850hPa, ua_djf_rcm_250hPa = import_rcm('ua', 'amz_neb', 'historical', '1986-2005')
lat, lon, ua_djf_gcm_850hPa, ua_djf_gcm_250hPa = import_gcm('ua', 'amz_neb', 'historical', '1986-2005')

lat, lon, va_djf_era5_850hPa, va_djf_era5_250hPa = import_obs('v', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, va_djf_rcm_850hPa, va_djf_rcm_250hPa = import_rcm('va', 'amz_neb', 'historical', '1986-2005')
lat, lon, va_djf_gcm_850hPa, va_djf_gcm_250hPa = import_gcm('va', 'amz_neb', 'historical', '1986-2005')

# Compute VIMF
vimf_era5 = comp_vimf(hus_djf_era5_850hPa, ua_djf_era5_850hPa, va_djf_era5_850hPa)

print(vimf_era5.shape)
exit()
print(np.min(vimf_era5))
print(np.max(vimf_era5))
exit()

# Plot models and obs database 
fig = plt.figure(figsize=(8,3))
levs = [0, 3, 6, 9, 12, 15, 18, 21]

ax1 = fig.add_subplot(1, 3, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimf_era5, levels=levs, latlon=True, cmap=cm.Greens)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')

ax2 = fig.add_subplot(1, 3, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimf_era5, levels=levs, latlon=True, cmap=cm.Greens)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

ax3 = fig.add_subplot(1, 3, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimf_era5, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_vimf.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






