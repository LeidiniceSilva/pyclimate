# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/09/2021"
__description__ = "This script plot vertical integrated moisture flux convergence"

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
	obs  = np.nanmean(var[:][:,27:36,:,:], axis=0)

	return lat, lon, obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2'
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm = np.nanmean(var[:][:,9:17,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp2'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm = np.nanmean(var[:][:,0:4,:,:], axis=0)

	return lat, lon, gcm
	

def comp_mfc(hus1, ua1, va1):
	
	p1 = ua1 + va1
	p2 = hus1 * p1
	p3 = -1 * p2
	
	return p3
	
		
def comp_vimfc(hus1, ua1, va1, hus2, ua2, va2):
	
	p1 = hus1 * (ua1 + va1) 
	p2 = hus2 * (ua2 + va2) 
	p3 = p1 - p2
	p4 = -0.1 * p3
	
	return p4

		
def comp_vimfc_rea(hus, ua, va):
	
	p1 = hus[0,:,:] * (ua[0,:,:] + va[0,:,:])
	p2 = hus[1,:,:] * (ua[1,:,:] + va[1,:,:])
	p3 = hus[2,:,:] * (ua[2,:,:] + va[2,:,:])
	p4 = hus[3,:,:] * (ua[3,:,:] + va[3,:,:])
	p5 = hus[4,:,:] * (ua[4,:,:] + va[4,:,:])
	p6 = hus[5,:,:] * (ua[5,:,:] + va[5,:,:])
	p7 = hus[6,:,:] * (ua[6,:,:] + va[6,:,:])
	p8 = hus[7,:,:] * (ua[7,:,:] + va[7,:,:])
	p9 = hus[8,:,:] * (ua[8,:,:] + va[8,:,:])
	
	p10 = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9
	p11 = 300 * p10
	
	return p11
		
	
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

	return map, xx, yy
		
	
# Import models and obs database 	   
lat, lon, hus_djf_rea = import_obs('q', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, hus_djf_rcm = import_rcm('hus', 'amz_neb', 'historical', '1986-2005')
lat, lon, hus_djf_gcm = import_gcm('hus', 'amz_neb', 'historical', '1986-2005')

lat, lon, ua_djf_rea = import_obs('u', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, ua_djf_rcm = import_rcm('ua', 'amz_neb', 'historical', '1986-2005')
lat, lon, ua_djf_gcm = import_gcm('ua', 'amz_neb', 'historical', '1986-2005')

lat, lon, va_djf_rea = import_obs('v', 'amz_neb', 'era5', '1986-2005')	   
lat, lon, va_djf_rcm = import_rcm('va', 'amz_neb', 'historical', '1986-2005')
lat, lon, va_djf_gcm = import_gcm('va', 'amz_neb', 'historical', '1986-2005')

print(hus_djf_rea.shape)
print(hus_djf_rcm.shape)
print(hus_djf_gcm.shape)

# Compute Vertical Integrated Moisture Flux Convergence between 1000hPa and 300hPa
vimfc_rea = comp_vimfc_rea(hus_djf_rea, ua_djf_rea, va_djf_rea)
print(vimfc_rea.shape)
print(np.min(vimfc_rea), np.max(vimfc_rea))
exit()

#~ # Compute Moisture Flux Convergence at 850hPa
#~ mfc_rea = comp_mfc(hus_djf_rea_850hPa, ua_djf_rea_850hPa, va_djf_rea_850hPa)
#~ mfc_rcm = comp_mfc(hus_djf_rcm_850hPa, ua_djf_rcm_850hPa, va_djf_rcm_850hPa)
#~ mfc_gcm = comp_mfc(hus_djf_gcm_850hPa, ua_djf_gcm_850hPa, va_djf_gcm_850hPa)

#~ # Compute Vertical Integrated Moisture Flux Convergence between 850hPa and 200hPa
#~ vimfc_rea = comp_vimfc(hus_djf_rea_850hPa, ua_djf_rea_850hPa, va_djf_rea_850hPa, hus_djf_rea_200hPa, ua_djf_rea_200hPa, va_djf_rea_200hPa)
#~ vimfc_rcm = comp_vimfc(hus_djf_rcm_850hPa, ua_djf_rcm_850hPa, va_djf_rcm_850hPa, hus_djf_rcm_200hPa, ua_djf_rcm_200hPa, va_djf_rcm_200hPa)
#~ vimfc_gcm = comp_vimfc(hus_djf_gcm_850hPa, ua_djf_gcm_850hPa, va_djf_gcm_850hPa, hus_djf_gcm_200hPa, ua_djf_gcm_200hPa, va_djf_gcm_200hPa)

# Plot models and obs database 
fig = plt.figure()
levs = [0, 0.001, 0.003, 0.006, 0.009, 0.01, 0.016, 0.019]

ax1 = fig.add_subplot(3, 1, 1)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimfc_rea, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rea_850hPa[::10,::10], va_djf_rea_850hPa[::10,::10])
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black') 
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)

ax2 = fig.add_subplot(3, 1, 2)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimfc_rcm, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_rcm_850hPa[::10,::10], va_djf_rcm_850hPa[::10,::10]) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)

ax3 = fig.add_subplot(3, 1, 3)
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, vimfc_gcm, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
map.quiver(xx[::10,::10], yy[::10,::10], ua_djf_gcm_850hPa[::10,::10], va_djf_gcm_850hPa[::10,::10]) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_vimfc.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






