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


def import_obs(var, area, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2'
	arq  = '{0}/{1}_{2}_era5_obs_djf_{3}_lonlat.nc'.format(path, var, area, dt)	
			
	data = netCDF4.Dataset(arq)
	lon  = data.variables['lon'][:]
	lev  = data.variables['level'][:]
	value = data.variables[var][:,:,0,:] 
	value = np.nanmean(value, axis=0)

	return lon, lev, value
	
		
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lon  = data.variables['lon'][:]
	lev  = data.variables['plev'][:]
	value = data.variables[var][:,:,0,:] 
	value = np.nanmean(value, axis=0)

	return lon, lev, value
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Downloads'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lev  = data.variables['plev'][:]
	value = data.variables[var][:,:,:,0] 
	value = np.nanmean(value, axis=0)

	return lat, lev, value
	
		
# Import models and obs database    
lon, hur_lev_rea, hur_var_rea = import_obs('r', 'amz_neb', '1986-2005')
lon, hur_lev_rcm, hur_var_rcm = import_rcm('rh', 'amz_neb', 'historical', '1986-2005')
lat, hur_lev_gcm, hur_var_gcm = import_gcm('wap', 'amz_neb', 'historical', '1986-2005')

print(lon.shape, hur_lev_rea.shape, hur_var_rea.shape)
print(lon.shape, hur_lev_rcm.shape, hur_var_rcm.shape)
print(lat.shape, hur_lev_gcm.shape, hur_var_gcm.shape)

# Plot models and obs database 
fig = plt.figure(figsize=(8,6))
time = np.arange(-85, -15, 5)
xtime = np.arange(-20, 15, 5)
ytime = np.arange(10000, 110000, 10000)

ax = fig.add_subplot(3, 1, 1)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lon, hur_lev_rea, hur_var_rea, np.arange(10, 90, 10), cmap=plt.cm.Blues, extend='max')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lon, hur_lev_rea, hur_var_rea, np.arange(10, 90, 10), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%d', fontsize=8, colors='black')
plt.ylabel('Pressure (hPa)', fontsize=8, fontweight='bold')
plt.ylim(100, 1000)
plt.xticks(time, ('85°W', '80°W', '75°W', '70°W', '65°W', '60°W', '55°W', '50°W', '45°W', '40°W', '35°W', '30°W', '25°W', '20°W', '15°W'), fontsize=8)
plt.yticks(np.arange(100, 1100, 100), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(3, 1, 2)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lon, hur_lev_rcm, hur_var_rcm, np.arange(10, 90, 10), cmap=plt.cm.Blues, extend='max')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lon, hur_lev_rcm, hur_var_rcm, np.arange(10, 90, 10), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%d', fontsize=8, colors='black')
plt.ylabel('Pressure (hPa)', fontsize=8, fontweight='bold')
plt.ylim(100, 1000)
plt.xticks(time, ('85°W', '80°W', '75°W', '70°W', '65°W', '60°W', '55°W', '50°W', '45°W', '40°W', '35°W', '30°W', '25°W', '20°W', '15°W'), fontsize=8)
plt.yticks(np.arange(100, 1100, 100), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(3, 1, 3)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lat, hur_lev_gcm, hur_var_gcm, np.arange(-0.10, 0.10, 0.02), cmap=plt.cm.seismic, extend='both')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lat, hur_lev_gcm, hur_var_gcm, np.arange(-0.10, 0.10, 0.02), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%.02f', fontsize=8, colors='black')
plt.ylabel('Vertical velocity (hPa)', fontsize=8, fontweight='bold')
plt.xlabel('Latitude', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '15°S', '10°S', '5°S', '0°', '5°N', '10°N'), fontsize=8)
plt.yticks(ytime, ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'), fontsize=8)
plt.gca().invert_yaxis()
																				
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_hur.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


