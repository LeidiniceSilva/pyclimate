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
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2/'
	arq  = '{0}/{1}_{2}_era5_obs_djf_{3}_lonlat.nc'.format(path, var, area, dt)	
			
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lev  = data.variables['level'][:]
	value = data.variables[var][:,:,:,0] 
	value = np.nanmean(value, axis=0)

	return lat, lev, value
	
		
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lev  = data.variables['plev'][:]
	value = data.variables[var][:,:,:,0] 
	value = np.nanmean(value, axis=0)

	return lat, lev, value
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_djf_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lev  = data.variables['plev'][:]
	value = data.variables[var][:,:,:,0] 
	value = np.nanmean(value, axis=0)

	return lat, lev, value
	
		
# Import models and obs database    
lat, w_lev_rea, w_var_rea = import_obs('w', 'amz_neb', '1986-2005')
lat, w_lev_rcm, w_var_rcm = import_rcm('omega', 'amz_neb', 'historical', '1986-2005')
lat, w_lev_gcm, w_var_gcm = import_gcm('wap', 'amz_neb', 'historical', '1986-2005')

print(lat.shape, w_lev_rea.shape, w_var_rea.shape)
print(lat.shape, w_lev_rcm.shape, w_var_rcm.shape)
print(lat.shape, w_lev_gcm.shape, w_var_gcm.shape)

# Plot models and obs database 
fig = plt.figure(figsize=(8,6))
xtime = np.arange(-20, 13, 3)

ax = fig.add_subplot(3, 1, 1)
ytime = np.arange(100, 1100, 100)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lat, w_lev_rea, w_var_rea*1000, np.arange(-60, 70, 10), cmap=plt.cm.PiYG, extend='both')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lat, w_lev_rea, w_var_rea*1000, np.arange(-60, 70, 10), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%d', fontsize=8, colors='black')
plt.ylabel('Pressure (hPa)', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '17°S', '14°S', '11°S', '8°S', '5°S', '2°S', '1°N', '4°N', '7°S', '10°N'), fontsize=8)
plt.yticks(ytime, ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.gca().invert_yaxis()
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(-10, linewidth=1., linestyle='--', color='black')

ax = fig.add_subplot(3, 1, 2)
ytime = np.arange(100, 1100, 100)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lat, w_lev_rcm, w_var_rcm*100000, np.arange(-60, 70, 10), cmap=plt.cm.PiYG, extend='both')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lat, w_lev_rcm, w_var_rcm*100000, np.arange(-60, 70, 10), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%d', fontsize=8, colors='black')
plt.ylabel('Pressure (hPa)', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '17°S', '14°S', '11°S', '8°S', '5°S', '2°S', '1°N', '4°N', '7°S', '10°N'), fontsize=8)
plt.yticks(ytime, ('100', '2plot_maps_clim_regcm47_obs_exp2_vert.py00', '300', '400', '500', '600', '700', '800', '900', '1000'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.gca().invert_yaxis()
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(-10, linewidth=1., linestyle='--', color='black')

ax = fig.add_subplot(3, 1, 3)
ytime = np.arange(10000, 110000, 10000)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt1 = ax.contourf(lat, w_lev_gcm, w_var_gcm*1000, np.arange(-60, 70, 10), cmap=plt.cm.PiYG, extend='both')
plt.colorbar(plt1, aspect=10, pad=0.025)
plt2 = ax.contour(lat, w_lev_gcm, w_var_gcm*1000, np.arange(-60, 70, 10), linewidths=0.5, colors='black')
ax.clabel(plt2, fmt='%d', fontsize=8, colors='black')
plt.ylabel('Pressure (hPa)', fontsize=8, fontweight='bold')
plt.xlabel('Latitude', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '17°S', '14°S', '11°S', '8°S', '5°S', '2°S', '1°N', '4°N', '7°S', '10°N'), fontsize=8)
plt.yticks(ytime, ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'), fontsize=8)
plt.gca().invert_yaxis()
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(-10, linewidth=1., linestyle='--', color='black')
																				
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_omega.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


