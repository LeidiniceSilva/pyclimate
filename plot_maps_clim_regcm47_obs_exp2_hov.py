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
	arq  = '{0}/{1}_{2}_obs_mon_{3}_lonlat.nc'.format(path, var, area, dt)	
			
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = data.variables[var][:,:,:] 
	
	obs = []
	for mon in range(1, 12 + 1):
		obs.append(np.nanmean(value[mon::12], axis=0))

	obs = np.nanmean(obs, axis=2)

	return lat, obs
	
		
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = data.variables[var][:,:,:] 

	rcm = []
	for mon in range(1, 12 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))

	rcm = np.nanmean(rcm, axis=2)
	
	return lat, rcm
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = data.variables[var][:,:,:] 

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))

	gcm = np.nanmean(gcm, axis=2)
	
	return lat, gcm
	
		
# Import models and obs database    
lat, pr_var_obs = import_obs('precip', 'amz_neb_gpcp_v2.2', '1986-2005')
lat, pr_var_rea = import_obs('mtpr', 'amz_neb_era5', '1986-2005')
lat, pr_var_rcm = import_rcm('pr', 'amz_neb', 'historical', '1986-2005')
lat, pr_var_gcm = import_gcm('pr', 'amz_neb', 'historical', '1986-2005')

print(lat.shape, pr_var_obs.shape)
print(lat.shape, pr_var_rea.shape)
print(lat.shape, pr_var_rcm.shape)
print(lat.shape, pr_var_gcm.shape)

# Plot models and obs database 
fig = plt.figure()
ytime = np.arange(0, 11+1)
xtime = np.arange(-20, 15, 5)

ax = fig.add_subplot(2, 2, 1)
plt1 = ax.contourf(lat, ytime, pr_var_obs, np.arange(0, 16, 2), cmap=plt.cm.Blues, extend='max')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Months', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '15°S', '10°S', '5°S', 'EQ', '5°S', '10°N', '15°N', '20°N'), fontsize=7)
plt.yticks(ytime, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linewidth=0.5, linestyle='--', color='white')
plt.colorbar(plt1, pad=0.040)

ax = fig.add_subplot(2, 2, 2)
plt1 = ax.contourf(lat, ytime, pr_var_rea, np.arange(0, 16, 2), cmap=plt.cm.Blues, extend='max')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
plt.xticks(xtime, ('20°S', '15°S', '10°S', '5°S', 'EQ', '5°S', '10°N', '15°N', '20°N'), fontsize=7)
plt.yticks(ytime, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.grid(True, which='major', linewidth=0.5, linestyle='--', color='white')
plt.colorbar(plt1, pad=0.040)

ax = fig.add_subplot(2, 2, 3)
plt1 = ax.contourf(lat, ytime, pr_var_rcm, np.arange(0, 16, 2), cmap=plt.cm.Blues, extend='max')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Latitude', fontsize=8, fontweight='bold')
plt.ylabel('Months', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '15°S', '10°S', '5°S', 'EQ', '5°S', '10°N', '15°N', '20°N'), fontsize=8)
plt.yticks(ytime, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.grid(True, which='major', linewidth=0.5, linestyle='--', color='white')
plt.colorbar(plt1, pad=0.040)

ax = fig.add_subplot(2, 2, 4)
plt1 = ax.contourf(lat, ytime, pr_var_gcm, np.arange(0, 16, 2), cmap=plt.cm.Blues, extend='max')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Latitude', fontsize=8, fontweight='bold')
plt.xticks(xtime, ('20°S', '15°S', '10°S', '5°S', 'EQ', '5°S', '10°N', '15°N', '20°N'), fontsize=7)
plt.yticks(ytime, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.grid(True, which='major', linewidth=0.5, linestyle='--', color='white')
plt.setp(ax.get_yticklabels(), visible=False)
plt.colorbar(plt1, pad=0.040)
																				
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_reg_exp2_pre_lat.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



