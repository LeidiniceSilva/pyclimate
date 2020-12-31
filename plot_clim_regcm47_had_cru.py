# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology graphics from Rec_EXP models end OBS basedata"


import os
import netCDF4
import statistics
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from sklearn import metrics
from scipy.stats import norm
from matplotlib.font_manager import FontProperties


def import_rcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_reg_had_{2}_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	rcm_clim = []
	for mon in range(1, 12 + 1):
		rcm = np.nanmean(value[mon::12], axis=0)
		rcm_clim.append(rcm)
	
	return rcm_clim


def import_gcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_amz_neb_Amon_HadGEM2-ES_{2}_r1i1p1_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm_clim = []
	for mon in range(1, 12 + 1):
		gcm = np.nanmean(value[mon::12], axis=0)
		gcm_clim.append(gcm)
	
	return gcm_clim

	
def import_obs(var, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_amz_neb_{2}_obs_mon_{3}_lonlat.nc'.format(path, var, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	obs_clim = []
	for mon in range(1, 12 + 1):
		obs = np.nanmean(value[mon::12], axis=0)
		obs_clim.append(obs)
	
	return obs_clim
	
	
               
# Import regcm exps model end obs database climatology
p_reg_clim = import_rcm('pr', 'hist', '1986-2005')
p_had_clim = import_gcm('pr', 'hist', '1986-2005')
p_cru_clim = import_obs('pre', 'cru_ts4.04', '1986-2005')
p_udel_clim = import_obs('pre', 'udel_v301', '1986-2005')
p_chirps_clim = import_obs('precip', 'chirps-v2.0', '1986-2005')
p_era5_clim = import_obs('mtpr', 'era5', '1986-2005')

t_reg_clim = import_rcm('pr', 'hist', '1986-2005')
t_had_clim = import_gcm('pr', 'hist', '1986-2005')
t_cru_clim = import_obs('pre', 'cru_ts4.04', '1986-2005')
t_udel_clim = import_obs('pre', 'udel_v301', '1986-2005')
t_chirps_clim = import_obs('precip', 'chirps-v2.0', '1986-2005')
t_era5_clim = import_obs('mtpr', 'era5', '1986-2005')

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
plt_clim1 = ax1.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim1
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

ax2 = fig.add_subplot(3, 2, 2)
plt_clim2 = ax2.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'B) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim2
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

legend = ['Reg','Had','CRU','UDEL','CHIRPS','ERA5']
plt.legend(plt_clim2, legend, loc=(1.019, -0.27))

ax3 = fig.add_subplot(3, 2, 3)
plt_clim3 = ax3.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'C) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim3
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

ax4 = fig.add_subplot(3, 2, 4)
plt_clim4 = ax4.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'D) ENEB', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim4
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

ax5 = fig.add_subplot(3, 2, 5)
plt_clim5 = ax5.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'E) MATOPIBA', fontweight='bold')
plt.xlabel('Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim5
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

ax6 = fig.add_subplot(3, 2, 6)
plt_clim6 = ax6.plot(time, p_reg_clim, time, p_had_clim, time, p_cru_clim, time, p_udel_clim, time, p_chirps_clim, time, p_era5_clim)
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xlabel('Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = plt_clim6
plt.setp(l1, linestyle='dashed', color='black')
plt.setp(l2, linestyle='dashed', color='gray')
plt.setp(l3, color='blue')
plt.setp(l4, color='green')
plt.setp(l5, color='red')
plt.setp(l6, color='yellow')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_maps_clim_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






