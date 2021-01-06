# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual cycle graphics from Reg and Had models end obs database"

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


def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	rcm = []
	for mon in range(1, 12 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))
	
	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
	
	return gcm

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	obs = []
	for mon in range(1, 12 + 1):
		obs.append(np.nanmean(value[mon::12], axis=0))
	
	return obs
	
               
# Import regcm exps model end obs database climatology
p_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
p_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
p_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
p_udel_samz = import_obs('pre', 'samz', 'udel_v301', '1986-2005')
p_chirps_samz = import_obs('precip', 'samz', 'chirps-v2.0', '1986-2005')
p_era5_samz = import_obs('mtpr', 'samz', 'era5', '1986-2005')

p_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
p_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
p_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
p_udel_eneb = import_obs('pre', 'eneb', 'udel_v301', '1986-2005')
p_chirps_eneb = import_obs('precip', 'eneb', 'chirps-v2.0', '1986-2005')
p_era5_eneb = import_obs('mtpr', 'eneb', 'era5', '1986-2005')

p_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
p_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
p_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
p_udel_matopiba = import_obs('pre', 'matopiba', 'udel_v301', '1986-2005')
p_chirps_matopiba = import_obs('precip', 'matopiba', 'chirps-v2.0', '1986-2005')
p_era5_matopiba = import_obs('mtpr', 'matopiba', 'era5', '1986-2005')

t_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
t_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
t_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
t_udel_samz = import_obs('temp', 'samz', 'udel_v301', '1986-2005')
t_era5_samz = import_obs('t2m', 'samz', 'era5', '1986-2005')

t_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
t_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
t_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
t_udel_eneb = import_obs('temp', 'eneb', 'udel_v301', '1986-2005')
t_era5_eneb = import_obs('t2m', 'eneb', 'era5', '1986-2005')

t_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
t_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
t_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
t_udel_matopiba = import_obs('temp', 'matopiba', 'udel_v301', '1986-2005')
t_era5_matopiba = import_obs('t2m', 'matopiba', 'era5', '1986-2005')

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
annual_cycle1 = ax1.plot(time, p_reg_samz, time, p_had_samz, time, p_cru_samz, time, p_udel_samz, time, p_chirps_samz, time, p_era5_samz)
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^')
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='red', marker='*')
plt.setp(l6, linewidth=1, color='yellow', marker='x')

ax2 = fig.add_subplot(3, 2, 2)
annual_cycle2 = ax2.plot(time, t_had_samz, time, t_had_samz, time, t_cru_samz, time, t_udel_samz, time, t_era5_samz)
plt.title(u'D) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 32, 2))
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle2
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^')
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='yellow', marker='*')

ax3 = fig.add_subplot(3, 2, 3)
annual_cycle3 = ax3.plot(time, p_reg_eneb, time, p_had_eneb, time, p_cru_eneb, time, p_udel_eneb, time, p_chirps_eneb, time, p_era5_eneb)
plt.title(u'B) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle3
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^')
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='red', marker='*')
plt.setp(l6, linewidth=1, color='yellow', marker='x')

ax4 = fig.add_subplot(3, 2, 4)
annual_cycle4 = ax4.plot(time, t_had_eneb, time, t_had_eneb, time, t_cru_eneb, time, t_udel_eneb, time, t_era5_eneb)
plt.title(u'E) ENEB', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 32, 2))
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle4
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^', )
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='yellow', marker='*')

ax5 = fig.add_subplot(3, 2, 5)
annual_cycle5 = ax5.plot(time, p_reg_matopiba, time, p_had_matopiba, time, p_cru_matopiba, time, p_udel_matopiba, time, p_chirps_matopiba, time, p_era5_matopiba)
plt.title(u'C) MATOPIBA', fontweight='bold')
plt.xlabel('Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 12, 2))
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle5
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^')
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='red', marker='*')
plt.setp(l6, linewidth=1, color='yellow', marker='x')

legend = ['Reg','Had','CRU','UDEL','CHIRPS','ERA5']
plt.legend(annual_cycle5, legend, loc='lower left', bbox_to_anchor=(-0.3, -0.8), shadow=True, ncol=6)

ax6 = fig.add_subplot(3, 2, 6)
annual_cycle6 = ax6.plot(time, t_had_matopiba, time, t_had_matopiba, time, t_cru_matopiba, time, t_udel_matopiba, time, t_era5_matopiba)
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xlabel('Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 32, 2))
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle6
plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
plt.setp(l2, linewidth=1, color='gray', marker='.', linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='^')
plt.setp(l4, linewidth=1, color='green', marker='+')
plt.setp(l5, linewidth=1, color='yellow', marker='*')

fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






