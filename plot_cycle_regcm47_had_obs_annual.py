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

def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
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

	              
# Import regcm exps model end obs database climatology
pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
pre_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
pre_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
pre_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
tas_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
tas_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

tas_cru_matopiba  = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
tas_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
annual_cycle1 = ax1.plot(time, pre_cru_samz, time, pre_reg_samz, time, pre_had_samz)
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax1.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')
legend = ['CRU','RegCM4','HadGEM2-ES']
plt.legend(annual_cycle1, legend, fontsize=6, loc=9, shadow=True, ncol=1)

ax2 = fig.add_subplot(3, 2, 2)
annual_cycle1 = ax2.plot(time, tas_cru_samz, time, np.nanmean(tas_reg_samz, axis=1), time, tas_had_samz)
plt.title(u'D) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(22, 34, 2))
plt.setp(ax2.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')

ax3 = fig.add_subplot(3, 2, 3)
annual_cycle1 = ax3.plot(time, pre_cru_eneb, time, pre_reg_eneb, time, pre_had_eneb)
plt.title(u'B) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')

ax4 = fig.add_subplot(3, 2, 4)
annual_cycle1 = ax4.plot(time, tas_cru_eneb, time, np.nanmean(tas_reg_eneb, axis=1), time, tas_had_eneb)
plt.title(u'E) ENEB', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(22, 34, 2))
plt.setp(ax4.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')

ax5 = fig.add_subplot(3, 2, 5)
annual_cycle1 = ax5.plot(time, pre_cru_matopiba, time, pre_reg_matopiba, time, pre_had_matopiba)
plt.title(u'C) MATOPIBA', fontweight='bold')
plt.xlabel(u'Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')

ax6 = fig.add_subplot(3, 2, 6)
annual_cycle1 = ax6.plot(time, tas_cru_matopiba, time, np.nanmean(tas_reg_matopiba, axis=1), time, tas_had_matopiba)
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xlabel(u'Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(22, 34, 2))
plt.grid(True, which='major', linestyle='--')
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1, color='black', marker='+', linestyle='dashed')
plt.setp(l2, linewidth=1, color='black', marker='o', markersize=4, linestyle='dashed')
plt.setp(l3, linewidth=1, color='black', marker='s', markersize=4, linestyle='dashed')

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






