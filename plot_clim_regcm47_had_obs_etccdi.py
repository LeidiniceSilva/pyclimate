# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual climatology from from extremes indices"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	

	dict_var = {u'precip': u'precip', 
	u'tmax': u'tmax',
	u'tmin': u'tmin'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	obs = []
	for mon in range(1, 12 + 1):
		obs.append(np.nanmean(value[mon::12], axis=0))
			
	return obs
	
	
def import_rcm(var, area, model, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_mon_{5}_lonlat_seamask.nc'.format(path, var, area, model, exp, dt)	

	dict_var = {u'pr': u'pr', 
	u'tasmax': u'tasmax',
	u'tasmin': u'tasmin'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	rcm = []
	for mon in range(1, 12 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))
			
	return rcm


def import_gcm(var, area, model, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_mon_{5}_lonlat_seamask.nc'.format(path, var, area, model, exp, dt)	

	dict_var = {u'pr': u'pr', 
	u'tasmax': u'tasmax',
	u'tasmin': u'tasmin'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
			
	return gcm
		
	            
# Import extreme indices 
mon_pre_cpc_samz = import_obs('precip', 'samz', 'cpc_obs', '1986-2005')
mon_pre_cpc_eneb = import_obs('precip', 'eneb', 'cpc_obs', '1986-2005')
mon_pre_cpc_matopiba = import_obs('precip', 'matopiba', 'cpc_obs', '1986-2005')
mon_pre_reg_samz = import_rcm('pr', 'samz', 'RegCM47_had', 'historical', '1986-2005')
mon_pre_reg_eneb = import_rcm('pr', 'eneb', 'RegCM47_had', 'historical', '1986-2005')
mon_pre_reg_matopiba = import_rcm('pr', 'matopiba', 'RegCM47_had', 'historical', '1986-2005')
mon_pre_had_samz = import_gcm('pr', 'samz', 'HadGEM2-ES', 'historical', '1986-2005')
mon_pre_had_eneb = import_gcm('pr', 'eneb', 'HadGEM2-ES', 'historical', '1986-2005')
mon_pre_had_matopiba = import_gcm('pr', 'matopiba', 'HadGEM2-ES', 'historical', '1986-2005')

mon_tmax_cpc_samz = import_obs('tmax', 'samz', 'cpc_obs', '1986-2005')
mon_tmax_cpc_eneb = import_obs('tmax', 'eneb', 'cpc_obs', '1986-2005')
mon_tmax_cpc_matopiba = import_obs('tmax', 'matopiba', 'cpc_obs', '1986-2005')
mon_tmax_reg_samz = import_rcm('tasmax', 'samz', 'RegCM47_had', 'historical', '1986-2005')
mon_tmax_reg_eneb = import_rcm('tasmax', 'eneb', 'RegCM47_had', 'historical', '1986-2005')
mon_tmax_reg_matopiba = import_rcm('tasmax', 'matopiba', 'RegCM47_had', 'historical', '1986-2005')
mon_tmax_had_samz = import_gcm('tasmax', 'samz', 'HadGEM2-ES', 'historical', '1986-2005')
mon_tmax_had_eneb = import_gcm('tasmax', 'eneb', 'HadGEM2-ES', 'historical', '1986-2005')
mon_tmax_had_matopiba = import_gcm('tasmax', 'matopiba', 'HadGEM2-ES', 'historical', '1986-2005')

mon_tmin_cpc_samz = import_obs('tmin', 'samz', 'cpc_obs', '1986-2005')
mon_tmin_cpc_eneb = import_obs('tmin', 'eneb', 'cpc_obs', '1986-2005')
mon_tmin_cpc_matopiba = import_obs('tmin', 'matopiba', 'cpc_obs', '1986-2005')
mon_tmin_reg_samz = import_rcm('tasmin', 'samz', 'RegCM47_had', 'historical', '1986-2005')
mon_tmin_reg_eneb = import_rcm('tasmin', 'eneb', 'RegCM47_had', 'historical', '1986-2005')
mon_tmin_reg_matopiba = import_rcm('tasmin', 'matopiba', 'RegCM47_had', 'historical', '1986-2005')
mon_tmin_had_samz = import_gcm('tasmin', 'samz', 'HadGEM2-ES', 'historical', '1986-2005')
mon_tmin_had_eneb = import_gcm('tasmin', 'eneb', 'HadGEM2-ES', 'historical', '1986-2005')
mon_tmin_had_matopiba = import_gcm('tasmin', 'matopiba', 'HadGEM2-ES', 'historical', '1986-2005')

# Plot extreme indices 
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 1, 1)
plt_clim1 = plt.bar(time, mon_pre_cpc_samz, color='gray', label='CPC', width=0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time + .25, mon_pre_reg_samz,  color='red', label='Reg', width=0.25, edgecolor='black', linewidth=1)
plt_clim3 = plt.bar(time + .50, mon_pre_had_samz,  color='blue', label='Had', width=0.25, edgecolor='black', linewidth=1)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.ylim(0, 500)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 600, 100), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)

ax1_1 = ax1.twinx()
line1_1, = ax1_1.plot(time + .25, mon_tmax_cpc_samz, color='gray', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line1_1, = ax1_1.plot(time + .25, np.nanmean(mon_tmax_reg_samz, axis=1), color='red', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line1_1, = ax1_1.plot(time + .25, mon_tmax_had_samz, color='blue', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line1_1, = ax1_1.plot(time + .25, mon_tmin_cpc_samz, color='gray', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line1_1, = ax1_1.plot(time + .25, np.nanmean(mon_tmin_reg_samz, axis=1), color='red', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line1_1, = ax1_1.plot(time + .25, mon_tmin_had_samz, color='blue', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
plt.ylim(15, 40)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(15, 45, 5), fontsize=8)
plt.setp(ax1_1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 1, 2)
plt_clim1 = plt.bar(time, mon_pre_cpc_eneb, color='gray', label='CPC', width=0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time + .25, mon_pre_reg_eneb,  color='red', label='Reg', width=0.25, edgecolor='black', linewidth=1)
plt_clim3 = plt.bar(time + .50, mon_pre_had_eneb,  color='blue', label='Had', width=0.25, edgecolor='black', linewidth=1)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation (mm d⁻¹)', fontsize=8)
plt.ylim(0, 500)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 600, 100), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.legend(fontsize=8, shadow=True, ncol=3, handlelength=0.75, handleheight=0.75)

ax2_1 = ax2.twinx()
line2_1, = ax2_1.plot(time + .25, mon_tmax_cpc_eneb, color='gray', label='Tmax', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line2_1, = ax2_1.plot(time + .25, np.nanmean(mon_tmax_reg_eneb, axis=1), color='red', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line2_1, = ax2_1.plot(time + .25, mon_tmax_had_eneb, color='blue', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line2_1, = ax2_1.plot(time + .25, mon_tmin_cpc_eneb, color='gray', label='Tmin', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line2_1, = ax2_1.plot(time + .25, np.nanmean(mon_tmin_reg_eneb, axis=1), color='red', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line2_1, = ax2_1.plot(time + .25, mon_tmin_had_eneb, color='blue', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
plt.ylabel(u'Temperature (°C)', fontsize=8)
plt.ylim(15, 40)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(15, 45, 5), fontsize=8)
plt.setp(ax1_1.get_xticklabels(), visible=False)
plt.legend(fontsize=8, shadow=True, ncol=3, loc=2)

ax3 = fig.add_subplot(3, 1, 3)
plt_clim1 = plt.bar(time, mon_pre_cpc_matopiba, color='gray', label='CPC', width=0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time + .25, mon_pre_reg_matopiba,  color='red', label='Reg', width=0.25, edgecolor='black', linewidth=1)
plt_clim3 = plt.bar(time + .50, mon_pre_had_matopiba,  color='blue', label='Had', width=0.25, edgecolor='black', linewidth=1)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Months', fontsize=8)
plt.ylim(0, 500)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 600, 100), fontsize=8)

ax3_1 = ax3.twinx()
line3_1, = ax3_1.plot(time + .25, mon_tmax_cpc_matopiba, color='gray', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line3_1, = ax3_1.plot(time + .25, np.nanmean(mon_tmax_reg_matopiba, axis=1), color='red', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line3_1, = ax3_1.plot(time + .25, mon_tmax_had_matopiba, color='blue', linestyle='-', marker='^', markersize=4, markerfacecolor='white', linewidth=1)
line3_1, = ax3_1.plot(time + .25, mon_tmin_cpc_matopiba, color='gray', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line3_1, = ax3_1.plot(time + .25, np.nanmean(mon_tmin_reg_matopiba, axis=1), color='red', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
line3_1, = ax3_1.plot(time + .25, mon_tmin_had_matopiba, color='blue', linestyle='-', marker='v', markersize=4, markerfacecolor='white', linewidth=1)
plt.ylim(15, 40)
plt.xticks(time + .25, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(15, 45, 5), fontsize=8)
plt.setp(ax1_1.get_xticklabels(), visible=False)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_had_obs_1986-2005_etccdi.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






