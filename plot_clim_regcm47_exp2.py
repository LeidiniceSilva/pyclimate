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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	obs = []
	conf_int = []		
	for mon in range(1, 12 + 1):
		obs.append(np.nanmean(value[mon::12], axis=0))

	return obs
	
	
def import_rcm_exp1(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/{0}'.format(exp)	
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


def import_rcm_exp2(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2/{0}'.format(exp)	
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
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'	
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


def compute_bias(sim, obs):
	
	bias = []
	for mon in range(0, 12):
		diff = sim[mon] - obs[mon]
		bias.append(diff)
	
	return bias
	
	          
# Import regcm exps model end obs database climatology
# Precipitation
mon_pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
mon_pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
mon_pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_pre_reg_exp1_samz = import_rcm_exp1('pr', 'samz', 'hist', '1986-2005')
mon_pre_reg_exp1_eneb = import_rcm_exp1('pr', 'eneb', 'hist', '1986-2005')
mon_pre_reg_exp1_matopiba = import_rcm_exp1('pr', 'matopiba', 'hist', '1986-2005')

mon_pre_reg_exp2_samz = import_rcm_exp2('pr', 'samz', 'hist', '1986-2005')
mon_pre_reg_exp2_eneb = import_rcm_exp2('pr', 'eneb', 'hist', '1986-2005')
mon_pre_reg_exp2_matopiba = import_rcm_exp2('pr', 'matopiba', 'hist', '1986-2005')

mon_pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
mon_tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
mon_tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
mon_tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_tas_reg_exp1_samz = import_rcm_exp1('tas', 'samz', 'hist', '1986-2005')
mon_tas_reg_exp1_eneb = import_rcm_exp1('tas', 'eneb', 'hist', '1986-2005')
mon_tas_reg_exp1_matopiba = import_rcm_exp1('tas', 'matopiba', 'hist', '1986-2005')

mon_tas_reg_exp2_samz = import_rcm_exp2('tas', 'samz', 'hist', '1986-2005')
mon_tas_reg_exp2_eneb = import_rcm_exp2('tas', 'eneb', 'hist', '1986-2005')
mon_tas_reg_exp2_matopiba = import_rcm_exp2('tas', 'matopiba', 'hist', '1986-2005')

mon_tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Compute bias
# Precipitation
bias_pre_reg_exp1_cru_samz = compute_bias(mon_pre_reg_exp1_samz, mon_pre_cru_samz)
bias_pre_reg_exp1_cru_eneb = compute_bias(mon_pre_reg_exp1_eneb, mon_pre_cru_eneb)
bias_pre_reg_exp1_cru_matopiba = compute_bias(mon_pre_reg_exp1_matopiba, mon_pre_cru_matopiba)

bias_pre_reg_exp2_cru_samz = compute_bias(mon_pre_reg_exp2_samz, mon_pre_cru_samz)
bias_pre_reg_exp2_cru_eneb = compute_bias(mon_pre_reg_exp2_eneb, mon_pre_cru_eneb)
bias_pre_reg_exp2_cru_matopiba = compute_bias(mon_pre_reg_exp2_matopiba, mon_pre_cru_matopiba)

bias_pre_had_cru_samz = compute_bias(mon_pre_had_samz, mon_pre_cru_samz)
bias_pre_had_cru_eneb = compute_bias(mon_pre_had_eneb, mon_pre_cru_eneb)
bias_pre_had_cru_matopiba = compute_bias(mon_pre_had_matopiba, mon_pre_cru_matopiba)

# Temperature
bias_tas_reg_exp1_cru_samz = compute_bias(np.nanmean(mon_tas_reg_exp1_samz, axis=1), mon_tas_cru_samz)
bias_tas_reg_exp1_cru_eneb = compute_bias(np.nanmean(mon_tas_reg_exp1_eneb, axis=1), mon_tas_cru_eneb)
bias_tas_reg_exp1_cru_matopiba = compute_bias(np.nanmean(mon_tas_reg_exp1_matopiba, axis=1), mon_tas_cru_matopiba)

bias_tas_reg_exp2_cru_samz = compute_bias(np.nanmean(mon_tas_reg_exp2_samz, axis=1), mon_tas_cru_samz)
bias_tas_reg_exp2_cru_eneb = compute_bias(np.nanmean(mon_tas_reg_exp2_eneb, axis=1), mon_tas_cru_eneb)
bias_tas_reg_exp2_cru_matopiba = compute_bias(np.nanmean(mon_tas_reg_exp2_matopiba, axis=1), mon_tas_cru_matopiba)

bias_tas_had_cru_samz = compute_bias(mon_tas_had_samz, mon_tas_cru_samz)
bias_tas_had_cru_eneb = compute_bias(mon_tas_had_eneb, mon_tas_cru_eneb)
bias_tas_had_cru_matopiba = compute_bias(mon_tas_had_matopiba, mon_tas_cru_matopiba)

print('a')
print(np.nanmean(bias_pre_reg_exp1_cru_samz))
print(np.nanmean(bias_pre_reg_exp2_cru_samz))
print(np.nanmean(bias_pre_had_cru_samz))
print('b')
print(np.nanmean(bias_pre_reg_exp1_cru_eneb))
print(np.nanmean(bias_pre_reg_exp2_cru_eneb))
print(np.nanmean(bias_pre_had_cru_eneb))
print('c')
print(np.nanmean(bias_pre_reg_exp1_cru_matopiba))
print(np.nanmean(bias_pre_reg_exp2_cru_matopiba))
print(np.nanmean(bias_pre_had_cru_matopiba))
print('d')
print(np.nanmean(bias_tas_reg_exp1_cru_samz))
print(np.nanmean(bias_tas_reg_exp2_cru_samz))
print(np.nanmean(bias_tas_had_cru_samz))
print('e')
print(np.nanmean(bias_tas_reg_exp1_cru_eneb))
print(np.nanmean(bias_tas_reg_exp2_cru_eneb))
print(np.nanmean(bias_tas_had_cru_eneb))
print('f')
print(np.nanmean(bias_tas_reg_exp1_cru_matopiba))
print(np.nanmean(bias_tas_reg_exp2_cru_matopiba))
print(np.nanmean(bias_tas_had_cru_matopiba))

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax10 = fig.add_subplot(3, 2, 1)
annual_cycle1 = ax10.plot(time, bias_pre_reg_exp1_cru_samz, time, bias_pre_reg_exp2_cru_samz, time, bias_pre_had_cru_samz)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
plt.setp(ax10.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle1
plt.setp(l1, linewidth=1.5, color='blue', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.legend(annual_cycle1, ['RegCM4.7_EXP1', 'RegCM4.7_EXP2', 'HadGEM2-ES'], fontsize=6, loc=9, ncol=2, frameon=False)
plt.text(1, 5.7, u'MBE = {}'.format(-1.6), fontsize=6, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(-0.8), fontsize=6, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(0.3), fontsize=6, color='gray')

ax20 = fig.add_subplot(3, 2, 2)
annual_cycle2 = ax20.plot(time, bias_tas_reg_exp1_cru_samz, time, bias_tas_reg_exp2_cru_samz, time, bias_tas_had_cru_samz)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
plt.setp(ax20.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle2
plt.setp(l1, linewidth=1.5, color='blue', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.text(1, 5.7, u'MBE = {}'.format(0.9), fontsize=6, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(-1.5), fontsize=6, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(-1.0), fontsize=6, color='gray')

ax30 = fig.add_subplot(3, 2, 3)
annual_cycle3 = ax30.plot(time, bias_pre_reg_exp1_cru_eneb, time, bias_pre_reg_exp2_cru_eneb, time, bias_pre_had_cru_eneb)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation (mm d⁻¹)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
plt.setp(ax30.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle3
plt.setp(l1, linewidth=1.5, color='blue', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.text(1, 5.7, u'MBE = {}'.format(-0.8), fontsize=6, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(0.2), fontsize=6, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(-0.5), fontsize=6, color='gray')
 
ax40 = fig.add_subplot(3, 2, 4)
annual_cycle4 = ax40.plot(time, bias_tas_reg_exp1_cru_eneb, time, bias_tas_reg_exp2_cru_eneb, time, bias_tas_had_cru_eneb)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature (°C)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
plt.setp(ax40.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle4
plt.setp(l1, linewidth=1.5, color='blue', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.text(1, 5.7, u'MBE = {}'.format(1.1), fontsize=6, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(-0.1), fontsize=6, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(-0.9), fontsize=6, color='gray')
	 
ax50 = fig.add_subplot(3, 2, 5)
annual_cycle5 = ax50.plot(time, bias_pre_reg_exp1_cru_matopiba, time, bias_pre_reg_exp2_cru_matopiba, time, bias_pre_had_cru_matopiba)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Months', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
l1, l2, l3 = annual_cycle5
plt.setp(l1, linewidth=1.5, color='blue', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=5, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.text(1, 5.7, u'MBE = {}'.format(-1.8), fontsize=6, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(0.5), fontsize=6, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(0.3), fontsize=6, color='gray')

ax60 = fig.add_subplot(3, 2, 6)
annual_cycle6 = ax60.plot(time, bias_tas_reg_exp1_cru_matopiba, time, bias_tas_reg_exp2_cru_matopiba, time, bias_tas_had_cru_matopiba)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Months', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-5, 5)
plt.yticks(np.arange(-5, 6, 2), fontsize=8)
l1, l2, l3 = annual_cycle6
plt.setp(l1, linewidth=1.5, color='blue', markersize=6, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l2, linewidth=1.5, color='red', markersize=6, marker='.', markerfacecolor='white', linestyle='--')
plt.setp(l3, linewidth=1.5, color='gray', markersize=6, marker='.', markerfacecolor='white', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axhline(-2, linewidth=1., linestyle='-', color='black')
plt.axhline(2, linewidth=1., linestyle='-', color='black')
plt.text(1, 5.7, u'MBE = {}'.format(2.0), fontsize=7, color='blue')
plt.text(5, 5.7, u'MBE = {}'.format(-0.2), fontsize=7, color='red')
plt.text(9, 5.7, u'MBE = {}'.format(-1.3), fontsize=7, color='gray')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()
