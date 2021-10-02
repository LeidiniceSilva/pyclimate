# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot annual climatology from regcm47 and hadgem models and obs database"

import os
import netCDF4
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2'
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
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
		
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
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
	
	return gcm


def compute_linear_regression(obs, sim):

	m, c, r, p, se1 = stats.linregress(obs, sim)
	eq = 'y = %2.2fx+%2.2f'%(m, c)
	r2 = 'R² = %1.2f'%(r)
	p = 'p-value = %1.2f'%(p)
   
	return r2
	
	        
# Import models and obs database 
mon_pre_cru_amz = import_obs('pre', 'amz', 'cru_ts4.04', '1986-2005')
mon_pre_gpcp_amz = import_obs('precip', 'amz', 'gpcp_v2.2', '1986-2005')
mon_pre_era5_amz = import_obs('mtpr', 'amz', 'era5', '1986-2005')
mon_pre_reg_amz = import_rcm('pr', 'amz', 'historical', '1986-2005')
mon_pre_had_amz = import_gcm('pr', 'amz', 'historical', '1986-2005')

mon_pre_cru_neb = import_obs('pre', 'neb', 'cru_ts4.04', '1986-2005')
mon_pre_gpcp_neb = import_obs('precip', 'neb', 'gpcp_v2.2', '1986-2005')
mon_pre_era5_neb = import_obs('mtpr', 'neb', 'era5', '1986-2005')
mon_pre_reg_neb = import_rcm('pr', 'neb', 'historical', '1986-2005')
mon_pre_had_neb = import_gcm('pr', 'neb', 'historical', '1986-2005')

mon_tas_cru_amz = import_obs('tmp', 'amz', 'cru_ts4.04', '1986-2005')
mon_tas_era5_amz = import_obs('t2m', 'amz', 'era5', '1986-2005')
mon_tas_reg_amz = import_rcm('tas', 'amz', 'historical', '1986-2005')
mon_tas_had_amz = import_gcm('tas', 'amz', 'historical', '1986-2005')

mon_tas_cru_neb = import_obs('tmp', 'neb', 'cru_ts4.04', '1986-2005')
mon_tas_era5_neb = import_obs('t2m', 'neb', 'era5', '1986-2005')
mon_tas_reg_neb = import_rcm('tas', 'neb', 'historical', '1986-2005')
mon_tas_had_neb = import_gcm('tas', 'neb', 'historical', '1986-2005')

# Import linear regression
r2_pre_reg_cru_amz  = compute_linear_regression(mon_pre_cru_amz,  mon_pre_reg_amz)
r2_pre_reg_gpcp_amz = compute_linear_regression(mon_pre_gpcp_amz, mon_pre_reg_amz)
r2_pre_reg_era5_amz = compute_linear_regression(mon_pre_era5_amz, mon_pre_reg_amz)
r2_pre_reg_cru_neb  = compute_linear_regression(mon_pre_cru_neb,  mon_pre_reg_neb)
r2_pre_reg_gpcp_neb = compute_linear_regression(mon_pre_gpcp_neb, mon_pre_reg_neb)
r2_pre_reg_era5_neb = compute_linear_regression(mon_pre_era5_neb, mon_pre_reg_neb)
r2_pre_had_cru_amz  = compute_linear_regression(mon_pre_cru_amz,  mon_pre_had_amz)
r2_pre_had_gpcp_amz = compute_linear_regression(mon_pre_gpcp_amz, mon_pre_had_amz)
r2_pre_had_era5_amz = compute_linear_regression(mon_pre_era5_amz, mon_pre_had_amz)
r2_pre_had_cru_neb  = compute_linear_regression(mon_pre_cru_neb,  mon_pre_had_neb)
r2_pre_had_gpcp_neb = compute_linear_regression(mon_pre_gpcp_neb, mon_pre_had_neb)
r2_pre_had_era5_neb = compute_linear_regression(mon_pre_era5_neb, mon_pre_had_neb)

r2_tas_reg_cru_amz  = compute_linear_regression(mon_tas_cru_amz,  np.nanmean(mon_tas_reg_amz, axis=1))
r2_tas_reg_era5_amz = compute_linear_regression(mon_tas_era5_amz, np.nanmean(mon_tas_reg_amz, axis=1))
r2_tas_reg_cru_neb  = compute_linear_regression(mon_tas_cru_neb,  np.nanmean(mon_tas_reg_neb, axis=1))
r2_tas_reg_era5_neb = compute_linear_regression(mon_tas_era5_neb, np.nanmean(mon_tas_reg_neb, axis=1))
r2_tas_had_cru_amz  = compute_linear_regression(mon_tas_cru_amz,  mon_tas_had_amz)
r2_tas_had_era5_amz = compute_linear_regression(mon_tas_era5_amz, mon_tas_had_amz)
r2_tas_had_cru_neb  = compute_linear_regression(mon_tas_cru_neb,  mon_tas_had_neb)
r2_tas_had_era5_neb = compute_linear_regression(mon_tas_era5_neb, mon_tas_had_neb)

print(r2_tas_reg_cru_amz)
print(r2_tas_reg_era5_amz)
print(r2_tas_had_cru_amz)
print(r2_tas_had_era5_amz)
print(r2_tas_reg_cru_neb)
print(r2_tas_reg_era5_neb)
print(r2_tas_had_cru_neb)
print(r2_tas_had_era5_neb)

# Plot models and obs database 
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(2, 2, 1)
annual_cycle = ax.plot(time, mon_pre_cru_amz,  linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label='CRU')
annual_cycle = ax.plot(time, mon_pre_gpcp_amz, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='gray', label='GPCP')
annual_cycle = ax.plot(time, mon_pre_era5_amz, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label='ERA5')
annual_cycle = ax.plot(time, mon_pre_reg_amz,  linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label='RegCM4')
annual_cycle = ax.plot(time, mon_pre_had_amz,  linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='red', label='HadGEM2-ES')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.text(0.5, 2.5, 'R²=0.84 (CRU)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 1.5, 'R²=0.82 (GPCP)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 0.5, 'R²=0.82 (ERA5)', fontsize=6, color='green', fontweight='bold')
plt.text(7.5, 2.5, 'R²=0.98 (CRU)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 1.5, 'R²=0.97 (GPCP)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 0.5, 'R²=0.96 (ERA5)', fontsize=6, color='red', fontweight='bold')
plt.legend(fontsize=6, loc=9, ncol=1)

ax = fig.add_subplot(2, 2, 2)
annual_cycle = ax.plot(time, mon_tas_cru_amz,  linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label='CRU')
annual_cycle = ax.plot(time, mon_tas_era5_amz, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label='ERA5')
annual_cycle = ax.plot(time, np.nanmean(mon_tas_reg_amz, axis=1),  linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label='RegCM4')
annual_cycle = ax.plot(time, mon_tas_had_amz,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='red', label='HadGEM2-ES')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(20, 34, 2), fontsize=8)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis='both', which='major', labelsize=8)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.text(0.5, 22, 'R²=0.85 (CRU)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 21, 'R²=0.76 (ERA5)', fontsize=6, color='green', fontweight='bold')
plt.text(7.5, 22, 'R²=0.89 (CRU)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 21, 'R²=0.68 (ERA5)', fontsize=6, color='red', fontweight='bold')
plt.legend(fontsize=6, loc=9, ncol=1)

ax = fig.add_subplot(2, 2, 3)
annual_cycle = ax.plot(time, mon_pre_cru_neb,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='black')
annual_cycle = ax.plot(time, mon_pre_gpcp_neb, linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='gray')
annual_cycle = ax.plot(time, mon_pre_era5_neb, linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='blue')
annual_cycle = ax.plot(time, mon_pre_reg_neb,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='green')
annual_cycle = ax.plot(time, mon_pre_had_neb,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='red')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Months', fontsize=8, fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.text(0.5, 10, 'R²=0.96 (CRU)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 9, 'R²=0.96 (GPCP)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 8, 'R²=0.96 (ERA5)', fontsize=6, color='green', fontweight='bold')
plt.text(7.5, 10, 'R²=0.94 (CRU)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 9, 'R²=0.95 (GPCP)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 8, 'R²=0.95 (ERA5)', fontsize=6, color='red', fontweight='bold')

ax = fig.add_subplot(2, 2, 4)
annual_cycle = ax.plot(time, mon_tas_cru_neb,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='black')
annual_cycle = ax.plot(time, mon_tas_era5_neb, linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='blue')
annual_cycle = ax.plot(time, np.nanmean(mon_tas_reg_neb, axis=1),  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='green')
annual_cycle = ax.plot(time, mon_tas_had_neb,  linewidth=1.5, markersize=5, marker='.', markerfacecolor='white', linestyle='--', color='red')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Months', fontsize=8, fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(20, 34, 2), fontsize=8)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis='both', which='major', labelsize=8)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.text(0.5, 22, 'R²=0.81 (CRU)', fontsize=6, color='green', fontweight='bold')
plt.text(0.5, 21, 'R²=0.86 (ERA5)', fontsize=6, color='green', fontweight='bold')
plt.text(7.5, 22, 'R²=0.96 (CRU)', fontsize=6, color='red', fontweight='bold')
plt.text(7.5, 21, 'R²=0.97 (ERA5)', fontsize=6, color='red', fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
