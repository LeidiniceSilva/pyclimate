# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute CDF functions from Reg and Had models"

import os
import netCDF4
import numpy as np
import numpy as np
import matplotlib.pyplot as plt


def cdf_function(data):
	
	sortedtime = np.sort(data)
	cdf = 1. * np.arange(len(data))/(len(data) - 1)
	
	return sortedtime, cdf
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return gcm

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return obs
	
	
# Import regcm exps model end obs database climatology
pre_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
pre_udel_samz = import_obs('pre', 'samz', 'udel_v301', '1986-2005')
pre_chirps_samz = import_obs('precip', 'samz', 'chirps-v2.0', '1986-2005')
pre_era5_samz = import_obs('mtpr', 'samz', 'era5', '1986-2005')

pre_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
pre_udel_eneb = import_obs('pre', 'eneb', 'udel_v301', '1986-2005')
pre_chirps_eneb = import_obs('precip', 'eneb', 'chirps-v2.0', '1986-2005')
pre_era5_eneb = import_obs('mtpr', 'eneb', 'era5', '1986-2005')

pre_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
pre_udel_matopiba = import_obs('pre', 'matopiba', 'udel_v301', '1986-2005')
pre_chirps_matopiba = import_obs('precip', 'matopiba', 'chirps-v2.0', '1986-2005')
pre_era5_matopiba = import_obs('mtpr', 'matopiba', 'era5', '1986-2005')

tas_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
tas_udel_samz = import_obs('temp', 'samz', 'udel_v301', '1986-2005')
tas_era5_samz = import_obs('t2m', 'samz', 'era5', '1986-2005')

tas_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
tas_udel_eneb = import_obs('temp', 'eneb', 'udel_v301', '1986-2005')
tas_era5_eneb = import_obs('t2m', 'eneb', 'era5', '1986-2005')

tas_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
tas_udel_matopiba = import_obs('temp', 'matopiba', 'udel_v301', '1986-2005')
tas_era5_matopiba = import_obs('t2m', 'matopiba', 'era5', '1986-2005')

# Import cdf function
sort_pre_reg_samz, cdf_pre_reg_samz = cdf_function(pre_reg_samz)
sort_pre_had_samz, cdf_pre_had_samz = cdf_function(pre_had_samz)
sort_pre_cru_samz, cdf_pre_cru_samz = cdf_function(pre_cru_samz)
sort_pre_udel_samz, cdf_pre_udel_samz = cdf_function(pre_udel_samz)
sort_pre_chirps_samz, cdf_pre_chirps_samz = cdf_function(pre_chirps_samz)
sort_pre_era5_samz, cdf_pre_era5_samz = cdf_function(pre_era5_samz)

sort_pre_reg_eneb, cdf_pre_reg_eneb = cdf_function(pre_reg_eneb)
sort_pre_had_eneb, cdf_pre_had_eneb = cdf_function(pre_had_eneb)
sort_pre_cru_eneb, cdf_pre_cru_eneb = cdf_function(pre_cru_eneb)
sort_pre_udel_eneb, cdf_pre_udel_eneb = cdf_function(pre_udel_eneb)
sort_pre_chirps_eneb, cdf_pre_chirps_eneb = cdf_function(pre_chirps_eneb)
sort_pre_era5_eneb, cdf_pre_era5_eneb = cdf_function(pre_era5_eneb)

sort_pre_reg_matopiba, cdf_pre_reg_matopiba = cdf_function(pre_reg_matopiba)
sort_pre_had_matopiba, cdf_pre_had_matopiba = cdf_function(pre_had_matopiba)
sort_pre_cru_matopiba, cdf_pre_cru_matopiba = cdf_function(pre_cru_matopiba)
sort_pre_udel_matopiba, cdf_pre_udel_matopiba = cdf_function(pre_udel_matopiba)
sort_pre_chirps_matopiba, cdf_pre_chirps_matopiba = cdf_function(pre_chirps_matopiba)
sort_pre_era5_matopiba, cdf_pre_era5_matopiba = cdf_function(pre_era5_matopiba)

sort_tas_reg_samz, cdf_tas_reg_samz = cdf_function(np.nanmean(tas_reg_samz, axis=0))
sort_tas_had_samz, cdf_tas_had_samz = cdf_function(tas_had_samz)
sort_tas_cru_samz, cdf_tas_cru_samz = cdf_function(tas_cru_samz)
sort_tas_udel_samz, cdf_tas_udel_samz = cdf_function(tas_udel_samz)
sort_tas_era5_samz, cdf_tas_era5_samz = cdf_function(tas_era5_samz)

sort_tas_reg_eneb, cdf_tas_reg_eneb = cdf_function(np.nanmean(tas_reg_eneb, axis=0))
sort_tas_had_eneb, cdf_tas_had_eneb = cdf_function(tas_had_eneb)
sort_tas_cru_eneb, cdf_tas_cru_eneb = cdf_function(tas_cru_eneb)
sort_tas_udel_eneb, cdf_tas_udel_eneb = cdf_function(tas_udel_eneb)
sort_tas_era5_eneb, cdf_tas_era5_eneb = cdf_function(tas_era5_eneb)

sort_tas_reg_matopiba, cdf_tas_reg_matopiba = cdf_function(np.nanmean(tas_reg_matopiba, axis=0))
sort_tas_had_matopiba, cdf_tas_had_matopiba = cdf_function(tas_had_matopiba)
sort_tas_cru_matopiba, cdf_tas_cru_matopiba = cdf_function(tas_cru_matopiba)
sort_tas_udel_matopiba, cdf_tas_udel_matopiba = cdf_function(tas_udel_matopiba)
sort_tas_era5_matopiba, cdf_tas_era5_matopiba = cdf_function(tas_era5_matopiba)

# Plot model end obs data climatology
fig = plt.figure()

ax1 = fig.add_subplot(3, 2, 1)
plt.plot(sort_pre_reg_samz, cdf_pre_reg_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_had_samz, cdf_pre_had_samz, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_cru_samz, cdf_pre_cru_samz, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_udel_samz, cdf_pre_udel_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_chirps_samz, cdf_pre_chirps_samz, color='yellow', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_era5_samz, cdf_pre_era5_samz, color='magenta', linestyle='--', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 2, 2)
plt.plot(sort_tas_reg_samz, cdf_tas_reg_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_samz, cdf_tas_had_samz, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_samz, cdf_tas_cru_samz, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_samz, cdf_tas_udel_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_samz, cdf_tas_era5_samz, color='magenta', linestyle='--', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B) SAMZ', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = fig.add_subplot(3, 2, 3)
plt.plot(sort_pre_reg_eneb, cdf_pre_reg_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_had_eneb, cdf_pre_had_eneb, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_cru_eneb, cdf_pre_cru_eneb, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_udel_eneb, cdf_pre_udel_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_chirps_eneb, cdf_pre_chirps_eneb, color='yellow', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_era5_eneb, cdf_pre_era5_eneb, color='magenta', linestyle='--', linewidth=1.5)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C) ENEB', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax4 = fig.add_subplot(3, 2, 4)
plt.plot(sort_tas_reg_eneb, cdf_tas_reg_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_eneb, cdf_tas_had_eneb, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_eneb, cdf_tas_cru_eneb, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_eneb, cdf_tas_udel_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_eneb, cdf_tas_era5_eneb, color='magenta', linestyle='--', linewidth=1.5)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'D) ENEB', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax4.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax5 = fig.add_subplot(3, 2, 5)
ax51 = plt.plot(sort_pre_reg_matopiba, cdf_pre_reg_matopiba, 
sort_pre_had_matopiba, cdf_pre_had_matopiba, 
sort_pre_cru_matopiba, cdf_pre_cru_matopiba, 
sort_pre_udel_matopiba, cdf_pre_udel_matopiba, 
sort_pre_chirps_matopiba, cdf_pre_chirps_matopiba, 
sort_pre_era5_matopiba, cdf_pre_era5_matopiba)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'E) MATOPIBA', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

l1, l2, l3, l4, l5, l6 = ax51
plt.setp(l1, color='black', linestyle='--', linewidth=1.5)
plt.setp(l2, color='blue', linestyle='--', linewidth=1.5)
plt.setp(l3, color='red', linestyle='--', linewidth=1.5)
plt.setp(l4, color='green', linestyle='--', linewidth=1.5)
plt.setp(l5, color='yellow', linestyle='--', linewidth=1.5)
plt.setp(l6, color='magenta', linestyle='--', linewidth=1.5)

legend = ['Reg','Had','CRU','UDEL','CHIRPS','ERA5']
plt.legend(ax51, legend, loc='lower left', bbox_to_anchor=(-0.3, -0.8), shadow=True, ncol=6)

ax6 = fig.add_subplot(3, 2, 6)
plt.plot(sort_tas_reg_matopiba, cdf_tas_reg_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_matopiba, cdf_tas_had_matopiba, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_matopiba, cdf_tas_cru_matopiba, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_matopiba, cdf_tas_udel_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_matopiba, cdf_tas_era5_matopiba, color='magenta', linestyle='--', linewidth=1.5)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Temperature (°C)', fontweight='bold')

fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_cdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()

