# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot trends from extremes index"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from matplotlib.font_manager import FontProperties


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)		

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 
	
	return obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 

	return rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 

	return gcm


def import_function_trend(data):
	
	data_x = [i for i in range(0, len(data))]
	data_x = np.reshape(data_x, (len(data_x), 1))
	data_y = data
	model = LinearRegression()
	model.fit(data_x, data_y)
	data_trend = model.predict(data_x)
	data_z = np.polyfit(data_y.flatten(), data_trend.flatten(), 1)
	data_r2 = r2_score(data_y, data_trend)
	data_median = statistics.median(data_y)

	return data_y, data_trend, data_median
	

# Import regcm exp and cru databases 	
obs_cdd = import_obs('eca_cdd', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_cdd = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cdd = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_cwd = import_obs('eca_cwd', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_cwd = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cwd = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_r10mm = import_obs('eca_r10mm', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_r10mm = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r10mm = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_r20mm = import_obs('eca_r20mm', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_r20mm = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r20mm = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Import function trend 	 
obs_cdd, obs_cdd_trend, obs_cdd_median = import_function_trend(obs_cdd)
obs_cwd, obs_cwd_trend, obs_cwd_median = import_function_trend(obs_cwd)
obs_r10mm, obs_r10mm_trend, obs_r10mm_median = import_function_trend(obs_r10mm)
obs_r20mm, obs_r20mm_trend, obs_r20mm_median = import_function_trend(obs_r20mm)

rcm_cdd, rcm_cdd_trend, rcm_cdd_median = import_function_trend(rcm_cdd)
rcm_cwd, rcm_cwd_trend, rcm_cwd_median = import_function_trend(rcm_cwd)
rcm_r10mm, rcm_r10mm_trend, rcm_r10mm_median = import_function_trend(rcm_r10mm)
rcm_r20mm, rcm_r20mm_trend, rcm_r20mm_median = import_function_trend(rcm_r20mm)

gcm_cdd, gcm_cdd_trend, gcm_cdd_median = import_function_trend(gcm_cdd)
gcm_cwd, gcm_cwd_trend, gcm_cwd_median = import_function_trend(gcm_cwd)
gcm_r20mm, gcm_r10mm_trend, gcm_r10mm_median = import_function_trend(gcm_r10mm)
gcm_r20mm, gcm_r20mm_trend, gcm_r20mm_median = import_function_trend(gcm_r20mm)

# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 20 + 1)
objects = [u'1986', u'1990', u'1994', u'1998', u'2002']

ax = fig.add_subplot(4, 3, 1)
ax.plot(time, obs_cdd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_cdd_trend, linewidth=1., color='gray')
ax.axhline(obs_cdd_median, linewidth=1., linestyle='dashed', color='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', direction='in', length=2, which='both', top=False, labelbottom=False)
plt.title(u'A) Obs y={0}x+{1}'.format(round(obs_cdd_median,2),round(obs_cdd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'CDD \n (days)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 50, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 2)
ax.plot(time, rcm_cdd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_cdd_trend, linewidth=1., color='gray')
ax.axhline(rcm_cdd_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'B) Reg y={0}x+{1}'.format(round(rcm_cdd_median,2),round(rcm_cdd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 50, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 3)
ax.plot(time, gcm_cdd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_cdd_trend, linewidth=1., color='gray')
ax.axhline(gcm_cdd_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'C) Had y={0}x+{1}'.format(round(gcm_cdd_median,2),round(gcm_cdd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 50, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 4)
ax.plot(time, obs_cwd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_cwd_trend, linewidth=1., color='gray')
ax.axhline(obs_cwd_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title('D) Obs y={0}x+{1}'.format(round(obs_cwd_median,2),round(obs_cwd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'CWD \n (days)', fontsize=7,  fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 35, 7), fontsize=7)

ax = fig.add_subplot(4, 3, 5)
ax.plot(time, rcm_cwd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_cwd_trend, linewidth=1., color='gray')
ax.axhline(rcm_cwd_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'E) Reg y={0}x+{1}'.format(round(rcm_cwd_median,2),round(rcm_cwd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 35, 7), fontsize=7)

ax = fig.add_subplot(4, 3, 6)
ax.plot(time, gcm_cwd, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_cwd_trend, linewidth=1., color='gray')
ax.axhline(gcm_cwd_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'F) Had y={0}x+{1}'.format(round(gcm_cwd_median,2),round(gcm_cwd_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 35, 7), fontsize=7)

ax = fig.add_subplot(4, 3, 7)
ax.plot(time, obs_r10mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_r10mm_trend, linewidth=1., color='gray')
ax.axhline(obs_r10mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'G) Obs y={0}x+{1}'.format(round(obs_r10mm_median,2),round(obs_r10mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'R10mm \n (days)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 8)
ax.plot(time, rcm_r10mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_r10mm_trend, linewidth=1., color='gray')
ax.axhline(rcm_r10mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'H) Reg y={0}x+{1}'.format(round(rcm_r10mm_median,2),round(rcm_r10mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 9)
ax.plot(time, gcm_r10mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_r10mm_trend, linewidth=1., color='gray')
ax.axhline(gcm_r10mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'I) Had y={0}x+{1}'.format(round(gcm_r10mm_median,2),round(gcm_r10mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(4, 3, 10)
ax.plot(time, obs_r20mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_r20mm_trend, linewidth=1., color='gray')
ax.axhline(obs_r20mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'J) Obs y={0}x+{1}'.format(round(obs_r20mm_median,2),round(obs_r20mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.ylabel(u'R20mm \n (days)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 35, 7), fontsize=7)

ax = fig.add_subplot(4, 3, 11)
ax.plot(time, rcm_r20mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_r20mm_trend, linewidth=1., color='gray')
ax.axhline(rcm_r20mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'K) Obs y={0}x+{1}'.format(round(rcm_r20mm_median,2),round(rcm_r20mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 35, 7), fontsize=7)

ax = fig.add_subplot(4, 3, 12)
ax.plot(time, gcm_r20mm, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_r20mm_trend, linewidth=1., color='gray')
ax.axhline(gcm_r20mm_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'L) Had y={0}x+{1}'.format(round(gcm_r20mm_median,2),round(gcm_r20mm_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 35, 7), fontsize=7)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.70)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_clim_etccdi_pre_count_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






