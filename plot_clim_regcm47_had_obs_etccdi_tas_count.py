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

	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}


	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 
	
	return obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(np.nanmean(var[:][:20,:,:], axis=1), axis=1) 

	return rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}

	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(np.nanmean(var[:][20:,:,:], axis=1), axis=1) 

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
obs_su = import_obs('eca_su', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_su = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_su = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tr = import_obs('eca_tr', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tr = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tr = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tx10p = import_obs('eca_tx10p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tx10p = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx10p = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tx90p = import_obs('eca_tx90p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tx90p = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx90p = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tn10p = import_obs('eca_tn10p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tn10p = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn10p = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tn90p = import_obs('eca_tn90p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tn90p = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn90p = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Import function trend 	 
obs_su, obs_su_trend, obs_su_median = import_function_trend(obs_su)
obs_tr, obs_tr_trend, obs_tr_median = import_function_trend(obs_tr)
obs_tx10p, obs_tx10p_trend, obs_tx10p_median = import_function_trend(obs_tx10p)
obs_tx90p, obs_tx90p_trend, obs_tx90p_median = import_function_trend(obs_tx90p)
obs_tn10p, obs_tn10p_trend, obs_tn10p_median = import_function_trend(obs_tn10p)
obs_tn90p, obs_tn90p_trend, obs_tn90p_median = import_function_trend(obs_tn90p)

rcm_su, rcm_su_trend, rcm_su_median = import_function_trend(np.nanmean(rcm_su, axis=1))
rcm_tr, rcm_tr_trend, rcm_tr_median = import_function_trend(np.nanmean(rcm_tr, axis=1))
rcm_tx10p, rcm_tx10p_trend, rcm_tx10p_median = import_function_trend(np.nanmean(rcm_tx10p, axis=1))
rcm_tx90p, rcm_tx90p_trend, rcm_tx90p_median = import_function_trend(np.nanmean(rcm_tx90p, axis=1))
rcm_tn10p, rcm_tn10p_trend, rcm_tn10p_median = import_function_trend(np.nanmean(rcm_tn10p, axis=1))
rcm_tn90p, rcm_tn90p_trend, rcm_tn90p_median = import_function_trend(np.nanmean(rcm_tn90p, axis=1))

gcm_su, gcm_su_trend, gcm_su_median = import_function_trend(gcm_su)
gcm_tr, gcm_tr_trend, gcm_tr_median = import_function_trend(gcm_tr)
gcm_tx10p, gcm_tx10p_trend, gcm_tx10p_median = import_function_trend(gcm_tx10p)
gcm_tx90p, gcm_tx90p_trend, gcm_tx90p_median = import_function_trend(gcm_tx90p)
gcm_tn10p, gcm_tn10p_trend, gcm_tn10p_median = import_function_trend(gcm_tn10p)
gcm_tn90p, gcm_tn90p_trend, gcm_tn90p_median = import_function_trend(gcm_tn90p)

# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 20 + 1)
objects = [u'1986', u'1990', u'1994', u'1998', u'2002']

ax = fig.add_subplot(6, 3, 1)
ax.plot(time, obs_su, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_su_trend, linewidth=1., color='gray')
ax.axhline(obs_su_median, linewidth=1., linestyle='dashed', color='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', direction='in', length=2, which='both', top=False, labelbottom=False)
plt.title(u'A) Obs y={0}x+{1}'.format(round(obs_su_median,2),round(obs_su_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'SU (days)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 2)
ax.plot(time, rcm_su, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_su_trend, linewidth=1., color='gray')
ax.axhline(rcm_su_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'B) Reg y={0}x+{1}'.format(round(rcm_su_median,2),round(rcm_su_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 3)
ax.plot(time, gcm_su, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_su_trend, linewidth=1., color='gray')
ax.axhline(gcm_su_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'C) Had y={0}x+{1}'.format(round(gcm_su_median,2),round(gcm_su_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 4)
ax.plot(time, obs_tr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tr_trend, linewidth=1., color='gray')
ax.axhline(obs_tr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title('D) Obs y={0}x+{1}'.format(round(obs_tr_median,2),round(obs_tr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'TR (days)', fontsize=7,  fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 5)
ax.plot(time, rcm_tr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_tr_trend, linewidth=1., color='gray')
ax.axhline(rcm_tr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'E) Reg y={0}x+{1}'.format(round(rcm_tr_median,2),round(rcm_tr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 6)
ax.plot(time, gcm_tr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tr_trend, linewidth=1., color='gray')
ax.axhline(gcm_tr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'F) Had y={0}x+{1}'.format(round(gcm_tr_median,2),round(gcm_tr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(100, 300, 50), fontsize=7)

ax = fig.add_subplot(6, 3, 7)
ax.plot(time, obs_tx10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tx10p_trend, linewidth=1., color='gray')
ax.axhline(obs_tx10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'G) Obs y={0}x+{1}'.format(round(obs_tx10p_median,2),round(obs_tx10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'TX10p (%)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 8)
ax.plot(time, rcm_tx10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_tx10p_trend, linewidth=1., color='gray')
ax.axhline(rcm_tx10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'H) Reg y={0}x+{1}'.format(round(rcm_tx10p_median,2),round(rcm_tx10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 9)
ax.plot(time, gcm_tx10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tx10p_trend, linewidth=1., color='gray')
ax.axhline(gcm_tx10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'I) Had y={0}x+{1}'.format(round(gcm_tx10p_median,2),round(gcm_tx10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 10)
ax.plot(time, obs_tx90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tx90p_trend, linewidth=1., color='gray')
ax.axhline(obs_tx90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'J) Obs y={0}x+{1}'.format(round(obs_tx90p_median,2),round(obs_tx90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'TX90p (%)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 11)
ax.plot(time, rcm_tx90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_tx90p_trend, linewidth=1., color='gray')
ax.axhline(rcm_tx90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title('K) Reg y={0}x+{1}'.format(round(rcm_tx90p_median,2),round(rcm_tx90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 12)
ax.plot(time, gcm_tx90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tx90p_trend, linewidth=1., color='gray')
ax.axhline(gcm_tx90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'L) Had y={0}x+{1}'.format(round(gcm_tx90p_median,2),round(gcm_tx90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 13)
ax.plot(time, obs_tn10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tn10p_trend, linewidth=1., color='gray')
ax.axhline(obs_tn10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'M) Obs y={0}x+{1}'.format(round(obs_tn10p_median,2),round(obs_tn10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'TN10p (%)', fontsize=7, fontweight='bold')
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 14)
ax.plot(time, rcm_tn10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_tn10p_trend, linewidth=1., color='gray')
ax.axhline(rcm_tn10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'N) Reg y={0}x+{1}'.format(round(rcm_tn10p_median,2),round(rcm_tn10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 15)
ax.plot(time, gcm_tn10p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tn10p_trend, linewidth=1., color='gray')
ax.axhline(gcm_tn10p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'O) Had y={0}x+{1}'.format(round(gcm_tn10p_median,2),round(gcm_tn10p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1, 30, 5), fontsize=7)

ax = fig.add_subplot(6, 3, 16)
ax.plot(time, obs_tn90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tn90p_trend, linewidth=1., color='gray')
ax.axhline(obs_tn90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'P) Obs y={0}x+{1}'.format(round(obs_tn90p_median,2),round(obs_tn90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.ylabel(u'TN90p (%)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(0, 50, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 17)
ax.plot(time, rcm_tn90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_tn90p_trend, linewidth=1., color='gray')
ax.axhline(rcm_tn90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'Q) Reg y={0}x+{1}'.format(round(rcm_tn90p_median,2),round(rcm_tn90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(0, 50, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 18)
ax.plot(time, gcm_tn90p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tn90p_trend, linewidth=1., color='gray')
ax.axhline(gcm_tn90p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'R) Had y={0}x+{1}'.format(round(gcm_tn90p_median,2),round(gcm_tn90p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(0, 50, 10), fontsize=7)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.70)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_clim_etccdi_tas_count_reg_had_obs_1986-2005'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()
