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

	dict_var = {u'eca_prectot': u'prec', 
	u'eca_r95p': u'prec',
	u'eca_r99p': u'prec', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 
	
	return obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prectot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 

	return rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prectot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period'}
	
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
obs_prcptot = import_obs('eca_prectot', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_prcptot = import_rcm('eca_prectot', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_prcptot = import_gcm('eca_prectot', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_r95p = import_obs('eca_r95p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_r95p = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r95p = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_r99p = import_obs('eca_r99p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_r99p = import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r99p = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_rx1day = import_obs('eca_rx1day', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_rx1day = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx1day = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_rx5day = import_obs('eca_rx5day', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_rx5day = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx5day = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_sdii = import_obs('eca_sdii', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_sdii = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_sdii = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Import function trend 	 
obs_prcptot, obs_prcptot_trend, obs_prcptot_median = import_function_trend(obs_prcptot)
obs_r95p, obs_r95p_trend, obs_r95p_median = import_function_trend(obs_r95p)
obs_r99p, obs_r99p_trend, obs_r99p_median = import_function_trend(obs_r99p)
obs_rx1day, obs_rx1day_trend, obs_rx1day_median = import_function_trend(obs_rx1day)
obs_rx5day, obs_rx5day_trend, obs_rx5day_median = import_function_trend(obs_rx5day)
obs_sdii, obs_sdii_trend, obs_sdii_median = import_function_trend(obs_sdii)

rcm_prcptot, rcm_prcptot_trend, rcm_prcptot_median = import_function_trend(rcm_prcptot)
rcm_r95p, rcm_r95p_trend, rcm_r95p_median = import_function_trend(rcm_r95p)
rcm_r99p, rcm_r99p_trend, rcm_r99p_median = import_function_trend(rcm_r99p)
rcm_rx1day, rcm_rx1day_trend, rcm_rx1day_median = import_function_trend(rcm_rx1day)
rcm_rx5day, rcm_rx5day_trend, rcm_rx5day_median = import_function_trend(rcm_rx5day)
rcm_sdii, rcm_sdii_trend, rcm_sdii_median = import_function_trend(rcm_sdii)

gcm_prcptot, gcm_prcptot_trend, gcm_prcptot_median = import_function_trend(gcm_prcptot)
gcm_r95p, gcm_r95p_trend, gcm_r95p_median = import_function_trend(gcm_r95p)
gcm_r99p, gcm_r99p_trend, gcm_r99p_median = import_function_trend(gcm_r99p)
gcm_rx1day, gcm_rx1day_trend, gcm_rx1day_median = import_function_trend(gcm_rx1day)
gcm_rx5day, gcm_rx5day_trend, gcm_rx5day_median = import_function_trend(gcm_rx5day)
gcm_sdii, gcm_sdii_trend, gcm_sdii_median = import_function_trend(gcm_sdii)
	
# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 20 + 1)
objects = [u'1986', u'1990', u'1994', u'1998', u'2002']

ax = fig.add_subplot(6, 3, 1)
ax.plot(time, obs_prcptot, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_prcptot_trend, linewidth=1., color='gray')
ax.axhline(obs_prcptot_median, linewidth=1., linestyle='dashed', color='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', direction='in', length=2, which='both', top=False, labelbottom=False)
plt.title(u'A) Obs y={0}x+{1}'.format(round(obs_prcptot_median,2),round(obs_prcptot_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'PRCPTOT (mm)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1300, 2500, 300), fontsize=7)

ax = fig.add_subplot(6, 3, 2)
ax.plot(time, rcm_prcptot, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_prcptot_trend, linewidth=1., color='gray')
ax.axhline(rcm_prcptot_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'B) Reg y={0}x+{1}'.format(round(rcm_prcptot_median,2),round(rcm_prcptot_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1300, 2500, 300), fontsize=7)

ax = fig.add_subplot(6, 3, 3)
ax.plot(time, gcm_prcptot, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_prcptot_trend, linewidth=1., color='gray')
ax.axhline(gcm_prcptot_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'C) Had y={0}x+{1}'.format(round(gcm_prcptot_median,2),round(gcm_prcptot_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(1300, 2500, 300), fontsize=7)

ax = fig.add_subplot(6, 3, 4)
ax.plot(time, obs_r95p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_r95p_trend, linewidth=1., color='gray')
ax.axhline(obs_r95p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title('D) Obs y={0}x+{1}'.format(round(obs_r95p_median,2),round(obs_r95p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'R95p \n(mm)', fontsize=7,  fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 5)
ax.plot(time, rcm_r95p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_r95p_trend, linewidth=1., color='gray')
ax.axhline(rcm_r95p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'E) Reg y={0}x+{1}'.format(round(rcm_r95p_median,2),round(rcm_r95p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 6)
ax.plot(time, gcm_r95p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_r95p_trend, linewidth=1., color='gray')
ax.axhline(gcm_r95p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'F) Had y={0}x+{1}'.format(round(gcm_r95p_median,2),round(gcm_r95p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 70, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 7)
ax.plot(time, obs_r99p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_r99p_trend, linewidth=1., color='gray')
ax.axhline(obs_r99p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'G) Obs y={0}x+{1}'.format(round(obs_r99p_median,2),round(obs_r99p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'R99p \n(mm)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 45, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 8)
ax.plot(time, rcm_r99p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_r99p_trend, linewidth=1., color='gray')
ax.axhline(rcm_r99p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'H) Reg y={0}x+{1}'.format(round(rcm_r99p_median,2),round(rcm_r99p_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 45, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 9)
ax.plot(time, gcm_r99p, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_r99p_trend, linewidth=1., color='gray')
ax.axhline(gcm_r99p_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'I) Had y={0}x+{1}'.format(round(gcm_prcptot_median,2),round(gcm_prcptot_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(5, 45, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 10)
ax.plot(time, obs_rx1day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_rx1day_trend, linewidth=1., color='gray')
ax.axhline(obs_rx1day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'J) Obs y={0}x+{1}'.format(round(obs_rx1day_median,2),round(obs_rx1day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'Rx1day \n(mm d⁻¹)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 50, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 11)
ax.plot(time, rcm_rx1day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_rx1day_trend, linewidth=1., color='gray')
ax.axhline(rcm_rx1day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title('K) Reg y={0}x+{1}'.format(round(rcm_rx1day_median,2),round(rcm_rx1day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 50, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 12)
ax.plot(time, gcm_rx1day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_rx1day_trend, linewidth=1., color='gray')
ax.axhline(gcm_rx1day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'L) Had y={0}x+{1}'.format(round(gcm_rx1day_median,2),round(gcm_rx1day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 50, 10), fontsize=7)

ax = fig.add_subplot(6, 3, 13)
ax.plot(time, obs_rx5day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_rx5day_trend, linewidth=1., color='gray')
ax.axhline(obs_rx5day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'M) Obs y={0}x+{1}'.format(round(obs_rx5day_median,2),round(obs_rx5day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'Rx5day \n(mm d⁻¹)', fontsize=7, fontweight='bold')
plt.yticks(np.arange(20, 200, 40), fontsize=7)

ax = fig.add_subplot(6, 3, 14)
ax.plot(time, rcm_rx5day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_rx5day_trend, linewidth=1., color='gray')
ax.axhline(rcm_rx5day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'N) Reg y={0}x+{1}'.format(round(rcm_rx5day_median,2),round(rcm_rx5day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 200, 40), fontsize=7)

ax = fig.add_subplot(6, 3, 15)
ax.plot(time, gcm_rx5day, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_rx5day_trend, linewidth=1., color='gray')
ax.axhline(gcm_rx5day_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'O) Had y={0}x+{1}'.format(round(gcm_rx5day_median,2),round(gcm_rx5day_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(20, 200, 40), fontsize=7)

ax = fig.add_subplot(6, 3, 16)
ax.plot(time, obs_sdii, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_sdii_trend, linewidth=1., color='gray')
ax.axhline(obs_sdii_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'P) Obs y={0}x+{1}'.format(round(obs_sdii_median,2),round(obs_sdii_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.ylabel(u'SDII \n(mm d⁻¹)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(7, 13, 1), fontsize=7)

ax = fig.add_subplot(6, 3, 17)
ax.plot(time, rcm_sdii, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, rcm_sdii_trend, linewidth=1., color='gray')
ax.axhline(rcm_sdii_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'Q) Reg y={0}x+{1}'.format(round(rcm_sdii_median,2),round(rcm_sdii_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(7, 13, 1), fontsize=7)

ax = fig.add_subplot(6, 3, 18)
ax.plot(time, gcm_sdii, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_sdii_trend, linewidth=1., color='gray')
ax.axhline(gcm_sdii_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'R) Had y={0}x+{1}'.format(round(gcm_sdii_median,2),round(gcm_sdii_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(7, 13, 1), fontsize=7)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.70)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_clim_etccdi_pre_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






