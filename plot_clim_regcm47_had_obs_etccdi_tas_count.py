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

	dict_var = {u'eca_txx': u'Tmax', 
	u'eca_txn': u'Tmax',
	u'eca_tnx': u'Tmin', 
	u'eca_tnn': u'Tmin',
	u'eca_dtr': u'Tmax'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 
	
	return obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

	return rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

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
obs_txx = import_obs('eca_txx', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_txx = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txx = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_txn = import_obs('eca_txn', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_txn = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txn = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tnx = import_obs('eca_tnx', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tnx = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnx = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_tnn = import_obs('eca_tnn', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_tnn = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnn = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

obs_dtr = import_obs('eca_dtr', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
rcm_dtr = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_dtr = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Import function trend 	 
obs_txx, obs_txx_trend, obs_txx_median = import_function_trend(obs_txx)
obs_txn, obs_txn_trend, obs_txn_median = import_function_trend(obs_txn)
obs_tnx, obs_tnx_trend, obs_tnx_median = import_function_trend(obs_tnx)
obs_tnn, obs_tnn_trend, obs_tnn_median = import_function_trend(obs_tnn)
obs_dtr, obs_dtr_trend, obs_dtr_median = import_function_trend(obs_dtr)

#~ rcm_txx, rcm_txx_trend, rcm_txx_median = import_function_trend(rcm_txx)
#~ rcm_txn, rcm_txn_trend, rcm_txn_median = import_function_trend(rcm_txn)
#~ rcm_tnx, rcm_tnx_trend, rcm_tnx_median = import_function_trend(rcm_tnx)
#~ rcm_tnn, rcm_tnn_trend, rcm_tnn_median = import_function_trend(rcm_tnn)
#~ rcm_dtr, rcm_dtr_trend, rcm_dtr_median = import_function_trend(rcm_dtr)

gcm_txx, gcm_txx_trend, gcm_txx_median = import_function_trend(gcm_txx)
gcm_txn, gcm_txn_trend, gcm_txn_median = import_function_trend(gcm_txn)
gcm_tnx, gcm_tnx_trend, gcm_tnx_median = import_function_trend(gcm_tnx)
gcm_tnn, gcm_tnn_trend, gcm_tnn_median = import_function_trend(gcm_tnn)
gcm_dtr, gcm_dtr_trend, gcm_dtr_median = import_function_trend(gcm_dtr)
	
# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 20 + 1)
objects = [u'1986', u'1990', u'1994', u'1998', u'2002']

ax = fig.add_subplot(5, 3, 1)
ax.plot(time, obs_txx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_txx_trend, linewidth=1., color='gray')
ax.axhline(obs_txx_median, linewidth=1., linestyle='dashed', color='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', direction='in', length=2, which='both', top=False, labelbottom=False)
plt.title(u'A) Obs y={0}x+{1}'.format(round(obs_txx_median,2),round(obs_txx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'txx (°C)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(34, 42, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 2)
ax.plot(time, gcm_txx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_txx_trend, linewidth=1., color='gray')
ax.axhline(gcm_txx_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'B) Reg y={0}x+{1}'.format(round(gcm_txx_median,2),round(gcm_txx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(34, 42, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 3)
ax.plot(time, gcm_txx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_txx_trend, linewidth=1., color='gray')
ax.axhline(gcm_txx_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'C) Had y={0}x+{1}'.format(round(gcm_txx_median,2),round(gcm_txx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(34, 42, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 4)
ax.plot(time, obs_txn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_txn_trend, linewidth=1., color='gray')
ax.axhline(obs_txn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title('D) Obs y={0}x+{1}'.format(round(obs_txn_median,2),round(obs_txn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'txn (°C)', fontsize=7,  fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(22, 28, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 5)
ax.plot(time, gcm_txn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_txn_trend, linewidth=1., color='gray')
ax.axhline(gcm_txn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'E) Reg y={0}x+{1}'.format(round(gcm_txn_median,2),round(gcm_txn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(22, 28, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 6)
ax.plot(time, gcm_txn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_txn_trend, linewidth=1., color='gray')
ax.axhline(gcm_txn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'F) Had y={0}x+{1}'.format(round(gcm_txn_median,2),round(gcm_txn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
#~ plt.yticks(np.arange(20, 40, 5), fontsize=7)

ax = fig.add_subplot(5, 3, 7)
ax.plot(time, obs_tnx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tnx_trend, linewidth=1., color='gray')
ax.axhline(obs_tnx_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'G) Obs y={0}x+{1}'.format(round(obs_tnx_median,2),round(obs_tnx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'tnx (°C)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(22, 28, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 8)
ax.plot(time, gcm_tnx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tnx_trend, linewidth=1., color='gray')
ax.axhline(gcm_tnx_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'H) Reg y={0}x+{1}'.format(round(gcm_tnx_median,2),round(gcm_tnx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(22, 28, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 9)
ax.plot(time, gcm_tnx, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tnx_trend, linewidth=1., color='gray')
ax.axhline(gcm_tnx_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'I) Had y={0}x+{1}'.format(round(gcm_tnx_median,2),round(gcm_tnx_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(22, 28, 2), fontsize=7)

ax = fig.add_subplot(5, 3, 10)
ax.plot(time, obs_tnn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_tnn_trend, linewidth=1., color='gray')
ax.axhline(obs_tnn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'J) Obs y={0}x+{1}'.format(round(obs_tnn_median,2),round(obs_tnn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.ylabel(u'tnn (°C)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(12, 32, 10), fontsize=7)

ax = fig.add_subplot(5, 3, 11)
ax.plot(time, gcm_tnn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tnn_trend, linewidth=1., color='gray')
ax.axhline(gcm_tnn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title('K) Reg y={0}x+{1}'.format(round(gcm_tnn_median,2),round(gcm_tnn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(12, 32, 10), fontsize=7)

ax = fig.add_subplot(5, 3, 12)
ax.plot(time, gcm_tnn, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_tnn_trend, linewidth=1., color='gray')
ax.axhline(gcm_tnn_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2, which='both', top=False, labelbottom=False)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'L) Had y={0}x+{1}'.format(round(gcm_tnn_median,2),round(gcm_tnn_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(12, 32, 10), fontsize=7)

ax = fig.add_subplot(5, 3, 13)
ax.plot(time, obs_dtr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, obs_dtr_trend, linewidth=1., color='gray')
ax.axhline(obs_dtr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2)
plt.title(u'M) Obs y={0}x+{1}'.format(round(obs_dtr_median,2),round(obs_dtr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.ylabel(u'dtr (°C)', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 12, 1), fontsize=7)

ax = fig.add_subplot(5, 3, 14)
ax.plot(time, gcm_dtr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_dtr_trend, linewidth=1., color='gray')
ax.axhline(gcm_dtr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'N) Reg y={0}x+{1}'.format(round(gcm_dtr_median,2),round(gcm_dtr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 12, 1), fontsize=7)

ax = fig.add_subplot(5, 3, 15)
ax.plot(time, gcm_dtr, linestyle='', linewidth=1., color='black', marker='.', markerfacecolor='white', markersize=5)
ax.plot(time, gcm_dtr_trend, linewidth=1., color='gray')
ax.axhline(gcm_dtr_median, linewidth=1., linestyle='dashed', color='black', alpha=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', direction='in', length=2)
ax.tick_params(axis='y', direction='in', length=2, which='both', right=False, labelleft=False)
plt.title(u'O) Had y={0}x+{1}'.format(round(gcm_dtr_median,2),round(gcm_dtr_trend[1],2)), loc='left', fontsize=7, fontweight='bold')
plt.xlabel(u'Years', fontsize=7, fontweight='bold')
plt.xticks(np.arange(1, 20, 4), objects, fontsize=7)
plt.yticks(np.arange(10, 12, 1), fontsize=7)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.70)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_clim_etccdi_tas_reg_had_obs_1986-2005'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()
