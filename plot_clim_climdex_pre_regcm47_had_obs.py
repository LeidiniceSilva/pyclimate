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
	obs  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 
	
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
	rcm  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

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
objects = [u'1986',u'',u'1988',u'',u'1990',u'',u'1992',u'',u'1994',u'',u'1996',u'',u'1998',u'',u'2000',u'',u'2002',u'',u'2004',u'']

ax1 = fig.add_subplot(611)
ax1.plot(time, obs_prcptot, marker='.', lw=1., color='black')
ax1.plot(time, obs_prcptot_trend, lw=1., color='black')
ax1.axhline(obs_prcptot_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('A) prectot', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(1000, 2400, 400))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax2 = fig.add_subplot(612)
ax2.plot(time, obs_r95p, marker='.', lw=1., color='black')
ax2.plot(time, obs_r95p_trend, lw=1., color='black')
ax2.axhline(obs_r95p_median, lw=1., linestyle='dashed', color='black', alpha=2)
ax2.set_ylabel('Precipitation (mm)', fontweight='bold')
plt.title('B) r95p', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(25, 40, 5))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

ax3 = fig.add_subplot(613)
ax3.plot(time, obs_r99p, marker='.', lw=1., color='black')
ax3.plot(time, obs_r99p_trend, lw=1., color='black')
ax3.axhline(obs_r99p_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('C) r99p', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(8, 14, 2))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 8, u'R2={0} y={1}x+{2}'.format(round(r99p_r2,2),round(r99p_z[0],2),round(r99p_z[1],2)), fontweight='bold')

ax4 = fig.add_subplot(614)
ax4.plot(time, obs_rx1day, marker='.', lw=1., color='black')
ax4.plot(time, obs_rx1day_trend, lw=1., color='black')
ax4.axhline(obs_rx1day_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('D) rx1day', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(21, 25, 1))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax5 = fig.add_subplot(615)
ax5.plot(time, obs_rx5day, marker='.', lw=1., color='black')
ax5.plot(time, obs_rx5day_trend, lw=1., color='black')
ax5.axhline(obs_rx5day_median, lw=1., linestyle='dashed', color='black', alpha=2)
ax5.set_ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.title('E) rx5day', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(42, 62, 5))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax6 = fig.add_subplot(616)
ax6.plot(time, obs_sdii, marker='.', lw=1., color='black')
ax6.plot(time, obs_sdii_trend, lw=1., color='black')
ax6.axhline(obs_sdii_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('F) sdii', loc='left', fontsize=10, fontweight='bold')
ax6.set_xlabel('Years', fontweight='bold')
plt.xticks(np.arange(1, 20, 1), objects)
plt.yticks(np.arange(6, 14, 2))
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.60)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_climdex_pre_regcm47_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






