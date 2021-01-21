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
	
	
#~ def import_rcm(var, area, model, exp, freq, dt):
	
	#~ path = '/home/nice/Documents/dataset/rcm/eca'	
	#~ arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	#~ dict_var = {u'eca_txx': u'Tmax', 
	#~ u'eca_txn': u'Tmax',
	#~ u'eca_tnx': u'Tmin', 
	#~ u'eca_tnn': u'Tmin',
	#~ u'eca_dtr': u'Tmax'}
	
	#~ data = netCDF4.Dataset(arq)
	#~ var  = data.variables[var][:]
	#~ lat  = data.variables['lat'][:]
	#~ lon  = data.variables['lon'][:]
	#~ rcm  = np.nanmean(np.nanmean(var[:-1,:,:], axis=1), axis=1) 

	#~ return rcm


#~ def import_gcm(var, area, freq,  model, exp, dt):
	
	#~ path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	#~ arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	

	#~ dict_var = {u'eca_txx': u'Tmax', 
	#~ u'eca_txn': u'Tmax',
	#~ u'eca_tnx': u'Tmin', 
	#~ u'eca_tnn': u'Tmin',
	#~ u'eca_dtr': u'Tmax'}
	
	#~ data = netCDF4.Dataset(arq)
	#~ var  = data.variables[dict_var[var]][:]
	#~ lat  = data.variables['lat'][:]
	#~ lon  = data.variables['lon'][:]
	#~ gcm  = np.nanmean(np.nanmean(var[:-1,:,:], axis=1), axis=1) 

	#~ return gcm
	

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
#~ rcm_prcptot = import_rcm('eca_txx', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
#~ gcm_prcptot = import_gcm('eca_txx', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_txn = import_obs('eca_txn', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
#~ rcm_r95p = import_rcm('eca_txn', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
#~ gcm_r95p = import_gcm('eca_txn', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_tnx = import_obs('eca_tnx', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
#~ rcm_rx1day = import_rcm('eca_tnx', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
#~ gcm_rx1day = import_gcm('eca_tnx', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_tnn = import_obs('eca_tnn', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
#~ rcm_rx5day = import_rcm('eca_tnn', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
#~ gcm_rx5day = import_gcm('eca_tnn', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_dtr = import_obs('eca_dtr', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
#~ rcm_r99p = import_rcm('eca_dtr', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
#~ gcm_r99p = import_gcm('eca_dtr', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

# Import function trend 	 
txx, txx_trend, txx_median = import_function_trend(obs_txx)
txn, txn_trend, txn_median = import_function_trend(obs_txn)
tnx, tnx_trend, tnx_median = import_function_trend(obs_tnx)
tnn, tnn_trend, tnn_median = import_function_trend(obs_tnn)
dtr, dtr_trend, dtr_median = import_function_trend(obs_dtr)
	
# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 20 + 1)
objects = [u'1986',u'',u'1988',u'',u'1990',u'',u'1992',u'',u'1994',u'',u'1996',u'',u'1998',u'',u'2000',u'',u'2002',u'',u'2004',u'']

ax1 = fig.add_subplot(511)
ax1.plot(time, txx, marker='.', lw=1., color='black')
ax1.plot(time, txx_trend, lw=1., color='black')
ax1.axhline(txx_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('A) txx', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(34, 39, 1))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax2 = fig.add_subplot(512)
ax2.plot(time, txn, marker='.', lw=1., color='black')
ax2.plot(time, txn_trend, lw=1., color='black')
ax2.axhline(txn_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('B) txn', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(22, 27, 1))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

ax3 = fig.add_subplot(513)
ax3.plot(time, tnx, marker='.', lw=1., color='black')
ax3.plot(time, tnx_trend, lw=1., color='black')
ax3.axhline(tnx_median, lw=1., linestyle='dashed', color='black', alpha=2)
ax3.set_ylabel('Temperature (°C)', fontweight='bold')
plt.title('C) tnx', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(22, 26, 1))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 8, u'R2={0} y={1}x+{2}'.format(round(r99p_r2,2),round(r99p_z[0],2),round(r99p_z[1],2)), fontweight='bold')

ax4 = fig.add_subplot(514)
ax4.plot(time, tnn, marker='.', lw=1., color='black')
ax4.plot(time, tnn_trend, lw=1., color='black')
ax4.axhline(tnn_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('D) tnn', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(13, 18, 1))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax5 = fig.add_subplot(515)
ax5.plot(time, dtr, marker='.', lw=1., color='black')
ax5.plot(time, dtr_trend, lw=1., color='black')
ax5.axhline(dtr_median, lw=1., linestyle='dashed', color='black', alpha=2)
plt.title('E) dtr', loc='left', fontsize=10, fontweight='bold')
ax5.set_xlabel('Years', fontweight='bold')
plt.xticks(np.arange(1, 20, 1), objects)
plt.yticks(np.arange(9, 13, 1))
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.60)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_climdex_tas_regcm47_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






