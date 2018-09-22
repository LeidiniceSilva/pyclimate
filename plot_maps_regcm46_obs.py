# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "09/22/2018"
__description__ = "This script plot mensal and seasonal simulation and obs data"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import shiftgrid

from datetime import datetime, date
from PyFuncemeClimateTools import PlotMaps as pm
from hidropy.utils.hidropy_utils import create_path
from matplotlib.font_manager import FontProperties
#from pltsst import plotmap


sim = []
obs = []

for mon in range(0,12):

	# Open mod data	 
	mod_path = '/home/nice'
	arq1     = '{0}/pre_amz_neb_regcm_exp1_2005_monmean.nc'.format(mod_path)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	mod      = var1[mon::12][:,:,:]
	mod_ini1 = np.nanmean(mod, axis=1)
	mod_end1 = np.nanmean(mod_ini1, axis=1)
	sim.append(mod_end1)
	
	# Open pre data		 
	pre_path = '/home/nice'
	arq2     = '{0}/pre_amz_neb_cmap_obs_2005_monmean.nc'.format(pre_path)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables['precip'][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[mon::12][:,:,:]
	pre_ini1 = np.nanmean(pre, axis=1)
	pre_end1 = np.nanmean(pre_ini1, axis=1)
	obs.append(pre_end1)
	

data_sim = np.squeeze(sim)
data_obs =  np.squeeze(obs)

print "sim", data_sim
print
print "obs", data_obs

# Plot figure of time series per region
fig = plt.figure(figsize=(12,6))
time = np.arange(0, 12)

a1 = plt.plot(time, data_sim, time, data_obs)

plt.title(u'Precipitação média - AMZ_NEB - 2005', fontsize=16, fontweight='bold')

plt.xlabel(u'Meses', fontsize=16, fontweight='bold')
plt.ylabel(u'Precipitação', fontsize=16, fontweight='bold')

plt.ylim([0,40])

objects = ['JAN', 'FEV', 'MAR', 'ABR', 'MAI', 'JUN', 'JUL', 'AGO', 'SET', 'OUT', 'NOV', 'DEZ']
plt.xticks(time, objects, fontsize=12)
plt.grid(True, which='major', linestyle='-.', linewidth='0.5', zorder=0.5)

font = FontProperties(size=10)
plt.legend([u'SIM', u'OBS'], loc='best', ncol=2, prop=font)

path_out = '/home/nice/'

graph_ts = 'pre_med_mensal_amz_neb_2005.png'
plt.savefig(os.path.join(path_out, graph_ts), bbox_inches='tight')
plt.close()
raise SystemExit

exit()

