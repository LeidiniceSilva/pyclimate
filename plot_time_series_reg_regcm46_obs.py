# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "09/22/2018"
__description__ = "This script plot mensal and seasonal simulation and obs data"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import shiftgrid

from datetime import datetime, date
from PyFuncemeClimateTools import PlotMaps as pm
from hidropy.utils.hidropy_utils import create_path
from matplotlib.font_manager import FontProperties


def import_exp_model(exp):
	
	sim = []
	for mon in range(0,12):

		mod_path = '/home/nice'
		arq1     = '{0}/pre_amz_neb_regcm_{1}_2005_monmean.nc'.format(mod_path, exp)
		data1    = netCDF4.Dataset(arq1)
		var1     = data1.variables['pr'][:]
		lat      = data1.variables['lat'][:]
		lon      = data1.variables['lon'][:]
		mod      = var1[mon::12][:,:,:]
		mod_ini1 = np.nanmean(mod, axis=1)
		mod_end1 = np.nanmean(mod_ini1, axis=1)
		sim.append(mod_end1)
	return np.squeeze(sim)
	

def import_obs_data(obs):

	vars_dic = {u'cmap': u'precip', u'trmm': u'r'}
	
	data = []
	for mon in range(0,12):

		pre_path = '/home/nice'
		arq2     = '{0}/pre_amz_neb_{1}_obs_2005_monmean.nc'.format(pre_path, obs)
		data2    = netCDF4.Dataset(arq2)
		var2     = data2.variables[vars_dic[obs]][:]
		lat      = data2.variables['lat'][:]
		lon      = data2.variables['lon'][:]
		pre      = var2[mon::12][:,:,:]
		pre_ini1 = np.nanmean(pre, axis=1)
		pre_end1 = np.nanmean(pre_ini1, axis=1)
		data.append(pre_end1)

	return np.squeeze(data)
	

# Import simulations experiments and observed databases
exp  = 'exp1'
exp1 = import_exp_model(exp)

exp  = 'exp2'
exp2 = import_exp_model(exp)

obs  = 'cmap'
cmap = import_obs_data(obs)

obs  = 'trmm'
trmm = import_obs_data(obs)

print "exp1", exp1
print
print "exp2", exp2
print
print "cmap", cmap
print
print "trmm", trmm

# Plot figure of time series per region
fig = plt.figure(figsize=(12,6))
time = np.arange(0, 12)

a1 = plt.plot(time, exp1, time, exp2, time, cmap, time, trmm)

l1, l2, l3, l4 = a1
plt.setp(l1,  linewidth=2, markeredgewidth=1, marker='o', color='blue')
plt.setp(l2,  linewidth=2, markeredgewidth=1, marker='o', color='green')
plt.setp(l3,  linewidth=2, markeredgewidth=1, marker='o', color='red')
plt.setp(l4,  linewidth=2, markeredgewidth=1, marker='o', color='black')
				     
plt.title(u'Precipitação Média - AMZ_NEB (A1) - 2005', fontsize=16, fontweight='bold')

plt.xlabel(u'Meses', fontsize=16, fontweight='bold')
plt.ylabel(u'Precipitação', fontsize=16, fontweight='bold')

plt.ylim([0,40])

objects = ['JAN', 'FEV', 'MAR', 'ABR', 'MAI', 'JUN', 'JUL', 'AGO', 'SET', 'OUT', 'NOV', 'DEZ']
plt.xticks(time, objects, fontsize=12)
plt.grid(True, which='major', linestyle='--', linewidth='0.5', zorder=0.5)

font = FontProperties(size=10)
plt.legend([u'RegCM_EXP1', U'RegCM_EXP2', u'CMAP', u'TRMM'], loc='best', ncol=2, prop=font)

path_out = '/home/nice/'

graph_ts = 'pre_serie_temp_neb_A1_2005.png'
plt.savefig(os.path.join(path_out, graph_ts), bbox_inches='tight')
plt.show()
raise SystemExit


