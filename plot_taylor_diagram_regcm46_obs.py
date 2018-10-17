# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "09/22/2018"
__description__ = "This script plot box and time series simulation and obs databases per region"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
#mpl.use('Agg')
import matplotlib.pyplot as plt
import skill_metrics as sm

from mpl_toolkits.basemap import shiftgrid
from datetime import datetime, date
from matplotlib.font_manager import FontProperties
from comp_statist_indices import compute_corr, compute_rmse


def import_exp_model(exp):
	
	mod_path = '/home/nice'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2004-2005_yseasmean.nc'.format(mod_path, exp)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	mod      = var1[:][:,:,:]
	mod_ini1 = np.nanmean(mod, axis=1)
	mod_end1 = np.nanmean(mod_ini1, axis=1)
	
	return np.squeeze(mod_end1)
	

def import_obs_data(obs):

	vars_dict = {u'cmap': u'precip', u'trmm': u'r'}

	pre_path = '/home/nice'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2004-2005_yseasmean.nc'.format(pre_path, obs)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dict[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[:][:,:,:]
	pre_ini1 = np.nanmean(pre, axis=1)
	pre_end1 = np.nanmean(pre_ini1, axis=1)

	return np.squeeze(pre_end1)
sdev

# Import simulations experiments and observed databases
exp  = 'exp1'
exp1 = import_exp_model(exp)

obs  = 'cmap'
cmap = import_obs_data(obs)


print "exp1", exp1
print

print "cmap", cmap
print

taylor_stats1 = sm.taylor_statistics(exp1,cmap,'data')
taylor_stats2 = sm.taylor_statistics(exp1,cmap,'data')
taylor_stats3 = sm.taylor_statistics(exp1,cmap,'data')

# Store statistics in arrays
sdev = np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1], taylor_stats2['sdev'][1], taylor_stats3['sdev'][1]])
crmsd = np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1], taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1]])
ccoef = np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1], taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1]])

print sdev
print crmsd
print ccoef
exit()


label = ['Non-Dimensional Observation', 'M1', 'M2', 'M3']

sm.taylor_diagram(sdev, crmsd, ccoef, markerLabel = label, styleOBS = '-', colOBS = 'r', 
markerobs = 'o', titleOBS = 'observation')

plt.savefig('taylor1.png')
plt.show()
exit()
