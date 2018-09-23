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
from mpl_toolkits.basemap import shiftgrid

from datetime import datetime, date
from PyFuncemeClimateTools import PlotMaps as pm
from hidropy.utils.hidropy_utils import create_path
from matplotlib.font_manager import FontProperties


def import_exp_model(exp, area):
	
	mod_path = '/home/nice/subareas'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2005_monmean_{2}.nc'.format(mod_path, exp, area)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	mod      = var1[:][:,:,:]
	mod_ini1 = np.nanmean(mod, axis=1)
	mod_end1 = np.nanmean(mod_ini1, axis=1)
	
	return np.squeeze(mod_end1)
	

def import_obs_data(obs, area):

	vars_dict = {u'cmap': u'precip', u'trmm': u'r'}

	pre_path = '/home/nice/subareas'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2005_monmean_{2}.nc'.format(pre_path, obs, area)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dict[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[:][:,:,:]
	pre_ini1 = np.nanmean(pre, axis=1)
	pre_end1 = np.nanmean(pre_ini1, axis=1)

	return np.squeeze(pre_end1)
	

area_list = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12']
	
for area in area_list:
	
	area_dict = {u'A1': u'A1_Sudoeste_AMZ', u'A2': u'A2_Norte_AMZ', u'A3': u'A3_Noroeste_AMZ', u'A4': u'A4_Nordeste_AMZ',
	u'A5': u'A5_Sudeste_AMZ', u'A6': u'A6_Andes_sul', u'A7': u'A7_Andes_norte', u'A8': u'A8_Litoral_leste_NEB',
	u'A9': u'A9_Semiárido_norte_NEB', u'A10': u'A10_Semiárido_sul_NEB', u'A11': u'A11_Litoral_sul_NEB',
	u'A12': u'A12_Noroeste_NEB'}
	
	print "area -----------------> ", area_dict[area]
	print 
	
	# Import simulations experiments and observed databases
	exp  = 'exp1'
	exp1 = import_exp_model(exp, area)

	exp  = 'exp2'
	exp2 = import_exp_model(exp, area)

	obs  = 'cmap'
	cmap = import_obs_data(obs, area)

	obs  = 'trmm'
	trmm = import_obs_data(obs, area)

	print "exp1", exp1
	print
	print "exp2", exp2
	print
	print "cmap", cmap
	print
	print "trmm", trmm
	print

	# Precipitation boxplot per region
	fig = plt.figure(figsize=(12,6))
	time = np.arange(1, 5)

	data = [exp1, exp2, cmap, trmm]

	a2 = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)

	plt.title(u'Boxplot de Precipitação - {0} - 2005'.format(area_dict[area]), fontsize=16, fontweight='bold')

	plt.xlabel(u'Experimentos e Observação', fontsize=16, fontweight='bold')
	plt.ylabel(u'Precipitação (mm/m)', fontsize=16, fontweight='bold')
	plt.ylim([0,20])

	objects = [u'Reg_EXP1', u'Reg_EXP2', u'CMAP', u'TRMM']
	plt.xticks(time, objects, fontsize=12)

	path_out = '/home/nice/subareas/'
	plt.savefig(os.path.join(path_out, 'precip_boxplot_amz_neb_{0}_2005.png'.format(area)), dpi=100)

	# Plot precipitation time series per region
	# fig = plt.figure(figsize=(12,6))
	# time = np.arange(0, 11+1)

	# a1 = plt.plot(time, exp1, time, exp2, time, cmap, time, trmm)

	# l1, l2, l3, l4 = a1
	# plt.setp(l1,  linewidth=2, markeredgewidth=4, marker='+', color='blue')
	# plt.setp(l2,  linewidth=2, markeredgewidth=4, marker='+', color='green')
	# plt.setp(l3,  linewidth=2, markeredgewidth=4, marker='+', color='red')
	# plt.setp(l4,  linewidth=2, markeredgewidth=4, marker='+', color='black')
						 
	# plt.title(u'Precipitação Média - AMZ_NEB (A1) - 2005', fontsize=16, fontweight='bold')

	# plt.xlabel(u'Meses', fontsize=16, fontweight='bold')
	# plt.ylabel(u'Precipitação (mm/m)', fontsize=16, fontweight='bold')
	# plt.ylim([0,40])

	# objects = [u'JAN', u'FEV', u'MAR', u'ABR', u'MAI', u'JUN', u'JUL', u'AGO', u'SET', u'OUT', u'NOV', u'DEZ']
	# plt.xticks(time, objects, fontsize=12)
	# plt.grid(True, which='major', linestyle='--', linewidth='0.5', zorder=0.5)

	# font = FontProperties(size=10)
	# plt.legend([u'RegCM4.6_EXP1', U'RegCM4.6_EXP2', u'CMAP', u'TRMM'], loc='best', ncol=2, prop=font)

	# path_out = '/home/nice/'
	# plt.savefig(os.path.join(path_out, 'precip_serie_temp_amz_neb_a1_2005.png'), bbox_inches='tight')
	# raise SystemExit

exit()
