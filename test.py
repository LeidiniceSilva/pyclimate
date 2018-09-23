# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "05/26/2018"
__description__ = "This script plot mensal and seasonal anomaly simulation"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
# mpl.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
from matplotlib.font_manager import FontProperties
from PyFuncemeClimateTools import PlotMaps as pm
from datetime import datetime, date

exp = 'exp1'
obs = 'cmap'

for season in range(0,4):
	
	season_dic = {0: u'DJF 2004/2005', 1: u'MAM 2005', 2: u'JJA 2005', 3: u'SON 2005'}
	print season_dic[season]

	# Open exp model 			   
	mod_path = '/home/nice'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2004-2005_yseasmean.nc'.format(mod_path, exp)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	sim      = var1[season,:,:]
		
	# Open obs model 			   
	vars_dic = {u'cmap': u'precip', u'trmm': u'r'}

	pre_path = '/home/nice'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2004-2005_yseasmean.nc'.format(pre_path, obs)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dic[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[season,:,:]
	
	# RegCM4.6_EXP1 plotmaps	
	
	fig = plt.figure(figsize=(12, 8))
	
	fig_title = u'Precipitação Sazonal - RegCM4.6 - DJF 2004/2005'
	
	pm.plotmap(sim, lat, lon, fig_title=fig_title, barinf='both', barloc='right')

	plt.show()
	exit()
		   

		   

		   
	
