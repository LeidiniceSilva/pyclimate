# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "05/26/2018"
__description__ = "Compute Bias and write statistic indice"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
# mpl.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import shiftgrid
from datetime import datetime, date
from PyFuncemeClimateTools import PlotMaps as pm
from matplotlib.font_manager import FontProperties

from pltsst import plotmap

def import_exp_model(exp, season):

	mod_path = '/home/nice'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2004-2005_yseasmean.nc'.format(mod_path, exp)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	sim      = var1[season,:,:]
	
	return sim
	

def import_obs_data(obs, season):

	vars_dic = {u'cmap': u'precip', u'trmm': u'r'}
	
	pre_path = '/home/nice'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2004-2005_yseasmean.nc'.format(pre_path, obs)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dic[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[season,:,:]

	return pre
	

# Import simulations experiments and observed databases 3D			   
for season in range(0,4):
	
	season_dic = {0: u'DJF 2004/2005', 1: u'MAM 2005', 2: u'JJA 2005', 3: u'SON 2005'}
	print season_dic[season]
	
	exp  = 'exp1'
	exp1 = import_exp_model(exp, season)
	print "exp1", exp1
	print

	exp  = 'exp2'
	exp2 = import_exp_model(exp, season)
	print "exp2", exp2
	print
	
	obs  = 'cmap'
	cmap = import_obs_data(obs, season)
	print "cmap", cmap
	print
	
	obs  = 'trmm'
	trmm = import_obs_data(obs, season)
	print "trmm", trmm
		   

		   
	
