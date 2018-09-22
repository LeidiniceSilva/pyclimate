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
from mpl_toolkits.basemap import shiftgrid

from datetime import datetime, date
from PyFuncemeClimateTools import PlotMaps as pm
from hidropy.utils.hidropy_utils import create_path
from vol3.nice.codes.PyClimateTools.comp_climate_stats import compute_anomaly
from matplotlib.font_manager import FontProperties


def import_exp_model(exp, season):

	mod_path = '/home/nice'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2004-2005_yseasmean.nc'.format(mod_path, exp)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:]
	sim      = var1[season][:,:,:]
	
	return sim
	

def import_obs_data(obs, season):

	vars_dic = {u'cmap': u'precip', u'trmm': u'r'}
	
	pre_path = '/home/nice'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2004-2005_yseasmean.nc'.format(pre_path, obs)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dic[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre     = var2[season][:,:,:]

	return pre
	

# Import simulations experiments and observed databases 3D
season_dic = {u'0': u'DJF 2004/2005', u'1': u'MAM 2005', 
			  u'2': u'JJA 2005', u'3': u'SON 2005' }
			   
for in season in range(0,4):
	print season_list[season]
	
	exp  = 'exp1'
	exp1 = import_exp_model(exp, season)
	print "exp1"
	print

	exp  = 'exp2'
	exp2 = import_exp_model(exp, season)
	print "exp2"
	print
	
	obs  = 'cmap'
	cmap = import_obs_data(obs, season)
	print "cmap", cmap
	print
	
	obs  = 'trmm'
	trmm = import_obs_data(obs, season)
	print "trmm", trmm
	
	# RegCM4.6_EXP1 plotmaps
	title1 = u'Precipitação Sazonal - RegCM4.6_EXP1 - {0}'.format(season_list[season])
	figou1 = 'precip_season_map_amz_neb_regcm4.6_exp1_{0}.png'.format(season_list[season])
	
	pm.plotmap(exp1[season, :, :], lats, lons, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
		   loneastpoint=x2, ocean_mask=1, fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1,
		   barinf='max', barloc='right')
		   
	# RegCM4.6_EXP2 plotmaps
	title1 = u'Precipitação Sazonal - {0}'.format(season_list[season])
	figou1 = 'precip_season_map_amz_neb_{0}.png'.format(season_list[season])
	
	pm.plotmap(exp2[season, :, :], lats, lons, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
		   loneastpoint=x2, ocean_mask=1, fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1,
		   barinf='max', barloc='right')

	# CMAP plotmaps
	title1 = u'Precipitação Sazonal - {0}'.format(season_list[season])
	figou1 = 'precip_season_map_amz_neb_{0}.png'.format(season_list[season])
	
	pm.plotmap(cmap[season, :, :], lats, lons, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
		   loneastpoint=x2, ocean_mask=1, fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1,
		   barinf='max', barloc='right')

	# TRMM plotmaps   
	title1 = u'Precipitação Sazonal - {0}'.format(season_list[season])
	figou1 = 'precip_season_map_amz_neb_{0}.png'.format(season_list[season])
	
	pm.plotmap(trmm[season, :, :], lats, lons, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
		   loneastpoint=x2, ocean_mask=1, fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1,
		   barinf='max', barloc='right')
		   
	
