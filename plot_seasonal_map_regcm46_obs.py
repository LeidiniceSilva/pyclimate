# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "05/26/2018"
__description__ = "This script plot precipitation seasonal simulation"


import os
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
# mpl.use('Agg')

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from matplotlib.font_manager import FontProperties

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import interp

from os.path import expanduser


def import_exp_model(exp, season):

	path_in1 = '/home/nice/Documentos/ufrn/papers/regcm_exp/exp_pbl/datas'
	arq1     = '{0}/pre_amz_neb_regcm_{1}_2004-2005_yseasmean.nc'.format(path_in1, exp)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['pr'][:]
	lat      = data1.variables['lat'][:]
	lon      = data1.variables['lon'][:] -360
	sim      = var1[season,:,:]
	
	return lat, lon, sim
	

def import_obs_data(obs, season):

	vars_dic = {u'cmap': u'precip', u'trmm': u'r'}
	
	path_in2 = '/home/nice/Documentos/ufrn/papers/regcm_exp/exp_pbl/datas'
	arq2     = '{0}/pre_amz_neb_{1}_obs_2004-2005_yseasmean.nc'.format(path_in2, obs)
	data2    = netCDF4.Dataset(arq2)
	var2     = data2.variables[vars_dic[obs]][:]
	lat      = data2.variables['lat'][:]
	lon      = data2.variables['lon'][:]
	pre      = var2[season,:,:]

	return lat, lon, pre
	

# Import simulations experiments and observed databases 3D			   
for season in range(0,4):
	
	season_dic = {0: u'DJF-2004/2005', 1: u'MAM-2005', 2: u'JJA-2005', 3: u'SON-2005'}
	print season_dic[season]
	
	exp  = u'exp1'
	lat, lon, exp1 = import_exp_model(exp, season)
	print "exp1", lat, lon, exp1
	print

	exp  = u'exp2'
	exp2 = import_exp_model(exp, season)
	print "exp2", exp2
	print
	
	obs  = u'cmap'
	cmap = import_obs_data(obs, season)
	print "cmap", cmap
	print
	
	obs  = u'trmm'
	trmm = import_obs_data(obs, season)
	print "trmm", trmm

	# Plot precipitation first season amz_neb 
	fig = plt.figure(figsize=(12,10))    

	title = u'Precipitação média sazonal (mm/d) \n{0}'.format(season_dic[season])      
	plt.title(title, size=14)

	s1 = u'Modelo: RegCM4.6'

	plt.text(0, -4, s1, fontsize=9)
	
	levs   = (0, 1, 2, 3, 4, 6, 8, 10, 12, 20, 25, 30)
	colors = ('#2372c9', '#3498ed', '#4ba7ef', '#76bbf3','#93d3f6', '#b0f0f7', '#ffffff', '#fbe78a', '#ff9d37', '#ff5f26', '#ff2e1b', '#ff0219', '#ae000c')
	
	maps = Basemap(projection='cyl', llcrnrlat=-19.75, urcrnrlat=9.75, llcrnrlon=275., urcrnrlon=345.)
	maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)
	maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
	
	cpalunder = colors[0]
	cpalover = colors[-1]
	barcolor = colors[1:-1]
	my_cmap = c.ListedColormap(barcolor)
	my_cmap.set_under(cpalunder)
	my_cmap.set_over(cpalover)
	norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)
	
	plot_maps = plt.contourf(exp1, cmap=my_cmap, norm=norml, levels=levs, extend='both',)
			    
	bar = fig.colorbar(plot_maps, spacing='uniform', ticks=levs, extendfrac='auto', extend='both', drawedges=True)
	bar.set_ticklabels(levs)

	plt.show()
	exit()
	
	path_out = '/home/nice/Documentos/ufrn/papers/regcm_exp/exp_pbl/results/results_new/'
	fig_name = 'pre_season_amz_neb_{0}.png'.format(season_dic[season])

	if not os.path.exists(path_out):
		os.makedirs(path_out)
	
	plt.savefig(os.path.join(path_out, fig_name))
	
	raise SystemExit
	exit()
	

