# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "12/26/2018"
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
	print "exp1"
	print

	exp  = u'exp2'
	exp2 = import_exp_model(exp, season)
	print "exp2"
	print
	
	obs  = u'cmap'
	cmap = import_obs_data(obs, season)
	print "cmap"
	print
	
	obs  = u'trmm'
	trmm = import_obs_data(obs, season)
	print "trmm"

	# Plot precipitation first season amz_neb 
	map_type = 'fill' # or fill
	
	if map_type == 'fill': # Or fill

		deltalat = np.mean(np.diff(lat))/2.
		deltalon = np.mean(np.diff(lon))/2.
		lat = lat - deltalat
		lon = lon - deltalon
	
	fig = plt.figure(figsize=(12,10))
	
	plt.title(u'Precipitação média sazonal (mm/d) \n RegCM4.6 - {0}'.format(season_dic[season]), fontsize=12, fontweight='bold')
	
	s1 = u'Condições de contorno: ERA15 (1981-2010)'
	plt.text(-85, -24, s1, fontsize=8)
    
	plt.xlabel(u'Longitude', fontsize=12, labelpad=30)
	plt.ylabel(u'Latitude', fontsize=12, labelpad=30)

	levs   = [0, 1, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35]
	colors = ['#FFFFFF', '#D1E6E5', '#A6DDDF', '#9AC6FE', '#5059F9', '#451BB5', '#00FF83', '#00EC0F',
	'#00CD1E', '#F6F76D', '#F9D001', '#FF5600', '#E60000', '#FA394E']

	maps = Basemap(projection='cyl', llcrnrlat=np.min(lat), urcrnrlat=np.max(lat), llcrnrlon=np.min(lon), urcrnrlon=np.max(lon))

	maps.drawmeridians(np.arange(maps.lonmin+0.25, maps.lonmax+3+0.25, 10), labels=[0,0,0,1], linewidth=0.1)
	maps.drawparallels(np.arange(maps.latmin+3, maps.latmax+3, 5), labels=[1,0,0,0], linewidth=0.1)

	x, y = maps(lon, lat)
	cpalunder = colors[0]
	cpalover = colors[-1]
	barcolor = colors[1:-1]
	my_cmap = c.ListedColormap(barcolor)
	my_cmap.set_under(cpalunder)
	my_cmap.set_over(cpalover)
	norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)

	if map_type=='fill':
		plot_maps = plt.contourf(x, y, exp1, cmap=my_cmap, norm=norml, levels=levs, extend='both')

	if map_type=='shade':
		plot_maps = plt.pcolormesh(x, y, exp1, cmap=my_cmap, norm=norml)

	# 	Drawing the line boundries
	maps.drawcoastlines(linewidth=1, color='k')
	maps.drawcountries(linewidth=1, color='k')

	maps.readshapefile('/home/nice/Documentos/shp/states_brazil', 'states_brazil', drawbounds=True, linewidth=.5, color='k')

	cbar_ax = fig.add_axes([0.92, 0.14, 0.04, 0.7])

	bar = fig.colorbar(plot_maps, cax=cbar_ax, spacing='uniform', ticks=levs, extend='both',
	extendfrac='auto', pad=0.05, drawedges=True)
	
	bar.set_ticklabels(levs)

	path_out = '/home/nice/Documentos/ufrn/papers/regcm_exp/exp_pbl/results/'
	fig_name = 'pre_sasonal_amz_neb_exp1.png'

	if not os.path.exists(path_out):
		os.makedirs(path_out)

	plt.savefig(os.path.join(path_out, fig_name), dpi=300, bbox_inches='tight')
	
	plt.show()
	raise SystemExit
	exit()


