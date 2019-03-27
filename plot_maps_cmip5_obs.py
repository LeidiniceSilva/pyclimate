# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot maps from CMIP5 models end OBS basedata"


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


def import_cmip5_mean(model):
	
	param = 'pr' # pr or tas
	area  = 'amz_neb' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:] 
	sim   = var[:,:,:]
			
	return lat, lon, sim


def import_obs_mean(database):
	
	param = 'pre' # pre or tmp
	area  = 'amz_neb' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, 
	database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]	
	obs   = var[:,:,:]

	return lat, lon, obs
	

mdl_list = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM']
          
for mdl in mdl_list:
	print 'CMIP5 Model:', mdl
		
	lat, lon, sim_mean = import_cmip5_mean(mdl)
	
	obs  = u'cru_ts4.02'
	lat, lon, obs_mean = import_obs_mean(obs)
	
	# Plot maps cmip5 model end obs database 
	var_name  = 'pre'
	area_name = 'amz_neb' 
	map_type  = 'fill' # shade or fill

	if map_type == 'fill': 

		deltalat = np.mean(np.diff(lat))/2.
		deltalon = np.mean(np.diff(lon))/2.
		lat = lat - deltalat
		lon = lon - deltalon

	fig = plt.figure(figsize=(12,10))

	fig.suptitle(u'Rainfall Mean - AMZ_NEB \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Reference period: 1850-2005)')

	s1 = u'Condições de contorno: ERA15 (1981-2010)'
	plt.text(-85, -24, s1, fontsize=8)

	plt.xlabel(u'Longitude', fontsize=12, labelpad=30)
	plt.ylabel(u'Latitude', fontsize=12, labelpad=30)

	levs   = [0, 1, 2, 4, 6, 8, 10, 12, 15]
	colors = ['#FFFFFF', '#4734F7', '#1E69A1', '#009B4A', '#58BD2E', '#AFDF10', '#FFF700', '#FF8200', '#FF0300', '#9A0000']

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
		plot_maps = plt.contourf(x, y, sim_mean, cmap=my_cmap, norm=norml, levels=levs, extend='both')

	if map_type=='shade':
		plot_maps = plt.pcolormesh(x, y, sim_mean, cmap=my_cmap, norm=norml)

	# 	Drawing the line boundries
	maps.drawcoastlines(linewidth=1, color='k')
	maps.drawcountries(linewidth=1, color='k')

	maps.readshapefile('/home/nice/Documentos/shp/shp_brasil/br_estados_brasil/states_brazil', 'states_brazil', drawbounds=True, linewidth=.5, color='k')

	cbar_ax = fig.add_axes([0.92, 0.14, 0.04, 0.7])

	bar = fig.colorbar(plot_maps, cax=cbar_ax, spacing='uniform', ticks=levs, extend='both',
	extendfrac='auto', pad=0.05, drawedges=True)

	bar.set_ticklabels(levs)

	path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
	fig_name = 'pyplt_maps_pre_amz_neb_cmip5_cru_1975-2005.png'

	if not os.path.exists(path_out):
		os.makedirs(path_out)

	plt.savefig(os.path.join(path_out, fig_name), dpi=100, bbox_inches='tight')

	plt.show()
	raise SystemExit
	exit()
