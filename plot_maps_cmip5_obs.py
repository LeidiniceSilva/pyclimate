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

from pylab import *
from netCDF4 import Dataset
from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import interp


def import_cmip5_clim(model):
	
	param = 'pr'
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_amz_neb_Amon_{2}_{3}_{4}.nc'.format(path, param,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	mdl_data = np.nanmean(value)
	
	return mdl_data


def import_obs_clim(database):
	
	param = 'tmp'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_amz_neb_{2}_obs_mon_{3}.nc'.format(path,
	param, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	obs_data = np.nanmean(value)
	
	return obs_data
	
	
# Import cmip5 model end obs database 
model = u'HadGEM2-ES'
mdl1  = import_cmip5_clim(model)
print mdl1
exit()s

model = u'ensmean_cmip5'
mdl2  = import_cmip5_clim(model)

database = u'cru_ts4.02'
obs1     = import_obs_clim(database)

# Plot model end obs data maps 
map_type = 'fill' # shade or fill

if map_type == 'fill': 

	deltalat = np.mean(np.diff(lat))/2.
	deltalon = np.mean(np.diff(lon))/2.
	lat = lat - deltalat
	lon = lon - deltalon

fig = plt.figure(figsize=(26, 18))

plt.title(u'Média de precipitação do modelo HadGEM2-ES (mm/dia)', fontsize=12, fontweight='bold')

#~ s1 = u'Peíodo de referência: 1975-2005 /n CMIP5-hist'
#~ plt.text(-85, -24, s1, fontsize=8)

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

path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
fig_name = 'plt_maps_precip_amz_neb_hadgem-es_1975-2005.png'

if not os.path.exists(path_out):
	os.makedirs(path_out)

plt.savefig(os.path.join(path_out, fig_name), dpi=300, bbox_inches='tight')

plt.show()
raise SystemExit
exit()
