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

#Open file
path_in1 = '/home/sullyandro'
arq1     = 'tmax.1dg.cru.1981.2011.nc'
data1    = netCDF4.Dataset(arq1)
var1     = data1.variables['tmx'][:]
lat      = data1.variables['LAT'][:]
lon      = data1.variables['LON'][:] 
obs      = var1[0,:,:]


# PLot maps
map_type = 'shade' # or fill

if map_type == 'shade':

	deltalat = np.mean(np.diff(lat))/2.
	deltalon = np.mean(np.diff(lon))/2.
	lat = lat - deltalat
	lon = lon - deltalon

fig = plt.figure(figsize=(12,10))    

title = u'Temperatura média máxima mensal (ºC)'  
plt.title(title, size=14)

s1 = u'Base de dados: CRU'

plt.text(-90, -66, s1, fontsize=9)

levs   = (10, 15, 20, 25, 30, 35)
colors = ('#ffffff', '#fbe78a', '#ff9d37', '#ff5f26', '#ff2e1b', '#ff0219', '#ae000c')

maps = Basemap(projection='cyl', llcrnrlat=np.min(lat), urcrnrlat=np.max(lat), llcrnrlon=np.min(lon), urcrnrlon=np.max(lon))

maps.drawmeridians(np.arange(maps.lonmin+0.5, maps.lonmax+3+0.5, 10), labels=[0,0,0,1], linewidth=0.1)
maps.drawparallels(np.arange(maps.latmin+0.5, maps.latmax+3+0.5, 5), labels=[1,0,0,0], linewidth=0.1)

x, y = maps(lon, lat)
cpalunder = colors[0]
cpalover = colors[-1]
barcolor = colors[1:-1]
my_cmap = c.ListedColormap(barcolor)
my_cmap.set_under(cpalunder)
my_cmap.set_over(cpalover)
norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)

if map_type=='fill':
	plot_maps = plt.contourf(x, y, obs, cmap=my_cmap, norm=norml, levels=levs, extend='both',)

if map_type=='shade':
	plot_maps = plt.pcolormesh(x, y, obs, cmap=my_cmap, norm=norml)

# 	Drawing the line boundries
maps.drawcoastlines(linewidth=1, color='k')
maps.drawcountries(linewidth=1, color='k')

maps.readshapefile('/home/sullyandro/Downloads/estados_2010', 'estados_2010', drawbounds=True, linewidth=.5, color='k')
# maps.readshapefile('/home/sullyandro/Downloads/estados_2010', 'estados_2010', drawbounds=True, linewidth=.5, color='k')


bar = fig.colorbar(plot_maps, spacing='uniform', ticks=levs, extendfrac='auto', extend='both', drawedges=True)
bar.set_ticklabels(levs)

path_out = '/home/sullyandro/'
fig_name = 'tmx_month_sa.png'

if not os.path.exists(path_out):
	os.makedirs(path_out)

plt.savefig(os.path.join(path_out, fig_name))

raise SystemExit
exit()




	

