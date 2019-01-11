# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot precipitation seasonal cmip5 models"


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


var   = 'pr'
model = 'BCC-CSM1-1'
exp   = 'historical_r1i1p1'
date  = '198001-200512'

mod_path = '/home/nice/Downloads'
arq1     = '{0}/{1}_amz_neb_Amon_{2}_{3}_mon_{4}.nc'.format(mod_path, var, model, exp, date)	
data1    = netCDF4.Dataset(arq1)

var1     = data1.variables['pr'][:]
lat      = data1.variables['lat'][:]
lon      = data1.variables['lon'][:]
mod      = var1[1,:,:]

# Plot precipitation first month amz_neb 
map_type = 'fill' # or fill
		
if map_type == 'fill': # Or fill

	deltalat = np.mean(np.diff(lat))/2.
	deltalon = np.mean(np.diff(lon))/2.
	lat = lat - deltalat
	lon = lon - deltalon

fig = plt.figure(figsize=(12,10))

plt.title(u'Precipitação média mensal (mm/d) \n {0}'.format(model, fontsize=12, fontweight='bold')

s1 = u'CMIP5_historical_r1i1p1'
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
	plot_maps = plt.contourf(x, y, mod, cmap=my_cmap, norm=norml, levels=levs, extend='both')

if map_type=='shade':
	plot_maps = plt.pcolormesh(x, y, mod, cmap=my_cmap, norm=norml)

# Drawing the line boundries
maps.drawcoastlines(linewidth=1, color='k')
maps.drawcountries(linewidth=1, color='k')

maps.readshapefile('/home/nice/Documentos/shp/states_brazil', 'states_brazil', drawbounds=True, linewidth=.5, color='k')

cbar_ax = fig.add_axes([0.92, 0.14, 0.04, 0.7])

bar = fig.colorbar(plot_maps, cax=cbar_ax, spacing='uniform', ticks=levs, extend='both',
extendfrac='auto', pad=0.05, drawedges=True)

bar.set_ticklabels(levs)

path_out = '/home/nice/Downloads'
fig_name = '{0}_amz_neb_Amon_{1}_{2}_mon_{3}.png'.format(var, model, exp, date)

if not os.path.exists(path_out):
	os.makedirs(path_out)

plt.savefig(os.path.join(path_out, fig_name), dpi=300, bbox_inches='tight')

plt.show()
raise SystemExit
exit()


