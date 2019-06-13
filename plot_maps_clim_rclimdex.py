# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

# Open hadgem2-es model 
path  = '/home/nice/Documentos/data_file/cmip_data/cmip5/hadgem2-es_rclimdex/historical'
arq   = '{0}/tnnETCCDI_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path)

data  = netCDF4.Dataset(arq)
var   = data.variables['tnnETCCDI'][:] 
lat   = data.variables['lat'][:]
lon   = data.variables['lon'][:] 
idx   = var[:,:,:]

aux_lon1 = []
aux_lon2 = []
for l in lon:
	if l <= 180:
		aux_lon1.append(l)
	else:
		aux_lon2.append(l-360)
	
lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
new_lat = lat[::-1]
new_lon = lon[::-1]
    
# Plot maps hadgem2-es model 
fig = plt.figure(figsize=(12,6))

plt.title(u'Yearly Minimum of Daily Minimum Temperature ($^\circ$C)', fontsize=20)
plt.xlabel(u'Longitude', fontsize=20, labelpad=30)
plt.ylabel(u'Latitude', fontsize=20, labelpad=40)

map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15., urcrnrlat=10., resolution='c')
map.drawmeridians(np.arange(-85.,-5.,10.), labels=[0,0,0,1], linewidth=0.2)
map.drawparallels(np.arange(-20.,15.,5.), labels=[1,0,0,0], linewidth=0.2)
map.drawcoastlines(linewidth=1, color='k')
map.drawcountries(linewidth=1, color='k')

xin = np.linspace(map.xmin,map.xmax,20) # nx is the number of x points on the grid
yin = np.linspace(map.ymin,map.ymax,20) # ny in the number of y points on the grid
lons = np.arange(-85.,-5.,0.25) 
lats = np.arange(-20.,15.,-0.25) 
lons, lats = np.meshgrid(new_lon,new_lat)
xx, yy = map(lons,lats)

levs   = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
plot_maps = plt.contourf(xx, yy, idx[10,:,:], levels=levs, cmap=cm.coolwarm, extend='both')

bar = fig.colorbar(plot_maps, spacing='uniform', ticks=levs, drawedges=True)

plt.savefig("test_map_matopiba.png", dpi=300, bbox_inches='tight')
plt.show()
exit()


