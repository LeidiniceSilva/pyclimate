# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import os
import conda
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

# Open hadgem2-es model 
path  = '/home/nice/Documentos/data_file/cmip_data/cmip5/hadgem2-es_rclimdex/historical'
arq   = '{0}/prcptotETCCDI_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path)

data  = netCDF4.Dataset(arq)
var   = data.variables['prcptotETCCDI'][:] 
lat   = data.variables['lat'][:]
lon   = data.variables['lon'][:] 
idx   = np.nanmean(var[:][:,:,:], axis=0)

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
fig = plt.figure(figsize=(12,5))

plt.title(u'Annual Total Precipitation in Wet Days (mm)', fontsize=20)
plt.xlabel(u'Longitude', fontsize=20, labelpad=30)
plt.ylabel(u'Latitude', fontsize=20, labelpad=30)

map = Basemap(projection='robin', lon_0=0, lat_0=0, resolution='c')
map.drawmeridians(np.arange(-180.,181.,60.), labels=[0,0,0,1], linewidth=0.4)
map.drawparallels(np.arange(-90.,91.,20.), labels=[1,0,0,0], linewidth=0.4)
map.drawcoastlines(linewidth=1, color='k')
map.drawcountries(linewidth=1, color='k')

xin = np.linspace(map.xmin,map.xmax,20) # nx is the number of x points on the grid
yin = np.linspace(map.ymin,map.ymax,20) # ny in the number of y points on the grid
lons = np.arange(-180,180,0.25) 
lats  = np.arange(90,-90,-0.25) 
lons, lats = np.meshgrid(new_lon,new_lat)
xx, yy = map(lons,lats)

levs   = [10, 100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
plot_maps = plt.contourf(xx, yy, idx, levels=levs, cmap=cm.YlGnBu, extend='both')
bar = fig.colorbar(plot_maps, spacing='uniform', ticks=levs, drawedges=True)

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/PhD_project/results/rclimdex'
name_out = 'pyplt_maps_prcptotETCCDI_rclimdex_hadgem_1859-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
