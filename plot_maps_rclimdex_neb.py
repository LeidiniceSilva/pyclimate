# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import os
import conda
import netCDF4
import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.cm as cm
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from netCDF4 import Dataset
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap


def import_rclimdex(variable):
	
	path = '/home/nice/Documentos/ufrn/papers/wmrn/data/hadgem2-es'
	arq  = '{0}/{1}_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path, variable)

	data = Dataset(arq)
	var  = data.variables[variable][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	idx  = np.nanmean(var[:][127:146,:,:], axis=0)
	return lat, lon, idx
       
		                                  
# Import rclimdex database 
indice=u'altcwdETCCDI'
lat, lon, idx = import_rclimdex(indice)

# Plot rclimdex database
fig, ax = plt.subplots()

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
                                  
map = Basemap(projection='cyl', llcrnrlon=-50., llcrnrlat=-20., urcrnrlon=-34.,urcrnrlat=0., resolution='c')
map.drawmeridians(np.arange(-48.,-34.,2.), size=6, labels=[0,0,0,1], linewidth=0.4)
map.drawparallels(np.arange(-18.,2.,2.), size=6, labels=[1,0,0,0], linewidth=0.4)
map.drawcoastlines(linewidth=1, color='k')
map.drawcountries(linewidth=1, color='k')
map.drawstates(linewidth=1, color='k')

lons, lats = np.meshgrid(new_lon, new_lat)
xx, yy = map(lons,lats)

plt.title(u'A) CWD (days) - NEB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')

levs = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
plt_map = map.contourf(xx, yy, idx, levels=levs, latlon=True, cmap=cm.YlGn)
cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/papers/wmrn/results'
name_out = 'pyplt_maps_{0}_rclimdex_hadgem_neb_1859-2005.png'.format(indice)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




