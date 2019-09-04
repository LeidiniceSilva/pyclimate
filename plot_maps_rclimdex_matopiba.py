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

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap


def import_rclimdex(variable):
	
	path = '/home/nice/Documentos/data_file/cmip_data/cmip5/hadgem2-es_rclimdex/historical'
	arq  = '{0}/{1}_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path, variable)

	data = netCDF4.Dataset(arq)
	var  = data.variables[variable][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	idx  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, idx
	

def read_shapefile(sf):

    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    
    return df
        
       
def basemap(lat, lon):
	
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
	
	map = Basemap(projection='cyl', llcrnrlon=-52., llcrnrlat=-17., urcrnrlon=-40.,urcrnrlat=-1., resolution='c')
	map.drawmeridians(np.arange(-52.,-38.,2.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-17.,2.,2.), size=6, labels=[1,0,0,0], linewidth=0.4)
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-52.,-38.,0.25) 
	lats = np.arange(-17.,2.,0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy         
                               

def plot_maps_rclimdex(idx):
		
	fig = plt.figure()
	ax = fig.add_subplot(111)
	levs = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 32]

	plt.title(u'A) Txn (°C) - MATOPIBA', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, idx, levels=levs, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
	return ax, plt_map
	
		                                  
# Open and read shapefile from matopiba 
shp_path = '/home/nice/Documentos/github_projects/shp/shp_matopiba/shp_matopiba.shp'
sf = shp.Reader(shp_path)
df = read_shapefile(sf)
 
# Import rclimdex database 
lat, lon, idx = import_rclimdex('txnETCCDI')
# Plot rclimdex database
ax, plt_map = plot_maps_rclimdex(idx)

#~ print(sf)
#~ print 
#~ print(df)
#~ exit()

for shape_rec in sf.shapeRecords():
	x = [i[0] for i in shape_rec.shape.points[:]]
	y = [i[1] for i in shape_rec.shape.points[:]]
	ax.plot(x, y, 'k')
					
	vertices = []
	codes = []
	pts = shape_rec.shape.points
	prt = list(shape_rec.shape.parts) + [len(pts)]

	for i in range(len(prt) - 1):
		for j in range(prt[i], prt[i+1]):
			vertices.append((pts[j][0], pts[j][1]))
		codes += [Path.MOVETO]
		codes += [Path.LINETO] * (prt[i+1] - prt[i] -2)
		codes += [Path.CLOSEPOLY]
	clip = Path(vertices, codes)
	clip = PathPatch(clip, transform=ax.transData)
		
for contour in plt_map.collections:
        contour.set_clip_path(clip)

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/PhD_project/results/rclimdex'
name_out = 'pyplt_maps_txnETCCDI_rclimdex_hadgem_matopiba_1859-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()




