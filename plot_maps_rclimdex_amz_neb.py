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
import seaborn as sns
import shapefile as shp
import matplotlib.cm as cm
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

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
        
        
def plot_map_fill_multiples_ids(comuna, sf, x_lim=None, y_lim=None, color='purple'):

    fig, ax = plt.subplots(figsize=(6,8))
    plt.title(u'Shapefile MATOPIBA')
    plt.xlabel(u'Longitude')
    plt.ylabel(u'Latitude')
        
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        ax.plot(x, y, 'k')
            
    for id in comuna:
        shape_ex = sf.shape(id)
        x_lon = np.zeros((len(shape_ex.points),1))
        y_lat = np.zeros((len(shape_ex.points),1))
        for ip in range(len(shape_ex.points)):
            x_lon[ip] = shape_ex.points[ip][0]
            y_lat[ip] = shape_ex.points[ip][1]
        ax.fill(x_lon,y_lat, color)
             
        x0 = np.mean(x_lon)
        y0 = np.mean(y_lat)
        plt.text(x0, y0, id, fontsize=10)
    
    if (x_lim != None) & (y_lim != None):     
        plt.xlim(x_lim)
        plt.ylim(y_lim)


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
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	map.drawmeridians(np.arange(-85.,-5.,10.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.,15.,5.), size=6, labels=[1,0,0,0], linewidth=0.4)
	map.drawcoastlines(linewidth=1, color='k')
	map.drawcountries(linewidth=1, color='k')
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	return map, xx, yy         
                               

def plot_maps_rclimdex(idx):
		
	fig, ax = plt.subplots()
	levs = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 32]

	plt.title(u'A) Txn (°C)', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, idx, levels=levs, latlon=True, cmap=cm.YlOrRd)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	
		                                  
# Open and read shapefile from matopiba 
shp_path = '/home/nice/Documentos/github_projects/shp/shp_matopiba/shp_matopiba.shp'
sf = shp.Reader(shp_path)
df = read_shapefile(sf)
# Plot shapefile from matopiba 
plot_map_fill_multiples_ids([0, 1, 2, 3], sf, color='purple')

# Import rclimdex database 
lat, lon, idx = import_rclimdex('txnETCCDI')
# Plot rclimdex database
plot_maps_rclimdex(idx)

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/PhD_project/results/rclimdex'
name_out = 'pyplt_maps_txnETCCDI_rclimdex_hadgem_1859-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')





