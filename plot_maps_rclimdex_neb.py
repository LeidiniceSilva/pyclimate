# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import os
import conda
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap


def import_rclimdex_his(variable, date):
	
	path = '/home/nice/Documents/ufrn/papers/wmrn/data/hadgem2-es'
	arq  = '{0}/{1}_yr_HadGEM2-ES_historical_r1i1p1_{2}.nc'.format(path, variable, date)

	data = Dataset(arq)
	var  = data.variables[variable][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	idx  = np.nanmean(var[:][127:146,:,:], axis=0)
	
	return lat, lon, idx

def import_rclimdex_rcp(variable, date):
	
	path = '/home/nice/Documents/ufrn/papers/wmrn/data/hadgem2-es'
	arq  = '{0}/{1}_yr_HadGEM2-ES_rcp85_r1i1p1_{2}.nc'.format(path, variable, date)

	data = Dataset(arq)
	var  = data.variables[variable][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	idx  = np.nanmean(var[:][76:95,:,:], axis=0)
	
	return lat, lon, idx 

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
									  
	map = Basemap(projection='cyl', llcrnrlon=-50., llcrnrlat=-20., urcrnrlon=-34.,urcrnrlat=0., resolution='c')
	map.drawmeridians(np.arange(-48.,-34.,2.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-18.,2.,2.), size=6, labels=[1,0,0,0], linewidth=0.4)
	map.drawcoastlines(linewidth=2, color='k')
	map.drawcountries(linewidth=2, color='k')
	map.drawstates(linewidth=2, color='k')
	map.drawmapboundary()

	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	return map, xx, yy

                                
# Import rclimdex database 
lat, lon, idx_his_cdd = import_rclimdex_his(u'altcddETCCDI', u'1859-2005')
lat, lon, idx_rcp_cdd = import_rclimdex_rcp(u'altcddETCCDI', u'2005-2299')

lat, lon, idx_his_r95p = import_rclimdex_his(u'r95pETCCDI', u'1859-2005')
lat, lon, idx_rcp_r95p = import_rclimdex_rcp(u'r95pETCCDI', u'2005-2299')

lat, lon, idx_his_txx = import_rclimdex_his(u'txxETCCDI', u'1859-2005')
lat, lon, idx_rcp_txx = import_rclimdex_rcp(u'txxETCCDI', u'2005-2299')

# Compute difference into rcp and hist
diff_cdd = idx_rcp_cdd - idx_his_cdd
diff_r95p = idx_rcp_r95p - idx_his_r95p
diff_txx = idx_rcp_txx - idx_his_txx

lons = [-36.70, -44.61, -40.79, -41.86, -39.00, -47.48, -45.93, -37.26, -43.35, -40.46, -43.71, -37.04]
lats = [-9.44, -13.33, -14.88, -11.3, -4.28, -5.53, -9.1, -7.01, -4.86, -9.36, -8.41, -10.95]
trend_cdd = [0.872, 1.190, -0.654, 0.444, 3.16, -0.789, -1.198, 0.718, 0.68, -1.592, 0.596, 0.137]
trend_r95p = [12,14, -2.03, -1.515, 3.4922, 4.2, 9.478, 6.717, -4.784, 6.522, 2.04, 1.021, 14.152]
trend_txx = [0.077, 0.062, 0.087, 0.11, -0.121, 0.066, 0.12, 0.112, 0.169, 0.126, 0.138, -0.068]
    
# Plot trend from xavier dataset
fig = plt.figure(figsize=(12,6))

# Plot firt map 
ax = fig.add_subplot(131)
plt.title(u'A) CDD_Xavier (dias) Período: 1986-2005', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)    

map.plot(-36.70, -9.44, 'r^', markersize=4.36)
map.plot(-44.61, -13.33, 'r^', markersize=5.94)
map.plot(-40.79, -14.88, 'b^', markersize=3.27)
map.plot(-41.86, -11.30, 'r^', markersize=2.22)
map.plot(-39.00, -4.28, 'r^', markersize=14.8)
map.plot(-47.48, -5.53, 'b^', markersize=3.94)
map.plot(-45.93, -9.10, 'b^', markersize=5.99)
map.plot(-37.26, -7.01, 'r^', markersize=4.59)
map.plot(-43.35, -4.86, 'r^', markersize=4.4)
map.plot(-40.46, -9.36, 'b^', markersize=8.9)
map.plot(-43.71, -8.41, 'r^', markersize=3.98)
map.plot(-37.04, -10.95, 'b^', markersize=1.685)
    
# Plot second map 
ax = fig.add_subplot(132)
plt.title(u'B) R95p_Xavier (mm) Período: 1986-2005', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

map.plot(-36.70, -9.44, 'r^', markersize=12)
map.plot(-44.61, -13.33, 'r^', markersize=14)
map.plot(-40.79, -14.88, 'b^', markersize=3.03)
map.plot(-41.86, -11.30, 'b^', markersize=2.515)
map.plot(-39.00, -4.28, 'r^', markersize=3.4822)
map.plot(-47.48, -5.53, 'r^', markersize=4.2)
map.plot(-45.93, -9.10, 'r^', markersize=9.478)
map.plot(-37.26, -7.01, 'r^', markersize=6.717)
map.plot(-43.35, -4.86, 'b^', markersize=4.784)
map.plot(-40.46, -9.36, 'r^', markersize=6.522)
map.plot(-43.71, -8.41, 'r^', markersize=3.04)
map.plot(-37.04, -10.95, 'r^', markersize=2.021)

# Plot thirth map 
ax = fig.add_subplot(133)
plt.title(u'B) Txx_Xavier (ºC) Período: 1986-2005', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

map.plot(-36.70, -9.44, 'r^', markersize=4.85)
map.plot(-44.61, -13.33, 'r^', markersize=4.10)
map.plot(-40.79, -14.88, 'r^', markersize=5.35)
map.plot(-41.86, -11.30, 'r^', markersize=6.50)
map.plot(-39.00, -4.28, 'b^', markersize=7.05)
map.plot(-47.48, -5.53, 'r^', markersize=4.30)
map.plot(-45.93, -9.10, 'r^', markersize=7.00)
map.plot(-37.26, -7.01, 'r^', markersize=6.56)
map.plot(-43.35, -4.86, 'r^', markersize=9.4)
map.plot(-40.46, -9.36, 'r^', markersize=7.3)
map.plot(-43.71, -8.41, 'r^', markersize=6.9)
map.plot(-37.04, -10.95, 'b^', markersize=4.40)

#~ def get_marker_color(trend):
    #~ if trend > 0.0:
        #~ return ('r>')
    #~ elif trend < 0.0:
        #~ return ('b>')
    #~ else:
        #~ return ('ro')

#~ for n in range(len(lats)):
    #~ x,y = map(lons, lats)
    #~ marker_string = get_marker_color(trend_r95p[n])
    #~ map.plot(x[n], y[n], marker_string, markersize=trend_r95p[n]*1)
    #~ print(x[n], y[n], marker_string, trend_r95p[n]*1)
#~ exit()

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/wmrn/results'
name_out = 'pyplt_maps_trend_rclimdex_xavier.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

# Plot difference from rclimdex Far future (rcp85) and historical
fig = plt.figure(figsize=(12,6))

levs1 = [-14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14]
ax = fig.add_subplot(131)
plt.title(u'A) CDD_HadGEM2-ES (dias) RCP85 - Hist', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')

map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_cdd, levels=levs1, latlon=True, cmap=cm.bwr)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

lons = [-36.70, -44.61, -40.79, -41.86, -39.00, -47.48, -45.93, -37.26, -43.35, -40.46, -43.71, -37.04]
lats = [-9.44, -13.33, -14.88, -11.3, -4.28, -5.53, -9.1, -7.01, -4.86, -9.36, -8.41, -10.95]
x,y = map(lons, lats)
map.plot(x, y, 'k^', markersize=5)

cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

# Plot second map 
levs2 = [-500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500]
ax = fig.add_subplot(132)
plt.title(u'B) R95p_HadGEM2-ES (mm) RCP85 - Hist', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=20, fontweight='bold')

map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_r95p, levels=levs2, latlon=True, cmap=cm.RdBu)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

lons = [-36.70, -44.61, -40.79, -41.86, -39.00, -47.48, -45.93, -37.26, -43.35, -40.46, -43.71, -37.04]
lats = [-9.44, -13.33, -14.88, -11.3, -4.28, -5.53, -9.1, -7.01, -4.86, -9.36, -8.41, -10.95]
x,y = map(lons, lats)
map.plot(x, y, 'k^', markersize=5)

cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

# Plot three map 
levs3 = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
ax = fig.add_subplot(133)
plt.title(u'B) Txx_HadGEM2-ES (ºC) RCP85 - Hist', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=20, fontweight='bold')

map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_txx, levels=levs3, latlon=True, cmap=cm.Reds)
ax.text(-35, -18,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

lons = [-36.70, -44.61, -40.79, -41.86, -39.00, -47.48, -45.93, -37.26, -43.35, -40.46, -43.71, -37.04]
lats = [-9.44, -13.33, -14.88, -11.3, -4.28, -5.53, -9.1, -7.01, -4.86, -9.36, -8.41, -10.95]
x,y = map(lons, lats)
map.plot(x, y, 'k^', markersize=5)

cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/wmrn/results'
name_out = 'pyplt_maps_diff_rclimdex_hadgem_neb_his_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()






