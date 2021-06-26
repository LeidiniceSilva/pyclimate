# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot RCP2.6 and RCP8.5 maps from extremes index"

import os
import conda
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")
import matplotlib.cm as cm

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, gcm


def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
# Import regcm and gcm databases 	
# Historical period
lat, lon, rcm_prcptot_hist = import_rcm('eca_prcptot', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_r95p_hist = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_r99p_hist = import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_rx1day_hist = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_rx5day_hist = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_sdii_hist = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_cdd_hist = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_cwd_hist = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_r10mm_hist = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_r20mm_hist = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')

lat, lon, gcm_prcptot_hist = import_gcm('eca_prcptot', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r95p_hist = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r99p_hist = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx1day_hist = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx5day_hist = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_sdii_hist = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cdd_hist = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cwd_hist = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r10mm_hist = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r20mm_hist = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# RCP2.6
lat, lon, rcm_prcptot_rcp26 = import_rcm('eca_prcptot', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_r95p_rcp26 = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_r99p_rcp26= import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_rx1day_rcp26 = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_rx5day_rcp26 = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_sdii_rcp26 = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_cdd_rcp26 = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_cwd_rcp26 = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_r10mm_rcp26 = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_r20mm_rcp26 = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

lat, lon, gcm_prcptot_rcp26 = import_gcm('eca_prcptot', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_r95p_rcp26 = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_r99p_rcp26 = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_rx1day_rcp26 = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_rx5day_rcp26 = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_sdii_rcp26 = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_cdd_rcp26 = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_cwd_rcp26 = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_r10mm_rcp26 = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_r20mm_rcp26 = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

# RCP8.5
lat, lon, rcm_prcptot_rcp85 = import_rcm('eca_prcptot', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_r95p_rcp85 = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_r99p_rcp85 = import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_rx1day_rcp85 = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_rx5day_rcp85 = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_sdii_rcp85 = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_cdd_rcp85 = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_cwd_rcp85 = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_r10mm_rcp85 = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_r20mm_rcp85 = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

lat, lon, gcm_prcptot_rcp85 = import_gcm('eca_prcptot', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_r95p_rcp85 = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_r99p_rcp85 = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_rx1day_rcp85 = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_rx5day_rcp85 = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_sdii_rcp85 = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_cdd_rcp85 = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_cwd_rcp85 = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_r10mm_rcp85 = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_r20mm_rcp85 = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# Compute difference between periods
# RCP26
diff_rcm_prcptot_rcp26_hist = rcm_prcptot_rcp26 - rcm_prcptot_hist
diff_rcm_r95p_rcp26_hist = rcm_r95p_rcp26 - rcm_r95p_hist
diff_rcm_r99p_rcp26_hist = rcm_r99p_rcp26 - rcm_r99p_hist
diff_rcm_rx1day_rcp26_hist = rcm_rx1day_rcp26 - rcm_rx1day_hist
diff_rcm_rx5day_rcp26_hist = rcm_rx5day_rcp26 - rcm_rx5day_hist
diff_rcm_sdii_rcp26_hist = rcm_sdii_rcp26 - rcm_sdii_hist
diff_rcm_cdd_rcp26_hist = rcm_cdd_rcp26 - rcm_cdd_hist
diff_rcm_cwd_rcp26_hist = rcm_cwd_rcp26 - rcm_cwd_hist
diff_rcm_r10mm_rcp26_hist = rcm_r10mm_rcp26 - rcm_r10mm_hist
diff_rcm_r20mm_rcp26_hist = rcm_r20mm_rcp26 - rcm_r20mm_hist

diff_gcm_prcptot_rcp26_hist = gcm_prcptot_rcp26 - rcm_prcptot_hist
diff_gcm_r95p_rcp26_hist = gcm_r95p_rcp26 - rcm_r95p_hist
diff_gcm_r99p_rcp26_hist = gcm_r99p_rcp26 - rcm_r99p_hist
diff_gcm_rx1day_rcp26_hist = gcm_rx1day_rcp26 - rcm_rx1day_hist
diff_gcm_rx5day_rcp26_hist = gcm_rx5day_rcp26 - rcm_rx5day_hist
diff_gcm_sdii_rcp26_hist = gcm_sdii_rcp26 - rcm_sdii_hist
diff_gcm_cdd_rcp26_hist = gcm_cdd_rcp26 - rcm_cdd_hist
diff_gcm_cwd_rcp26_hist = gcm_cwd_rcp26 - rcm_cwd_hist
diff_gcm_r10mm_rcp26_hist = gcm_r10mm_rcp26 - rcm_r10mm_hist
diff_gcm_r20mm_rcp26_hist = gcm_r20mm_rcp26 - rcm_r20mm_hist

# RCP85
diff_rcm_prcptot_rcp85_hist = rcm_prcptot_rcp85 - rcm_prcptot_hist
diff_rcm_r95p_rcp85_hist = rcm_r95p_rcp85 - rcm_r95p_hist
diff_rcm_r99p_rcp85_hist = rcm_r99p_rcp85 - rcm_r99p_hist
diff_rcm_rx1day_rcp85_hist = rcm_rx1day_rcp85 - rcm_rx1day_hist
diff_rcm_rx5day_rcp85_hist = rcm_rx5day_rcp85 - rcm_rx5day_hist
diff_rcm_sdii_rcp85_hist = rcm_sdii_rcp85 - rcm_sdii_hist
diff_rcm_cdd_rcp85_hist = rcm_cdd_rcp85 - rcm_cdd_hist
diff_rcm_cwd_rcp85_hist = rcm_cwd_rcp85 - rcm_cwd_hist
diff_rcm_r10mm_rcp85_hist = rcm_r10mm_rcp85 - rcm_r10mm_hist
diff_rcm_r20mm_rcp85_hist = rcm_r20mm_rcp85 - rcm_r20mm_hist

diff_gcm_prcptot_rcp85_hist = gcm_prcptot_rcp85 - rcm_prcptot_hist
diff_gcm_r95p_rcp85_hist = gcm_r95p_rcp85- rcm_r95p_hist
diff_gcm_r99p_rcp85_hist = gcm_r99p_rcp85 - rcm_r99p_hist
diff_gcm_rx1day_rcp85_hist = gcm_rx1day_rcp85 - rcm_rx1day_hist
diff_gcm_rx5day_rcp85_hist = gcm_rx5day_rcp85 - rcm_rx5day_hist
diff_gcm_sdii_rcp85_hist = gcm_sdii_rcp85 - rcm_sdii_hist
diff_gcm_cdd_rcp85_hist = gcm_cdd_rcp85 - rcm_cdd_hist
diff_gcm_cwd_rcp85_hist = gcm_cwd_rcp85 - rcm_cwd_hist
diff_gcm_r10mm_rcp85_hist = gcm_r10mm_rcp85 - rcm_r10mm_hist
diff_gcm_r20mm_rcp85_hist = gcm_r20mm_rcp85 - rcm_r20mm_hist

# Plot maps with the function
fig = plt.figure(figsize=(8, 10))
levs1 = [-150, -100, -50, 50, 100, 150]
levs2 = [-11, -7, -3, 3, 7, 11]
levs3 = [-11, -7, -3, 3, 7, 11]
levs4 = [-11, -7, -3, 3, 7, 11]
levs5 = [-75, -50, -25, 25, 50, 75]
levs6 = [-5, -3, -1., 1, 3, 5]
levs7 = [-15, -10, -5, 5, 10, 15]
levs8 = [-10, -7, -4, 4, 7, 10]
levs9 = [-15, -10, -5, 5, 10, 15]
levs10 = [-10, -7, -4, 4, 7, 10]

levs11 = [-150, -100, -50, 50, 100, 150]
levs21 = [-90, -60, -30, 30, 60, 90]
levs31 = [-90, -60, -30, 30, 60, 90]
levs41 = [-90, -60, -30, 30, 60, 90]
levs51 = [-90, -60, -30, 30, 60, 90]
levs61 = [-9, -6, -3., 3, 6, 9]
levs71 = [-20, -15, -10, 10, 15, 20]
levs81 = [-15, -10, -5, 5, 10, 15]
levs91 = [-20, -15, -10, 10, 15, 20]
levs101 = [-15, -10, -5, 5, 10, 15]

ax = fig.add_subplot(10, 4, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_prcptot_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_prcptot_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs0, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_prcptot_rcp26_hist, levels=levs11, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_prcptot_rcp85_hist, levels=levs11, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs0, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_r95p_rcp26_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_r95p_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r95p_rcp26_hist, levels=levs21, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r95p_rcp85_hist, levels=levs21, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_r99p_rcp26_hist, levels=levs3, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 10)
map, xx, yy = basemap(lat, lon)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_r99p_rcp85_hist, levels=levs3, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r99p_rcp26_hist, levels=levs31, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r99p_rcp85_hist, levels=levs31, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 13)
map, xx, yy = basemap(lat, lon)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_rx1day_rcp26_hist, levels=levs4, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_rx1day_rcp85_hist, levels=levs4, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_rx1day_rcp26_hist, levels=levs41, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 16)
map, xx, yy = basemap(lat, lon)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_rx1day_rcp85_hist, levels=levs41, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_rx5day_rcp26_hist, levels=levs5, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_rx5day_rcp85_hist, levels=levs5, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_rx5day_rcp26_hist, levels=levs51, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_rx5day_rcp85_hist, levels=levs51, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_sdii_rcp26_hist, levels=levs6, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_sdii_rcp85_hist, levels=levs6, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_sdii_rcp26_hist, levels=levs61, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_sdii_rcp85_hist, levels=levs61, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 25)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_cdd_rcp26_hist, levels=levs7, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_cdd_rcp85_hist, levels=levs7, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_cdd_rcp26_hist, levels=levs71, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 28)
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_cdd_rcp85_hist, levels=levs71, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_cwd_rcp26_hist, levels=levs8, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_cwd_rcp85_hist, levels=levs8, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 31)
map, xx, yy = basemap(lat, lon)
plt.title(u'E.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_cwd_rcp26_hist, levels=levs81, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 32)
map, xx, yy = basemap(lat, lon)
plt.title(u'F.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_cwd_rcp85_hist, levels=levs81, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 33)
map, xx, yy = basemap(lat, lon)
plt.title(u'G.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_r10mm_rcp26_hist, levels=levs9, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 34)
map, xx, yy = basemap(lat, lon)
plt.title(u'H.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_r10mm_rcp85_hist, levels=levs9, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 35)
map, xx, yy = basemap(lat, lon)
plt.title(u'I.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r10mm_rcp26_hist, levels=levs91, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 36)
map, xx, yy = basemap(lat, lon)
plt.title(u'J.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_r10mm_rcp85_hist, levels=levs91, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(10, 4, 37)
map, xx, yy = basemap(lat, lon)
plt.title(u'K.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_r20mm_rcp26_hist, levels=levs10, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(10, 4, 38)
map, xx, yy = basemap(lat, lon)
plt.title(u'L.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_r20mm_rcp85_hist, levels=levs10, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 4, 39)
map, xx, yy = basemap(lat, lon)
plt.title(u'M.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_r20mm_rcp26_hist, levels=levs101, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(10, 4, 40)
map, xx, yy = basemap(lat, lon)
plt.title(u'N.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_r20mm_rcp85_hist, levels=levs101, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_etccdi_pre_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.close('all')
plt.cla()
exit()	
	
	
	

#~ ax = fig.add_subplot(10, 4, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_prcptot_rcp26_hist, levels=levs0, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_prcptot_rcp85_hist, levels=levs0, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs0, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_prcptot_rcp85_hist, levels=levs0, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_prcptot_rcp85_hist, levels=levs0, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs0, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_r95p_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 6)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r95p_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 7)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_r95p_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 8)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r95p_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 9)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_r99p_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 10)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r99p_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 11)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_r99p_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 12)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r99p_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 13)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_rx1day_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 14)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_rx1day_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 15)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_rx1day_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 16)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_rx1day_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 17)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_rx5day_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 18)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_rx5day_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 19)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_rx5day_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 20)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_rx5day_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 21)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_sdii_rcp26_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 22)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_sdii_rcp26_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 23)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_sdii_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 24)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_sdii_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 25)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_cdd_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 26)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_cdd_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 27)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_cdd_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 28)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_cdd_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 29)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_cwd_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 30)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_cwd_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 31)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'E.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_cwd_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 32)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'F.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_cwd_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 33)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'G.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_r10mm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 34)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'H.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r10mm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 35)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'I.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_rcm_r10mm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 36)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'J.1)', loc='left', fontsize=8, fontweight='bold')
#~ map.contourf(xx, yy, diff_gcm_r10mm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)

#~ ax = fig.add_subplot(10, 4, 37)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'K.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
#~ plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_r20mm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
#~ ax = fig.add_subplot(10, 4, 38)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'L.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_gcm_r20mm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6) 

#~ ax = fig.add_subplot(10, 4, 39)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'M.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_rcm_r20mm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(10, 4, 40)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'N.1)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
#~ map.contourf(xx, yy, diff_gcm_r20mm_rcp85_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
	
	
