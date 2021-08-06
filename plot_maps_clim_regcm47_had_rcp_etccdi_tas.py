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

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}

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
	
	
# Import regcm exp and cru databases 
# Historical period	
lat, lon, rcm_txx_hist = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_txn_hist = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tnx_hist = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tnn_hist = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_dtr_hist = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_su_hist = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tr_hist = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tx10p_hist = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tx90p_hist = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tn10p_hist = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, rcm_tn90p_hist = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')

lat, lon, gcm_txx_hist = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_txn_hist = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnx_hist = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnn_hist = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_dtr_hist = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_su_hist = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tr_hist = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx10p_hist = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx90p_hist = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn10p_hist = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn90p_hist = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# RCP26
lat, lon, rcm_txx_rcp26 = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_txn_rcp26 = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tnx_rcp26 = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tnn_rcp26 = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_dtr_rcp26 = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_su_rcp26 = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tr_rcp26 = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tx10p_rcp26 = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tx90p_rcp26 = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tn10p_rcp26 = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
lat, lon, rcm_tn90p_rcp26 = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

lat, lon, gcm_txx_rcp26 = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_txn_rcp26 = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tnx_rcp26 = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tnn_rcp26 = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_dtr_rcp26 = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_su_rcp26 = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tr_rcp26 = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tx10p_rcp26 = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tx90p_rcp26 = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tn10p_rcp26 = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
lat, lon, gcm_tn90p_rcp26 = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

# RCP85
lat, lon, rcm_txx_rcp85 = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_txn_rcp85 = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tnx_rcp85 = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tnn_rcp85 = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_dtr_rcp85 = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_su_rcp85 = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tr_rcp85 = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tx10p_rcp85 = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tx90p_rcp85 = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tn10p_rcp85 = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
lat, lon, rcm_tn90p_rcp85 = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

lat, lon, gcm_txx_rcp85 = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_txn_rcp85 = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tnx_rcp85 = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tnn_rcp85 = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_dtr_rcp85 = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_su_rcp85 = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tr_rcp85 = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tx10p_rcp85 = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tx90p_rcp85 = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tn10p_rcp85 = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
lat, lon, gcm_tn90p_rcp85 = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# Compute difference between periods
# RCP26
diff_rcm_txx_rcp26_hist = np.nanmean(rcm_txx_rcp26, axis=0) - np.nanmean(rcm_txx_hist, axis=0)
diff_rcm_txn_rcp26_hist = np.nanmean(rcm_txn_rcp26, axis=0) - np.nanmean(rcm_txn_hist, axis=0)
diff_rcm_tnx_rcp26_hist = np.nanmean(rcm_tnx_rcp26, axis=0) - np.nanmean(rcm_tnx_hist, axis=0)
diff_rcm_tnn_rcp26_hist = np.nanmean(rcm_tnn_rcp26, axis=0) - np.nanmean(rcm_tnn_hist, axis=0)
diff_rcm_dtr_rcp26_hist = np.nanmean(rcm_dtr_rcp26, axis=0) - np.nanmean(rcm_dtr_hist, axis=0)
diff_rcm_su_rcp26_hist = np.nanmean(rcm_su_rcp26, axis=0) - np.nanmean(rcm_su_hist, axis=0)
diff_rcm_tr_rcp26_hist = np.nanmean(rcm_tr_rcp26, axis=0) - np.nanmean(rcm_tr_hist, axis=0)
diff_rcm_tx10p_rcp26_hist = np.nanmean(rcm_tx10p_rcp26, axis=0) - np.nanmean(rcm_tx10p_hist, axis=0)
diff_rcm_tx90p_rcp26_hist = np.nanmean(rcm_tx90p_rcp26, axis=0) - np.nanmean(rcm_tx90p_hist, axis=0)
diff_rcm_tn10p_rcp26_hist = np.nanmean(rcm_tn10p_rcp26, axis=0) - np.nanmean(rcm_tn10p_hist, axis=0)
diff_rcm_tn90p_rcp26_hist = np.nanmean(rcm_tn90p_rcp26, axis=0) - np.nanmean(rcm_tn90p_hist, axis=0)

diff_gcm_txx_rcp26_hist = gcm_txx_rcp26 - gcm_txx_hist
diff_gcm_txn_rcp26_hist = gcm_txn_rcp26 - gcm_txn_hist
diff_gcm_tnx_rcp26_hist = gcm_tnx_rcp26 - gcm_tnx_hist
diff_gcm_tnn_rcp26_hist = gcm_tnn_rcp26 - gcm_tnn_hist
diff_gcm_dtr_rcp26_hist = gcm_dtr_rcp26 - gcm_dtr_hist
diff_gcm_su_rcp26_hist = gcm_su_rcp26 - gcm_su_hist
diff_gcm_tr_rcp26_hist = gcm_tr_rcp26 - gcm_tr_hist
diff_gcm_tx10p_rcp26_hist = gcm_tx10p_rcp26 - gcm_tx10p_hist
diff_gcm_tx90p_rcp26_hist = gcm_tx90p_rcp26 - gcm_tx90p_hist
diff_gcm_tn10p_rcp26_hist = gcm_tn10p_rcp26 - gcm_tn10p_hist
diff_gcm_tn90p_rcp26_hist = gcm_tn90p_rcp26 - gcm_tn90p_hist

# RCP85
diff_rcm_txx_rcp85_hist = np.nanmean(rcm_txx_rcp85, axis=0) - np.nanmean(rcm_txx_hist, axis=0)
diff_rcm_txn_rcp85_hist = np.nanmean(rcm_txn_rcp85, axis=0)- np.nanmean(rcm_txn_hist, axis=0)
diff_rcm_tnx_rcp85_hist = np.nanmean(rcm_tnx_rcp85, axis=0) - np.nanmean(rcm_tnx_hist, axis=0)
diff_rcm_tnn_rcp85_hist = np.nanmean(rcm_tnn_rcp85, axis=0) - np.nanmean(rcm_tnn_hist, axis=0)
diff_rcm_dtr_rcp85_hist = np.nanmean(rcm_dtr_rcp85, axis=0) - np.nanmean(rcm_dtr_hist, axis=0)
diff_rcm_su_rcp85_hist = np.nanmean(rcm_su_rcp85, axis=0) - np.nanmean(rcm_su_hist, axis=0)
diff_rcm_tr_rcp85_hist = np.nanmean(rcm_tr_rcp85, axis=0) - np.nanmean(rcm_tr_hist, axis=0)
diff_rcm_tx10p_rcp85_hist = np.nanmean(rcm_tx10p_rcp85, axis=0) - np.nanmean(rcm_tx10p_hist, axis=0)
diff_rcm_tx90p_rcp85_hist = np.nanmean(rcm_tx90p_rcp85, axis=0) - np.nanmean(rcm_tx90p_hist, axis=0)
diff_rcm_tn10p_rcp85_hist = np.nanmean(rcm_tn10p_rcp85, axis=0) - np.nanmean(rcm_tn10p_hist, axis=0)
diff_rcm_tn90p_rcp85_hist = np.nanmean(rcm_tn90p_rcp85, axis=0) - np.nanmean(rcm_tn90p_hist, axis=0)

diff_gcm_txx_rcp85_hist = gcm_txx_rcp85 - gcm_txx_hist
diff_gcm_txn_rcp85_hist = gcm_txn_rcp85 - gcm_txn_hist
diff_gcm_tnx_rcp85_hist = gcm_tnx_rcp85 - gcm_tnx_hist
diff_gcm_tnn_rcp85_hist = gcm_tnn_rcp85 - gcm_tnn_hist
diff_gcm_dtr_rcp85_hist = gcm_dtr_rcp85 - gcm_dtr_hist
diff_gcm_su_rcp85_hist = gcm_su_rcp85 - gcm_su_hist
diff_gcm_tr_rcp85_hist = gcm_tr_rcp85 - gcm_tr_hist
diff_gcm_tx10p_rcp85_hist = gcm_tx10p_rcp85 - gcm_tx10p_hist
diff_gcm_tx90p_rcp85_hist = gcm_tx90p_rcp85 - gcm_tx90p_hist
diff_gcm_tn10p_rcp85_hist = gcm_tn10p_rcp85 - gcm_tn10p_hist
diff_gcm_tn90p_rcp85_hist = gcm_tn90p_rcp85 - gcm_tn90p_hist

# Plot maps with the function
fig = plt.figure(figsize=(8, 11))

levs1 = [0.5, 1, 2, 4, 6, 8]
levs2 = [0.1, 0.5, 1, 1.5, 2, 2.5]
levs3 = [10, 30, 50, 70, 90, 110]
levs4 = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]

ax = fig.add_subplot(11, 4, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_txx_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_txx_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_txx_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_txx_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_txn_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_txn_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_txn_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_txn_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_txn_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 10)
map, xx, yy = basemap(lat, lon)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_tnx_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tnx_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tnx_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 13)
map, xx, yy = basemap(lat, lon)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tnn_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_tnn_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tnn_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 16)
map, xx, yy = basemap(lat, lon)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tnn_rcp85_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_dtr_rcp26_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_dtr_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_dtr_rcp26_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_dtr_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_su_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_su_rcp85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_su_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_su_rcp85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 25)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tr_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_tr_rcp85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tr_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 28)
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tr_rcp85_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tx10p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_tx10p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 31)
map, xx, yy = basemap(lat, lon)
plt.title(u'E.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tx10p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 32)
map, xx, yy = basemap(lat, lon)
plt.title(u'F.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tx10p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 33)
map, xx, yy = basemap(lat, lon)
plt.title(u'G.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tx90p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 34)
map, xx, yy = basemap(lat, lon)
plt.title(u'H.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_rcm_tx90p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 35)
map, xx, yy = basemap(lat, lon)
plt.title(u'I.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tx90p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 36)
map, xx, yy = basemap(lat, lon)
plt.title(u'J.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, diff_gcm_tx90p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 37)
map, xx, yy = basemap(lat, lon)
plt.title(u'K.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tn10p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 38)
map, xx, yy = basemap(lat, lon)
plt.title(u'L.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tn10p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 39)
map, xx, yy = basemap(lat, lon)
plt.title(u'M.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_tn10p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 40)
map, xx, yy = basemap(lat, lon)
plt.title(u'N.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_tn10p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 4, 41)
map, xx, yy = basemap(lat, lon)
plt.title(u'K.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tn90p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(11, 4, 42)
map, xx, yy = basemap(lat, lon)
plt.title(u'L.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_rcm_tn90p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 4, 43)
map, xx, yy = basemap(lat, lon)
plt.title(u'M.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_tn90p_rcp26_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(11, 4, 44)
map, xx, yy = basemap(lat, lon)
plt.title(u'N.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, diff_gcm_tn90p_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_etccdi_tas_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.close('all')
plt.cla()
exit()	
	

