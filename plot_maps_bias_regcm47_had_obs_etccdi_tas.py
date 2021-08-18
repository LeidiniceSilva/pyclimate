# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot climatology maps from extremes index"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch
from comp_statist_indices import compute_av


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_txx': u'tmax', 
	u'eca_txn': u'tmax',
	u'eca_tnx': u'tmin', 
	u'eca_tnn': u'tmin',
	u'eca_dtr': u'tmax',
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
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	std   = np.std(var[:][:,:,:], axis=0)
	
	return lat, lon, obs, std
	
	
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
	std   = np.std(var[:][:,:,:], axis=0)

	return lat, lon, rcm, std


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
	std   = np.std(var[:][:,:,:], axis=0)

	return lat, lon, gcm, std


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
	

def ttest(mean_sample1, mean_sample2, std_sample1, std_sample2):

	# Calculate t statistics	
	p1 = mean_sample1 - mean_sample2 
	p2 = (std_sample1 - std_sample2) / np.sqrt(240)

	ttest = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(ttest, df=240)

	return p_value
	
		
# Import extreme indices 
lat, lon, obs_txx, obs_txx_std = import_obs('eca_txx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_txx, rcm_txx_std = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_txx, gcm_txx_std = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_txn, obs_txn_std = import_obs('eca_txn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_txn, rcm_txn_std = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_txn, gcm_txn_std = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tnx, obs_tnx_std = import_obs('eca_tnx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tnx, rcm_tnx_std = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnx, gcm_tnx_std = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tnn, obs_tnn_std = import_obs('eca_tnn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tnn, rcm_tnn_std = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnn, gcm_tnn_std = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_dtr, obs_dtr_std = import_obs('eca_dtr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_dtr, rcm_dtr_std = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_dtr, gcm_dtr_std = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_su, obs_su_std = import_obs('eca_su', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_su, rcm_su_std = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_su, gcm_su_std = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tr, obs_tr_std = import_obs('eca_tr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tr, rcm_tr_std = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tr, gcm_tr_std = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx10p, obs_tx10p_std = import_obs('eca_tx10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx10p, rcm_tx10p_std = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx10p, gcm_tx10p_std = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx90p, obs_tx90p_std = import_obs('eca_tx90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx90p, rcm_tx90p_std = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx90p, gcm_tx90p_std = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn10p, obs_tn10p_std = import_obs('eca_tn10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn10p, rcm_tn10p_std = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn10p, gcm_tn10p_std = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn90p, obs_tn90p_std = import_obs('eca_tn90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn90p, rcm_tn90p_std = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn90p, gcm_tn90p_std = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Compute bias from extreme indices 
bias_rcm_txx = np.nanmean(rcm_txx, axis=0) - obs_txx
bias_rcm_txn = np.nanmean(rcm_txn, axis=0) - obs_txn
bias_rcm_tnx = np.nanmean(rcm_tnx, axis=0) - obs_tnx
bias_rcm_tnn = np.nanmean(rcm_tnn, axis=0) - obs_tnn
bias_rcm_dtr = np.nanmean(rcm_dtr, axis=0) - obs_dtr
bias_rcm_su = np.nanmean(rcm_su, axis=0) - obs_su
bias_rcm_tr = np.nanmean(rcm_tr, axis=0) - obs_tr
bias_rcm_tx10p = np.nanmean(rcm_tx10p, axis=0) - obs_tx10p
bias_rcm_tx90p = np.nanmean(rcm_tx90p, axis=0) - obs_tx90p
bias_rcm_tn10p = np.nanmean(rcm_tn10p, axis=0) - obs_tn10p
bias_rcm_tn90p = np.nanmean(rcm_tn90p, axis=0) - obs_tn90p

bias_gcm_txx = gcm_txx - obs_txx
bias_gcm_txn = gcm_txn - obs_txn
bias_gcm_tnx = gcm_tnx - obs_tnx
bias_gcm_tnn = gcm_tnn - obs_tnn
bias_gcm_dtr = gcm_dtr - obs_dtr
bias_gcm_su = gcm_su - obs_su
bias_gcm_tr = gcm_tr - obs_tr
bias_gcm_tx10p = gcm_tx10p - obs_tx10p
bias_gcm_tx90p = gcm_tx90p - obs_tx90p
bias_gcm_tn10p = gcm_tn10p - obs_tn10p
bias_gcm_tn90p = gcm_tn90p - obs_tn90p

av_txx = compute_av(gcm_txx, np.nanmean(rcm_txx, axis=0), obs_txx)
av_txn = compute_av(gcm_txn, np.nanmean(rcm_txn, axis=0), obs_txn)
av_tnx = compute_av(gcm_tnx, np.nanmean(rcm_tnx, axis=0), obs_tnx)
av_tnn = compute_av(gcm_tnn, np.nanmean(rcm_tnn, axis=0), obs_tnn)
av_dtr = compute_av(gcm_dtr, np.nanmean(rcm_dtr, axis=0), obs_dtr)
av_su = compute_av(gcm_su, np.nanmean(rcm_su, axis=0), obs_su)
av_tr = compute_av(gcm_tr, np.nanmean(rcm_tr, axis=0), obs_tr)
av_tx10p = compute_av(gcm_tx10p, np.nanmean(rcm_tx10p, axis=0), obs_tx10p)
av_tx90p = compute_av(gcm_tx90p, np.nanmean(rcm_tx90p, axis=0), obs_tx90p)
av_tn10p = compute_av(gcm_tn10p, np.nanmean(rcm_tn10p, axis=0), obs_tn10p)
av_tn90p = compute_av(gcm_tn90p, np.nanmean(rcm_tn90p, axis=0), obs_tn90p)

# Compute ttest from models and obs database 
p_value_rcm_txx = ttest(rcm_txx, obs_txx, rcm_txx_std, obs_txx_std)
p_value_rcm_txn = ttest(rcm_txn, obs_txn, rcm_txn_std, obs_txn_std)
p_value_rcm_tnx = ttest(rcm_tnx, obs_tnx, rcm_tnx_std, obs_tnx_std)
p_value_rcm_tnn = ttest(rcm_tnn, obs_tnn, rcm_tnn_std, obs_tnn_std)
p_value_rcm_dtr = ttest(rcm_dtr, obs_dtr, rcm_dtr_std, obs_dtr_std)
p_value_rcm_su = ttest(rcm_su, obs_su, rcm_su_std, obs_su_std)
p_value_rcm_tr = ttest(rcm_tr, obs_tr, rcm_tr_std, obs_tr_std)
p_value_rcm_tx10p = ttest(rcm_tx10p, obs_tx10p, rcm_tx10p_std, obs_tx10p_std)
p_value_rcm_tx90p = ttest(rcm_tx90p, obs_tx90p, rcm_tx90p_std, obs_tx90p_std)
p_value_rcm_tn10p = ttest(rcm_tn10p, obs_tn10p, rcm_tn10p_std, obs_tn10p_std)
p_value_rcm_tn90p = ttest(rcm_tn90p, obs_tn90p, rcm_tn90p_std, obs_tn90p_std)

p_value_gcm_txx = ttest(gcm_txx, obs_txx, gcm_txx_std, obs_txx_std)
p_value_gcm_txn = ttest(gcm_txn, obs_txn, gcm_txn_std, obs_txn_std)
p_value_gcm_tnx = ttest(gcm_tnx, obs_tnx, gcm_tnx_std, obs_tnx_std)
p_value_gcm_tnn = ttest(gcm_tnn, obs_tnn, gcm_tnn_std, obs_tnn_std)
p_value_gcm_dtr = ttest(gcm_dtr, obs_dtr, gcm_dtr_std, obs_dtr_std)
p_value_gcm_su = ttest(gcm_su, obs_su, gcm_su_std, obs_su_std)
p_value_gcm_tr = ttest(gcm_tr, obs_tr, gcm_tr_std, obs_tr_std)
p_value_gcm_tx10p = ttest(gcm_tx10p, obs_tx10p, gcm_tx10p_std, obs_tx10p_std)
p_value_gcm_tx90p = ttest(gcm_tx90p, obs_tx90p, gcm_tx90p_std, obs_tx90p_std)
p_value_gcm_tn10p = ttest(gcm_tn10p, obs_tn10p, gcm_tn10p_std, obs_tn10p_std)
p_value_gcm_tn90p = ttest(gcm_tn90p, obs_tn90p, gcm_tn90p_std, obs_tn90p_std)
# Plot extreme indices 
fig = plt.figure(figsize=(6, 11))

levs1 = [-6, -4, -2, 2, 4, 6]
levs2 = [-60, -40, -20, 20, 40, 60]
levs3 = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
levs11 = [-1, -0.5, -0.1, 0.1, 0.5, 1]

ax = fig.add_subplot(11, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_txx, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_txx = ma.masked_where(p_value_rcm_txx >= 0.05, p_value_rcm_txx) 
map.contourf(xx, yy, p_value_rcm_txx[0,:,:], colors='none', hatches=['....'])
	
ax = fig.add_subplot(11, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_txx, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_txx = ma.masked_where(p_value_gcm_txx >= 0.05, p_value_gcm_txx) 
map.contourf(xx, yy, p_value_gcm_txx, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_txx, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(11, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_txn, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_txn = ma.masked_where(p_value_rcm_txn >= 0.05, p_value_rcm_txn) 
map.contourf(xx, yy, p_value_rcm_txn[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_txn, levels=levs1, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_txn = ma.masked_where(p_value_gcm_txn >= 0.05, p_value_gcm_txn) 
map.contourf(xx, yy, p_value_gcm_txn, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_txn, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tnx, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tnx = ma.masked_where(p_value_rcm_tnx >= 0.05, p_value_rcm_tnx) 
map.contourf(xx, yy, p_value_rcm_tnx[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tnx, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tnx = ma.masked_where(p_value_gcm_tnx >= 0.05, p_value_gcm_tnx) 
map.contourf(xx, yy, p_value_gcm_tnx, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tnx, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tnn, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tnn = ma.masked_where(p_value_rcm_tnn >= 0.05, p_value_rcm_tnn) 
map.contourf(xx, yy, p_value_rcm_tnn[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tnn, levels=levs1, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], BrBG=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tnn = ma.masked_where(p_value_gcm_tnn >= 0.05, p_value_gcm_tnn) 
map.contourf(xx, yy, p_value_gcm_tnn, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tnn, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 13) 
map, xx, yy = basemap(lat, lon)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_dtr, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_dtr = ma.masked_where(p_value_rcm_dtr >= 0.05, p_value_rcm_dtr) 
map.contourf(xx, yy, p_value_rcm_dtr[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_dtr, levels=levs1, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_dtr = ma.masked_where(p_value_gcm_dtr >= 0.05, p_value_gcm_dtr) 
map.contourf(xx, yy, p_value_gcm_dtr, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_dtr, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 16) 
map, xx, yy = basemap(lat, lon)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_su, levels=levs2, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_su = ma.masked_where(p_value_rcm_su >= 0.05, p_value_rcm_su) 
map.contourf(xx, yy, p_value_rcm_su[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_su, levels=levs2, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_su = ma.masked_where(p_value_gcm_su >= 0.05, p_value_gcm_su) 
map.contourf(xx, yy, p_value_gcm_su, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_su, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tr, levels=levs2, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tr = ma.masked_where(p_value_rcm_tr >= 0.05, p_value_rcm_tr) 
map.contourf(xx, yy, p_value_rcm_tr[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tr, levels=levs2, latlon=True, cmap=cm.bwr, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tr = ma.masked_where(p_value_gcm_tr >= 0.05, p_value_gcm_tr) 
map.contourf(xx, yy, p_value_gcm_tr, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tr, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tx10p, levels=levs3, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tx10p = ma.masked_where(p_value_rcm_tx10p >= 0.05, p_value_rcm_tx10p) 
map.contourf(xx, yy, p_value_rcm_tx10p[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tx10p, levels=levs3, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tx10p = ma.masked_where(p_value_gcm_tx10p >= 0.05, p_value_gcm_tx10p) 
map.contourf(xx, yy, p_value_gcm_tx10p, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tx10p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 25)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tx90p, levels=levs3, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tx90p = ma.masked_where(p_value_rcm_tx90p >= 0.05, p_value_rcm_tx90p) 
map.contourf(xx, yy, p_value_rcm_tx90p[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tx90p, levels=levs3, latlon=True, cmap=cm.bwr, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tx90p = ma.masked_where(p_value_gcm_tx90p >= 0.05, p_value_gcm_tx90p) 
map.contourf(xx, yy, p_value_gcm_tx90p, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tx90p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(11, 3, 28) 
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_tn10p, levels=levs3, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tn10p = ma.masked_where(p_value_rcm_tn10p >= 0.05, p_value_rcm_tn10p) 
map.contourf(xx, yy, p_value_rcm_tn10p[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_tn10p, levels=levs3, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tn10p = ma.masked_where(p_value_gcm_tn10p >= 0.05, p_value_gcm_tn10p) 
map.contourf(xx, yy, p_value_gcm_tn10p, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_tn10p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(11, 3, 31) 
map, xx, yy = basemap(lat, lon)
plt.title(u'E.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, bias_rcm_tn90p, levels=levs3, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_tn90p = ma.masked_where(p_value_rcm_tn90p >= 0.05, p_value_rcm_tn90p) 
map.contourf(xx, yy, p_value_rcm_tn90p[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 32)
map, xx, yy = basemap(lat, lon)
plt.title(u'F.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, bias_gcm_tn90p, levels=levs3, latlon=True, cmap=cm.bwr, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_tn90p = ma.masked_where(p_value_gcm_tn90p >= 0.05, p_value_gcm_tn90p) 
map.contourf(xx, yy, p_value_gcm_tn90p, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 33)
map, xx, yy = basemap(lat, lon)
plt.title(u'G.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, av_tn90p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_etccdi_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()	
	
	
	
	
	
