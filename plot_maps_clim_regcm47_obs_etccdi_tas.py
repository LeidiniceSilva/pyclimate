# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot climatology maps from extremes indices"

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
	sample  = np.nanmean(var[:][0:2,:,:], axis=0)
		
	return lat, lon, obs, std, sample
	
	
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
	sample  = np.nanmean(var[:][0:2,:,:], axis=0)
	
	return lat, lon, rcm, std, sample


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
	sample  = np.nanmean(var[:][0:2,:,:], axis=0)
	
	return lat, lon, gcm, std, sample


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
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=0.8)

	return map, xx, yy
	

def ttest(mean, std, sample):

	# Calculate t statistics
	p1 = mean - sample
	p2= std / np.sqrt(240)
	p3 = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(p3, df=240)
	
	return p_value
	
		
# Import extreme indices 
lat, lon, obs_txx, obs_txx_std, obs_txx_sample = import_obs('eca_txx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_txx, rcm_txx_std, rcm_txx_sample = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_txx, gcm_txx_std, gcm_txx_sample = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_txn, obs_txn_std, obs_txn_sample = import_obs('eca_txn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_txn, rcm_txn_std, rcm_txn_sample = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_txn, gcm_txn_std, gcm_txn_sample = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tnx, obs_tnx_std, obs_tnx_sample = import_obs('eca_tnx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tnx, rcm_tnx_std, rcm_tnx_sample = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnx, gcm_tnx_std, gcm_tnx_sample = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tnn, obs_tnn_std, obs_tnn_sample = import_obs('eca_tnn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tnn, rcm_tnn_std, rcm_tnn_sample = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tnn, gcm_tnn_std, gcm_tnn_sample = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_dtr, obs_dtr_std, obs_dtr_sample = import_obs('eca_dtr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_dtr, rcm_dtr_std, rcm_dtr_sample = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_dtr, gcm_dtr_std, gcm_dtr_sample = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_su, obs_su_std, obs_su_sample = import_obs('eca_su', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_su, rcm_su_std, rcm_su_sample = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_su, gcm_su_std, gcm_su_sample = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tr, obs_tr_std, obs_tr_sample = import_obs('eca_tr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tr, rcm_tr_std, rcm_tr_sample = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tr, gcm_tr_std, gcm_tr_sample = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx10p, obs_tx10p_std, obs_tx10p_sample = import_obs('eca_tx10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx10p, rcm_tx10p_std, rcm_tx10p_sample = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx10p, gcm_tx10p_std, gcm_tx10p_sample = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx90p, obs_tx90p_std, obs_tx90p_sample = import_obs('eca_tx90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx90p, rcm_tx90p_std, rcm_tx90p_sample = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx90p, gcm_tx90p_std, gcm_tx90p_sample = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn10p, obs_tn10p_std, obs_tn10p_sample = import_obs('eca_tn10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn10p, rcm_tn10p_std, rcm_tn10p_sample = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn10p, gcm_tn10p_std, gcm_tn10p_sample = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn90p, obs_tn90p_std, obs_tn90p_sample = import_obs('eca_tn90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn90p, rcm_tn90p_std, rcm_tn90p_sample = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn90p, gcm_tn90p_std, gcm_tn90p_sample = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Compute ttest from models and obs database 
p_value_txx_obs = ttest(obs_txx, obs_txx_std, obs_txx_sample)
p_value_txn_obs = ttest(obs_txn, obs_txn_std, obs_txn_sample)
p_value_tnx_obs = ttest(obs_tnx, obs_tnx_std, obs_tnx_sample)
p_value_tnn_obs = ttest(obs_tnn, obs_tnn_std, obs_tnn_sample)
p_value_dtr_obs = ttest(obs_dtr, obs_dtr_std, obs_dtr_sample)
p_value_su_obs = ttest(obs_su, obs_su_std, obs_su_sample)
p_value_tr_obs = ttest(obs_tr, obs_tr_std, obs_tr_sample)
p_value_tx10p_obs = ttest(obs_tx10p, obs_tx10p_std, obs_tx10p_sample)
p_value_tx90p_obs = ttest(obs_tx90p, obs_tx90p_std, obs_tx90p_sample)
p_value_tn10p_obs = ttest(obs_tn10p, obs_tn10p_std, obs_tn10p_sample)
p_value_tn90p_obs = ttest(obs_tn90p, obs_tn90p_std, obs_tn90p_sample)

p_value_txx_rcm = ttest(rcm_txx, rcm_txx_std, rcm_txx_sample)
p_value_txn_rcm = ttest(rcm_txn, rcm_txn_std, rcm_txn_sample)
p_value_tnx_rcm = ttest(rcm_tnx, rcm_tnx_std, rcm_tnx_sample)
p_value_tnn_rcm = ttest(rcm_tnn, rcm_tnn_std, rcm_tnn_sample)
p_value_dtr_rcm = ttest(rcm_dtr, rcm_dtr_std, rcm_dtr_sample)
p_value_su_rcm = ttest(rcm_su, rcm_su_std, rcm_su_sample)
p_value_tr_rcm = ttest(rcm_tr, rcm_tr_std, rcm_tr_sample)
p_value_tx10p_rcm = ttest(rcm_tx10p, rcm_tx10p_std, rcm_tx10p_sample)
p_value_tx90p_rcm = ttest(rcm_tx90p, rcm_tx90p_std, rcm_tx90p_sample)
p_value_tn10p_rcm = ttest(rcm_tn10p, rcm_tn10p_std, rcm_tn10p_sample)
p_value_tn90p_rcm = ttest(rcm_tn90p, rcm_tn90p_std, rcm_tn90p_sample)

p_value_txx_gcm = ttest(gcm_txx, gcm_txx_std, gcm_txx_sample)
p_value_txn_gcm = ttest(gcm_txn, gcm_txn_std, gcm_txn_sample)
p_value_tnx_gcm = ttest(gcm_tnx, gcm_tnx_std, gcm_tnx_sample)
p_value_tnn_gcm = ttest(gcm_tnn, gcm_tnn_std, gcm_tnn_sample)
p_value_dtr_gcm = ttest(gcm_dtr, gcm_dtr_std, gcm_dtr_sample)
p_value_su_gcm = ttest(gcm_su, gcm_su_std, gcm_su_sample)
p_value_tr_gcm = ttest(gcm_tr, gcm_tr_std, gcm_tr_sample)
p_value_tx10p_gcm = ttest(gcm_tx10p, gcm_tx10p_std, gcm_tx10p_sample)
p_value_tx90p_gcm = ttest(gcm_tx90p, gcm_tx90p_std, gcm_tx90p_sample)
p_value_tn10p_gcm = ttest(gcm_tn10p, gcm_tn10p_std, gcm_tn10p_sample)
p_value_tn90p_gcm = ttest(gcm_tn90p, gcm_tn90p_std, gcm_tn90p_sample)

# Plot extreme indices 
fig = plt.figure(figsize=(6, 11))

levs1 = [28, 31, 34, 37, 40, 43]
levs2 = [15, 18, 21, 24, 27, 30]
levs3 = [17, 20, 23, 26, 29, 32]
levs4 = [9, 12, 15, 18, 21, 24]
levs5 = [1, 4, 7, 10, 13, 16, 19]
levs6 = [260, 280, 300, 320, 340, 366]
levs7 = [160, 200, 240, 280, 320, 366]
levs8 = [0, 4, 8, 12, 16, 20]

ax = fig.add_subplot(11, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_txx, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_txx_obs = ma.masked_where(p_value_txx_obs >= 0.05, p_value_txx_obs) 
map.contourf(xx, yy, p_value_txx_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_txx[0,:,:], levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_txx_rcm = ma.masked_where(p_value_txx_rcm >= 0.05, p_value_txx_rcm) 
map.contourf(xx, yy, p_value_txx_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_txx, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black') 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.set_label('TXx (°C)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_txx_gcm = ma.masked_where(p_value_txx_gcm >= 0.05, p_value_txx_gcm) 
map.contourf(xx, yy, p_value_txx_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_txn, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_txn_obs = ma.masked_where(p_value_txn_obs >= 0.05, p_value_txn_obs) 
map.contourf(xx, yy, p_value_txn_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_txn[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_txn_rcm = ma.masked_where(p_value_txn_rcm >= 0.05, p_value_txn_rcm) 
map.contourf(xx, yy, p_value_txn_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_txn, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.set_label('TXn (°C)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_txn_gcm = ma.masked_where(p_value_txn_gcm >= 0.05, p_value_txn_gcm) 
map.contourf(xx, yy, p_value_txn_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tnx, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tnx_obs = ma.masked_where(p_value_tnx_obs >= 0.05, p_value_tnx_obs) 
map.contourf(xx, yy, p_value_tnx_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tnx[0,:,:], levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tnx_rcm = ma.masked_where(p_value_tnx_rcm >= 0.05, p_value_tnx_rcm) 
map.contourf(xx, yy, p_value_tnx_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tnx, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.set_label('TNx (°C)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tnx_gcm = ma.masked_where(p_value_tnx_gcm >= 0.05, p_value_tnx_gcm) 
map.contourf(xx, yy, p_value_tnx_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tnn, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tnn_obs = ma.masked_where(p_value_tnn_obs >= 0.05, p_value_tnn_obs) 
map.contourf(xx, yy, p_value_tnn_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tnn[0,:,:], levels=levs4, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tnn_rcm = ma.masked_where(p_value_tnn_rcm >= 0.05, p_value_tnn_rcm) 
map.contourf(xx, yy, p_value_tnn_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tnn, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.set_label('TNn (°C)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tnn_gcm = ma.masked_where(p_value_tnn_gcm >= 0.05, p_value_tnn_gcm) 
map.contourf(xx, yy, p_value_tnn_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 13) 
map, xx, yy = basemap(lat, lon)
plt.title(u'M) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_dtr, levels=levs5, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_dtr_obs = ma.masked_where(p_value_dtr_obs >= 0.05, p_value_dtr_obs) 
map.contourf(xx, yy, p_value_dtr_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_dtr[0,:,:], levels=levs5, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_dtr_rcm = ma.masked_where(p_value_dtr_rcm >= 0.05, p_value_dtr_rcm) 
map.contourf(xx, yy, p_value_dtr_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_dtr, levels=levs5, latlon=True, cmap=cm.YlOrRd, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs5, drawedges=True, ax=ax)
cbar.set_label('DTR (°C)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_dtr_gcm = ma.masked_where(p_value_dtr_gcm >= 0.05, p_value_dtr_gcm) 
map.contourf(xx, yy, p_value_dtr_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 16)
map, xx, yy = basemap(lat, lon)
plt.title(u'P) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_su, levels=levs6, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_su_obs = ma.masked_where(p_value_su_obs >= 0.05, p_value_su_obs) 
map.contourf(xx, yy, p_value_su_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_su[0,:,:], levels=levs6, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_su_rcm = ma.masked_where(p_value_su_rcm >= 0.05, p_value_su_rcm) 
map.contourf(xx, yy, p_value_su_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_su, levels=levs6, latlon=True, cmap=cm.YlOrRd, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs6, drawedges=True, ax=ax)
cbar.set_label('SU (dias)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_su_gcm = ma.masked_where(p_value_su_gcm >= 0.05, p_value_su_gcm) 
map.contourf(xx, yy, p_value_su_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tr, levels=levs7, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tr_obs = ma.masked_where(p_value_tr_obs >= 0.05, p_value_tr_obs) 
map.contourf(xx, yy, p_value_tr_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tr[0,:,:], levels=levs7, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tr_rcm = ma.masked_where(p_value_tr_rcm >= 0.05, p_value_tr_rcm) 
map.contourf(xx, yy, p_value_tr_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tr, levels=levs7, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs7, drawedges=True, ax=ax)
cbar.set_label('TR (dias)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tr_gcm = ma.masked_where(p_value_tr_gcm >= 0.05, p_value_tr_gcm) 
map.contourf(xx, yy, p_value_tr_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tx10p, levels=levs8, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tx10p_obs = ma.masked_where(p_value_tx10p_obs >= 0.05, p_value_tx10p_obs) 
map.contourf(xx, yy, p_value_tx10p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tx10p[0,:,:], levels=levs8, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tx10p_rcm = ma.masked_where(p_value_tx10p_rcm >= 0.05, p_value_tx10p_rcm) 
map.contourf(xx, yy, p_value_tx10p_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tx10p, levels=levs8, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.set_label('TX10p (%)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tx10p_gcm = ma.masked_where(p_value_tx10p_gcm >= 0.05, p_value_tx10p_gcm) 
map.contourf(xx, yy, p_value_tx10p_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 25) 
map, xx, yy = basemap(lat, lon)
plt.title(u'Y) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tx90p, levels=levs8, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tx90p_obs = ma.masked_where(p_value_tx90p_obs >= 0.05, p_value_tx90p_obs) 
map.contourf(xx, yy, p_value_tx90p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tx90p[0,:,:], levels=levs8, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tx90p_rcm = ma.masked_where(p_value_tx90p_rcm >= 0.05, p_value_tx90p_rcm) 
map.contourf(xx, yy, p_value_tx90p_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tx90p, levels=levs8, latlon=True, cmap=cm.YlOrRd, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.set_label('TX90p (%)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tx90p_gcm = ma.masked_where(p_value_tx90p_gcm >= 0.05, p_value_tx90p_gcm) 
map.contourf(xx, yy, p_value_tx90p_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 28)
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1) CPC', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tn10p, levels=levs8, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tn10p_obs = ma.masked_where(p_value_tn10p_obs >= 0.05, p_value_tn10p_obs) 
map.contourf(xx, yy, p_value_tn10p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_tn10p[0,:,:], levels=levs8, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tn10p_rcm = ma.masked_where(p_value_tn10p_rcm >= 0.05, p_value_tn10p_rcm) 
map.contourf(xx, yy, p_value_tn10p_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_tn10p, levels=levs8, latlon=True, cmap=cm.YlOrRd, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.set_label('TN10p (%)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tn10p_gcm = ma.masked_where(p_value_tn10p_gcm >= 0.05, p_value_tn10p_gcm) 
map.contourf(xx, yy, p_value_tn10p_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 31) 
map, xx, yy = basemap(lat, lon)
plt.title(u'E.1) CPC', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_tn90p, levels=levs8, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tn90p_obs = ma.masked_where(p_value_tn90p_obs >= 0.05, p_value_tn90p_obs) 
map.contourf(xx, yy, p_value_tn90p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 32)
map, xx, yy = basemap(lat, lon)
plt.title(u'F.1) RegCM4.7', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, rcm_tn90p[0,:,:], levels=levs8, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tn90p_rcm = ma.masked_where(p_value_tn90p_rcm >= 0.05, p_value_tn90p_rcm) 
map.contourf(xx, yy, p_value_tn90p_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(11, 3, 33)
map, xx, yy = basemap(lat, lon)
plt.title(u'G.1) HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, gcm_tn90p, levels=levs8, latlon=True, cmap=cm.YlOrRd, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.set_label('TN90p (%)', size=6, labelpad=10, fontweight='bold', rotation=90)
cbar.ax.tick_params(labelsize=6) 
p_value_tn90p_gcm = ma.masked_where(p_value_tn90p_gcm >= 0.05, p_value_tn90p_gcm) 
map.contourf(xx, yy, p_value_tn90p_gcm, colors='none', hatches=['....'])

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_etccdi_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()	
	
	
	
	
	
