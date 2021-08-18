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

	dict_var = {u'eca_prcptot': u'precip', 
	u'eca_r95p': u'precip',
	u'eca_r99p': u'precip', 
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
	obs  = np.nanmean(var[:][:,:,:], axis=0)

	std   = np.std(var[:][:,:,:], axis=0)
	sample  = np.nanmean(var[:][0:2,:,:], axis=0)
	
	return lat, lon, obs, std, sample
	
	
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

	std   = np.std(var[:][:,:,:], axis=0)
	sample  = np.nanmean(var[:][0:2,:,:], axis=0)

	return lat, lon, rcm, std, sample


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
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy


def ttest(mean, std, sample):

	# Calculate t statistics
	p1 = mean - sample
	p2= std / np.sqrt(240)
	p3 = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(p3, df=239)
	
	return p_value
		
	
# Import extreme indices 
lat, lon, obs_prcptot, obs_prcptot_std, obs_prcptot_sample = import_obs('eca_prcptot', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_prcptot, rcm_prcptot_std, rcm_prcptot_sample = import_rcm('eca_prcptot', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_prcptot, gcm_prcptot_std, gcm_prcptot_sample = import_gcm('eca_prcptot', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r95p, obs_r95p_std, obs_r95p_sample = import_obs('eca_r95p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r95p, rcm_r95p_std, rcm_r95p_sample = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r95p, gcm_r95p_std, gcm_r95p_sample = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r99p, obs_r99p_std, obs_r99p_sample = import_obs('eca_r99p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r99p, rcm_r99p_std, rcm_r99p_sample = import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r99p, gcm_r99p_std, gcm_r99p_sample = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_rx1day, obs_rx1day_std, obs_rx1day_sample = import_obs('eca_rx1day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_rx1day, rcm_rx1day_std, rcm_rx1day_sample = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx1day, gcm_rx1day_std, gcm_rx1day_sample = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_rx5day, obs_rx5day_std, obs_rx5day_sample = import_obs('eca_rx5day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_rx5day, rcm_rx5day_std, rcm_rx5day_sample = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx5day, gcm_rx5day_std, gcm_rx5day_sample = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_sdii, obs_sdii_std, obs_sdii_sample = import_obs('eca_sdii', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_sdii, rcm_sdii_std, rcm_sdii_sample = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_sdii, gcm_sdii_std, gcm_sdii_sample = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_cdd, obs_cdd_std, obs_cdd_sample = import_obs('eca_cdd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_cdd, rcm_cdd_std, rcm_cdd_sample = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cdd, gcm_cdd_std, gcm_cdd_sample = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_cwd, obs_cwd_std, obs_cwd_sample = import_obs('eca_cwd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_cwd, rcm_cwd_std, rcm_cwd_sample = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cwd, gcm_cwd_std, gcm_cwd_sample = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r10mm, obs_r10mm_std, obs_r10mm_sample = import_obs('eca_r10mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r10mm, rcm_r10mm_std, rcm_r10mm_sample = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r10mm, gcm_r10mm_std, gcm_r10mm_sample = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r20mm, obs_r20mm_std, obs_r20mm_sample = import_obs('eca_r20mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r20mm, rcm_r20mm_std, rcm_r20mm_sample = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r20mm, gcm_r20mm_std, gcm_r20mm_sample = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Compute ttest from models and obs database 
p_value_prcptot_obs = ttest(obs_prcptot, obs_prcptot_std, obs_prcptot_sample)
p_value_r95p_obs = ttest(obs_r95p, obs_r95p_std, obs_r95p_sample)
p_value_r99p_obs = ttest(obs_r99p, obs_r99p_std, obs_r99p_sample)
p_value_rx1day_obs = ttest(obs_rx1day, obs_rx1day_std, obs_rx1day_sample)
p_value_rx5day_obs = ttest(obs_rx5day, obs_rx5day_std, obs_rx5day_sample)
p_value_sdii_obs = ttest(obs_sdii, obs_sdii_std, obs_sdii_sample)
p_value_cdd_obs = ttest(obs_cdd, obs_cdd_std, obs_cdd_sample)
p_value_cwd_obs = ttest(obs_cwd, obs_cwd_std, obs_cwd_sample)
p_value_r10mm_obs = ttest(obs_r10mm, obs_r10mm_std, obs_r10mm_sample)
p_value_r20mm_obs = ttest(obs_r20mm, obs_r20mm_std, obs_r20mm_sample)

p_value_prcptot_rcm = ttest(rcm_prcptot, rcm_prcptot_std, rcm_prcptot_sample)
p_value_r95p_rcm = ttest(rcm_r95p, rcm_r95p_std, rcm_r95p_sample)
p_value_r99p_rcm = ttest(rcm_r99p, rcm_r99p_std, rcm_r99p_sample)
p_value_rx1day_rcm = ttest(rcm_rx1day, rcm_rx1day_std, rcm_rx1day_sample)
p_value_rx5day_rcm = ttest(rcm_rx5day, rcm_rx5day_std, rcm_rx5day_sample)
p_value_sdii_rcm = ttest(rcm_sdii, rcm_sdii_std, rcm_sdii_sample)
p_value_cdd_rcm = ttest(rcm_cdd, rcm_cdd_std, rcm_cdd_sample)
p_value_cwd_rcm = ttest(rcm_cwd, rcm_cwd_std, rcm_cwd_sample)
p_value_r10mm_rcm = ttest(rcm_r10mm, rcm_r10mm_std, rcm_r10mm_sample)
p_value_r20mm_rcm = ttest(rcm_r20mm, rcm_r20mm_std, rcm_r20mm_sample)

p_value_prcptot_gcm = ttest(gcm_prcptot, gcm_prcptot_std, gcm_prcptot_sample)
p_value_r95p_gcm = ttest(gcm_r95p, gcm_r95p_std, gcm_r95p_sample)
p_value_r99p_gcm = ttest(gcm_r99p, gcm_r99p_std, gcm_r99p_sample)
p_value_rx1day_gcm = ttest(gcm_rx1day, gcm_rx1day_std, gcm_rx1day_sample)
p_value_rx5day_gcm = ttest(gcm_rx5day, gcm_rx5day_std, gcm_rx5day_sample)
p_value_sdii_gcm = ttest(gcm_sdii, gcm_sdii_std, gcm_sdii_sample)
p_value_cdd_gcm = ttest(gcm_cdd, gcm_cdd_std, gcm_cdd_sample)
p_value_cwd_gcm = ttest(gcm_cwd, gcm_cwd_std, gcm_cwd_sample)
p_value_r10mm_gcm = ttest(gcm_r10mm, gcm_r10mm_std, gcm_r10mm_sample)
p_value_r20mm_gcm = ttest(gcm_r20mm, gcm_r20mm_std, gcm_r20mm_sample)

# Plot extreme indices 
fig = plt.figure(figsize=(6, 10))

levs1 = [100, 500, 1000, 1500, 2000, 3000]
levs2 = [5, 10, 20, 40, 60, 80]
levs3 = [1, 5, 10, 20, 30, 40]
levs4 = [1, 10, 20, 30, 40, 50]
levs5 = [10, 40, 70, 100, 130, 160]
levs6 = [1, 3, 6, 9, 12, 15]
levs7 = [10, 40, 70, 100, 130, 160]
levs8 = [1, 10, 20, 40, 60, 80]
levs9 = [5, 20, 40, 60, 80, 100]
levs10 = [1, 5, 10, 20, 30, 40]

ax = fig.add_subplot(10, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_prcptot_obs = ma.masked_where(p_value_prcptot_obs >= 0.05, p_value_prcptot_obs) 
map.contourf(xx, yy, p_value_prcptot_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_prcptot_rcm = ma.masked_where(p_value_prcptot_rcm >= 0.05, p_value_prcptot_rcm) 
map.contourf(xx, yy, p_value_prcptot_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)  
p_value_prcptot_gcm = ma.masked_where(p_value_prcptot_gcm >= 0.05, p_value_prcptot_gcm) 
map.contourf(xx, yy, p_value_prcptot_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_r95p, levels=levs2, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_r95p_obs = ma.masked_where(p_value_r95p_obs >= 0.05, p_value_r95p_obs) 
map.contourf(xx, yy, p_value_r95p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_r95p, levels=levs2, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_r95p_rcm = ma.masked_where(p_value_r95p_rcm >= 0.05, p_value_r95p_rcm) 
map.contourf(xx, yy, p_value_r95p_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_r95p, levels=levs2, latlon=True, cmap=cm.YlGnBu, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_r95p_gcm = ma.masked_where(p_value_r95p_gcm >= 0.05, p_value_r95p_gcm) 
map.contourf(xx, yy, p_value_r95p_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_r99p, levels=levs3, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_r99p_obs = ma.masked_where(p_value_r99p_obs >= 0.05, p_value_r99p_obs) 
map.contourf(xx, yy, p_value_r99p_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_r99p, levels=levs3, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_r99p_rcm = ma.masked_where(p_value_r99p_rcm >= 0.05, p_value_r99p_rcm) 
map.contourf(xx, yy, p_value_r99p_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_r99p, levels=levs3, latlon=True, cmap=cm.YlGnBu, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_r99p_gcm = ma.masked_where(p_value_r99p_gcm >= 0.05, p_value_r99p_gcm) 
map.contourf(xx, yy, p_value_r99p_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_rx1day, levels=levs4, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rx1day_obs = ma.masked_where(p_value_rx1day_obs >= 0.05, p_value_rx1day_obs) 
map.contourf(xx, yy, p_value_rx1day_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_rx1day, levels=levs4, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_rx1day_rcm = ma.masked_where(p_value_rx1day_rcm >= 0.05, p_value_rx1day_rcm) 
map.contourf(xx, yy, p_value_rx1day_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_rx1day, levels=levs4, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_rx1day_gcm = ma.masked_where(p_value_rx1day_gcm >= 0.05, p_value_rx1day_gcm) 
map.contourf(xx, yy, p_value_rx1day_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 13) 
map, xx, yy = basemap(lat, lon)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_rx5day, levels=levs5, latlon=True, cmap=cm.YlGnBu, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rx5day_obs = ma.masked_where(p_value_rx5day_obs >= 0.05, p_value_rx5day_obs) 
map.contourf(xx, yy, p_value_rx5day_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_rx5day, levels=levs5, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_rx5day_rcm = ma.masked_where(p_value_rx5day_rcm >= 0.05, p_value_rx5day_rcm) 
map.contourf(xx, yy, p_value_rx5day_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_rx5day, levels=levs5, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs5, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_rx5day_gcm = ma.masked_where(p_value_rx5day_gcm >= 0.05, p_value_rx5day_gcm) 
map.contourf(xx, yy, p_value_rx5day_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 16) 
map, xx, yy = basemap(lat, lon)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_sdii, levels=levs6, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_sdii_obs = ma.masked_where(p_value_sdii_obs >= 0.05, p_value_sdii_obs) 
map.contourf(xx, yy, p_value_sdii_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_sdii, levels=levs6, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_sdii_rcm = ma.masked_where(p_value_sdii_rcm >= 0.05, p_value_sdii_rcm) 
map.contourf(xx, yy, p_value_sdii_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_sdii, levels=levs6, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs6, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_sdii_gcm = ma.masked_where(p_value_sdii_gcm >= 0.05, p_value_sdii_gcm) 
map.contourf(xx, yy, p_value_sdii_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_cdd, levels=levs7, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_cdd_obs = ma.masked_where(p_value_cdd_obs >= 0.05, p_value_cdd_obs) 
map.contourf(xx, yy, p_value_cdd_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_cdd, levels=levs7, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_cdd_rcm = ma.masked_where(p_value_cdd_rcm >= 0.05, p_value_cdd_rcm) 
map.contourf(xx, yy, p_value_cdd_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_cdd, levels=levs7, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs7, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_cdd_gcm = ma.masked_where(p_value_cdd_gcm >= 0.05, p_value_cdd_gcm) 
map.contourf(xx, yy, p_value_cdd_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_cwd, levels=levs8, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_cwd_obs = ma.masked_where(p_value_cwd_obs >= 0.05, p_value_cwd_obs) 
map.contourf(xx, yy, p_value_cwd_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_cwd, levels=levs8, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_cwd_rcm = ma.masked_where(p_value_cwd_rcm >= 0.05, p_value_cwd_rcm) 
map.contourf(xx, yy, p_value_cwd_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_cwd, levels=levs8, latlon=True, cmap=cm.YlGnBu, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_cwd_gcm = ma.masked_where(p_value_cwd_gcm >= 0.05, p_value_cwd_gcm) 
map.contourf(xx, yy, p_value_cwd_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 25)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, obs_r10mm, levels=levs9, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_r10mm_obs = ma.masked_where(p_value_r10mm_obs >= 0.05, p_value_r10mm_obs) 
map.contourf(xx, yy, p_value_r10mm_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, rcm_r10mm, levels=levs9, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_r10mm_rcm = ma.masked_where(p_value_r10mm_rcm >= 0.05, p_value_r10mm_rcm) 
map.contourf(xx, yy, p_value_r10mm_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, gcm_r10mm, levels=levs9, latlon=True, cmap=cm.YlGnBu, extend='max')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs9, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_r10mm_gcm = ma.masked_where(p_value_r10mm_gcm >= 0.05, p_value_r10mm_gcm) 
map.contourf(xx, yy, p_value_r10mm_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 28) 
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, obs_r20mm, levels=levs10, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_r20mm_obs = ma.masked_where(p_value_r20mm_obs >= 0.05, p_value_r20mm_obs) 
map.contourf(xx, yy, p_value_r20mm_obs, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, rcm_r20mm, levels=levs10, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_r20mm_rcm = ma.masked_where(p_value_r20mm_rcm >= 0.05, p_value_r20mm_rcm) 
map.contourf(xx, yy, p_value_r20mm_rcm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, gcm_r20mm, levels=levs10, latlon=True, cmap=cm.YlGnBu, extend='max') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs10, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_r20mm_gcm = ma.masked_where(p_value_r20mm_gcm >= 0.05, p_value_r20mm_gcm) 
map.contourf(xx, yy, p_value_r20mm_gcm, colors='none', hatches=['....'])

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_etccdi_pre_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()	
	
	
	
	
	
