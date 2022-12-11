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

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch
from comp_statist_indices import compute_av


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documentos/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_prcptot': u'precip', 
	u'eca_r95': u'precip',
	u'eca_r99': u'precip', 
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
	
	return lat, lon, obs, std
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documentos/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95': u'pr',
	u'eca_r99': u'pr', 
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

	return lat, lon, rcm, std


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documentos/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95': u'pr',
	u'eca_r99': u'pr', 
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
	
	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)
	
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
lat, lon, obs_prcptot, obs_prcptot_std = import_obs('eca_prcptot', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_prcptot, rcm_prcptot_std = import_rcm('eca_prcptot', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_prcptot, gcm_prcptot_std = import_gcm('eca_prcptot', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r95p, obs_r95p_std = import_obs('eca_r95', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r95p, rcm_r95p_std = import_rcm('eca_r95', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r95p, gcm_r95p_std = import_gcm('eca_r95', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r99p, obs_r99p_std = import_obs('eca_r99', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r99p, rcm_r99p_std = import_rcm('eca_r99', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r99p, gcm_r99p_std = import_gcm('eca_r99', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_rx1day, obs_rx1day_std = import_obs('eca_rx1day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_rx1day, rcm_rx1day_std = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx1day, gcm_rx1day_std = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_rx5day, obs_rx5day_std = import_obs('eca_rx5day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_rx5day, rcm_rx5day_std = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_rx5day, gcm_rx5day_std = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_sdii, obs_sdii_std = import_obs('eca_sdii', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_sdii, rcm_sdii_std = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_sdii, gcm_sdii_std = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_cdd, obs_cdd_std = import_obs('eca_cdd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_cdd, rcm_cdd_std = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cdd, gcm_cdd_std = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_cwd, obs_cwd_std = import_obs('eca_cwd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_cwd, rcm_cwd_std = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cwd, gcm_cwd_std = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r10mm, obs_r10mm_std = import_obs('eca_r10mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r10mm, rcm_r10mm_std = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r10mm, gcm_r10mm_std = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r20mm, obs_r20mm_std = import_obs('eca_r20mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
lat, lon, rcm_r20mm, rcm_r20mm_std = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r20mm, gcm_r20mm_std = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Compute bias from extreme indices 
bias_rcm_prcptot = rcm_prcptot - obs_prcptot
bias_rcm_r95p = rcm_r95p - obs_r95p
bias_rcm_r99p = rcm_r99p - obs_r99p
bias_rcm_rx1day = rcm_rx1day - obs_rx1day
bias_rcm_rx5day = rcm_rx5day - obs_rx5day
bias_rcm_sdii = rcm_sdii - obs_sdii
bias_rcm_cdd = rcm_cdd - obs_cdd
bias_rcm_cwd = rcm_cwd - obs_cwd
bias_rcm_r10mm = rcm_r10mm - obs_r10mm
bias_rcm_r20mm = rcm_r20mm - obs_r20mm

bias_gcm_prcptot = gcm_prcptot - obs_prcptot
bias_gcm_r95p = gcm_r95p - obs_r95p
bias_gcm_r99p = gcm_r99p - obs_r99p
bias_gcm_rx1day = gcm_rx1day - obs_rx1day
bias_gcm_rx5day = gcm_rx5day - obs_rx5day
bias_gcm_sdii = gcm_sdii - obs_sdii
bias_gcm_cdd = gcm_cdd - obs_cdd
bias_gcm_cwd = gcm_cwd - obs_cwd
bias_gcm_r10mm = gcm_r10mm - obs_r10mm
bias_gcm_r20mm = gcm_r20mm - obs_r20mm

av_prcptot = compute_av(gcm_prcptot, rcm_prcptot, obs_prcptot)
av_r95p = compute_av(gcm_r95p, rcm_r95p, obs_r95p)
av_r99p = compute_av(gcm_r99p, rcm_r99p, obs_r99p)
av_rx1day = compute_av(gcm_rx1day, rcm_rx1day, obs_rx1day)
av_rx5day = compute_av(gcm_rx5day, rcm_rx5day, obs_rx5day)
av_sdii = compute_av(gcm_sdii, rcm_sdii, obs_sdii)
av_cdd = compute_av(gcm_cdd, rcm_cdd, obs_cdd)
av_cwd = compute_av(gcm_cwd, rcm_cwd, obs_cwd)
av_r10mm = compute_av(gcm_r10mm, rcm_r10mm, obs_r10mm)
av_r20mm = compute_av(gcm_r20mm, rcm_r20mm, obs_r20mm)

# Compute ttest from models and obs database 
p_value_rcm_prcptot = ttest(rcm_prcptot, obs_prcptot, rcm_prcptot_std, obs_prcptot_std)
p_value_rcm_r95p = ttest(rcm_r95p, obs_r95p, rcm_r95p_std, obs_r95p_std)
p_value_rcm_r99p = ttest(rcm_r99p, obs_r99p, rcm_r99p_std, obs_r99p_std)
p_value_rcm_rx1day = ttest(rcm_rx1day, obs_rx1day, rcm_rx1day_std, obs_rx1day_std)
p_value_rcm_rx5day = ttest(rcm_rx5day, obs_rx5day, rcm_rx5day_std, obs_rx5day_std)
p_value_rcm_sdii = ttest(rcm_sdii, obs_sdii, rcm_sdii_std, obs_sdii_std)
p_value_rcm_cdd = ttest(rcm_cdd, obs_cdd, rcm_cdd_std, obs_cdd_std)
p_value_rcm_cwd = ttest(rcm_cwd, obs_cwd, rcm_cwd_std, obs_cwd_std)
p_value_rcm_r10mm = ttest(rcm_r10mm, obs_r10mm, rcm_r10mm_std, obs_r10mm_std)
p_value_rcm_r20mm = ttest(rcm_r20mm, obs_r20mm, rcm_r20mm_std, obs_r20mm_std)

p_value_gcm_prcptot = ttest(gcm_prcptot, obs_prcptot, gcm_prcptot_std, obs_prcptot_std)
p_value_gcm_r95p = ttest(gcm_r95p, obs_r95p, gcm_r95p_std, obs_r95p_std)
p_value_gcm_r99p = ttest(gcm_r99p, obs_r99p, gcm_r99p_std, obs_r99p_std)
p_value_gcm_rx1day = ttest(gcm_rx1day, obs_rx1day, gcm_rx1day_std, obs_rx1day_std)
p_value_gcm_rx5day = ttest(gcm_rx5day, obs_rx5day, gcm_rx5day_std, obs_rx5day_std)
p_value_gcm_sdii = ttest(gcm_sdii, obs_sdii, gcm_sdii_std, obs_sdii_std)
p_value_gcm_cdd = ttest(gcm_cdd, obs_cdd, gcm_cdd_std, obs_cdd_std)
p_value_gcm_cwd = ttest(gcm_cwd, obs_cwd, gcm_cwd_std, obs_cwd_std)
p_value_gcm_r10mm = ttest(gcm_r10mm, obs_r10mm, gcm_r10mm_std, obs_r10mm_std)
p_value_gcm_r20mm = ttest(gcm_r20mm, obs_r20mm, gcm_r20mm_std, obs_r20mm_std)

# Plot extreme indices 
fig = plt.figure(figsize=(6, 10))
levs1 = [-600, -400, -200, 200, 400, 600]
levs2 = [-25, -15, -5, 5, 15, 25]
levs3 = [-25, -15, -5, 5, 15, 25]
levs4 = [-40, -30, -10, 10, 30, 40]
levs5 = [-75, -50, -25, 25, 50, 75]
levs6 = [-6, -4, -2, 2, 4, 6]
levs7 = [-75, -50, -25, 25, 50, 75]
levs8 = [-75, -50, -25, 25, 50, 75]
levs9 = [-75, -50, -25, 25, 50, 75]
levs10 = [-25, -15, -5, 5, 15, 25]
levs11 = [-1, -0.5, -0.1, 0.1, 0.5, 1]

ax = fig.add_subplot(10, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_prcptot, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_prcptot = ma.masked_where(p_value_rcm_prcptot >= 0.05, p_value_rcm_prcptot) 
map.contourf(xx, yy, p_value_rcm_prcptot, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_prcptot, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_prcptot = ma.masked_where(p_value_gcm_prcptot >= 0.05, p_value_gcm_prcptot) 
map.contourf(xx, yy, p_value_gcm_prcptot, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) AV PRCPTOT', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_prcptot, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(10, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_r95p, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_r95p = ma.masked_where(p_value_rcm_r95p >= 0.05, p_value_rcm_r95p) 
map.contourf(xx, yy, p_value_rcm_r95p, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_r95p, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_r95p = ma.masked_where(p_value_gcm_r95p >= 0.05, p_value_gcm_r95p) 
map.contourf(xx, yy, p_value_gcm_r95p, colors='none', hatches=['....'])
		
ax = fig.add_subplot(10, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) AV R95p', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_r95p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_r99p, levels=levs3, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_r99p = ma.masked_where(p_value_rcm_r99p >= 0.05, p_value_rcm_r99p) 
map.contourf(xx, yy, p_value_rcm_r99p, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_r99p, levels=levs3, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_r99p = ma.masked_where(p_value_gcm_r99p >= 0.05, p_value_gcm_r99p) 
map.contourf(xx, yy, p_value_gcm_r99p, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) AV R99p', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_r99p, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_rx1day, levels=levs4, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_rx1day = ma.masked_where(p_value_rcm_rx1day >= 0.05, p_value_rcm_rx1day) 
map.contourf(xx, yy, p_value_rcm_rx1day, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_rx1day, levels=levs4, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], BrBG=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_rx1day = ma.masked_where(p_value_gcm_rx1day >= 0.05, p_value_gcm_rx1day) 
map.contourf(xx, yy, p_value_gcm_rx1day, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) AV Rx1day', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_rx1day, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 13) 
map, xx, yy = basemap(lat, lon)
plt.title(u'M) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_rx5day, levels=levs5, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_rx5day = ma.masked_where(p_value_rcm_rx5day >= 0.05, p_value_rcm_rx5day) 
map.contourf(xx, yy, p_value_rcm_rx5day, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_rx5day, levels=levs5, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs5, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_rx5day = ma.masked_where(p_value_gcm_rx5day >= 0.05, p_value_gcm_rx5day) 
map.contourf(xx, yy, p_value_gcm_rx5day, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O) AV Rx5day', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_rx5day, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 16) 
map, xx, yy = basemap(lat, lon)
plt.title(u'P) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_sdii, levels=levs6, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_sdii = ma.masked_where(p_value_rcm_sdii >= 0.05, p_value_rcm_sdii) 
map.contourf(xx, yy, p_value_rcm_sdii, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_sdii, levels=levs6, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs6, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_sdii = ma.masked_where(p_value_gcm_sdii >= 0.05, p_value_gcm_sdii) 
map.contourf(xx, yy, p_value_gcm_sdii, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R) AV SDII', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_sdii, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 19)
map, xx, yy = basemap(lat, lon)
plt.title(u'S) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_cdd, levels=levs7, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_cdd = ma.masked_where(p_value_rcm_cdd >= 0.05, p_value_rcm_cdd) 
map.contourf(xx, yy, p_value_rcm_cdd, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 20)
map, xx, yy = basemap(lat, lon)
plt.title(u'T) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_cdd, levels=levs7, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs7, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_cdd = ma.masked_where(p_value_gcm_cdd >= 0.05, p_value_gcm_cdd) 
map.contourf(xx, yy, p_value_gcm_cdd, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 21)
map, xx, yy = basemap(lat, lon)
plt.title(u'U) AV CDD', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_cdd, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 22)
map, xx, yy = basemap(lat, lon)
plt.title(u'V) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_cwd, levels=levs8, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_cwd = ma.masked_where(p_value_rcm_cwd >= 0.05, p_value_rcm_cwd) 
map.contourf(xx, yy, p_value_rcm_cwd, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 23)
map, xx, yy = basemap(lat, lon)
plt.title(u'W) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_cwd, levels=levs8, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs8, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_cwd = ma.masked_where(p_value_gcm_cwd >= 0.05, p_value_gcm_cwd) 
map.contourf(xx, yy, p_value_gcm_cwd, colors='none', hatches=['....'])
		
ax = fig.add_subplot(10, 3, 24)
map, xx, yy = basemap(lat, lon)
plt.title(u'X) AV CWD', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_cwd, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 25)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
map.contourf(xx, yy, bias_rcm_r10mm, levels=levs9, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_r10mm = ma.masked_where(p_value_rcm_r10mm >= 0.05, p_value_rcm_r10mm) 
map.contourf(xx, yy, p_value_rcm_r10mm, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 26)
map, xx, yy = basemap(lat, lon)
plt.title(u'Z) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, bias_gcm_r10mm, levels=levs9, latlon=True, cmap=cm.BrBG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs9, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_r10mm = ma.masked_where(p_value_gcm_r10mm >= 0.05, p_value_gcm_r10mm) 
map.contourf(xx, yy, p_value_gcm_r10mm, colors='none', hatches=['....'])

ax = fig.add_subplot(10, 3, 27)
map, xx, yy = basemap(lat, lon)
plt.title(u'A.1) AV R10mm', loc='left', fontsize=8, fontweight='bold')
map.contourf(xx, yy, av_r10mm, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(10, 3, 28) 
map, xx, yy = basemap(lat, lon)
plt.title(u'B.1) MBE RegCM 4.7', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, bias_rcm_r20mm, levels=levs10, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_rcm_r20mm = ma.masked_where(p_value_rcm_r20mm >= 0.05, p_value_rcm_r20mm) 
map.contourf(xx, yy, p_value_rcm_r20mm, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 29)
map, xx, yy = basemap(lat, lon)
plt.title(u'C.1) MBE HadGEM2-ES', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, bias_gcm_r20mm, levels=levs10, latlon=True, cmap=cm.BrBG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs10, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
p_value_gcm_r20mm = ma.masked_where(p_value_gcm_r20mm >= 0.05, p_value_gcm_r20mm) 
map.contourf(xx, yy, p_value_gcm_r20mm, colors='none', hatches=['....'])
	
ax = fig.add_subplot(10, 3, 30)
map, xx, yy = basemap(lat, lon)
plt.title(u'D.1) AV R20mm', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
map.contourf(xx, yy, av_r20mm, levels=levs11, latlon=True, cmap=cm.PiYG, extend='both') 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs11, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30, hspace=0.30)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_etccdi_pre_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()	
	
	
	
	
	
