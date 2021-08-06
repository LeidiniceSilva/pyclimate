# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal climatology bias maps from regcm47 and hadgem models and obs database"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib as mpl 
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from comp_statist_indices import compute_av


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_obs = value[2:240:3,:,:]

	std_djf = np.std(season_obs[3:80:4], axis=0)
	std_mam = np.std(season_obs[0:80:4], axis=0)
	std_jja = np.std(season_obs[1:80:4], axis=0)
	std_son = np.std(season_obs[2:80:4], axis=0)
	std_annual = np.std(value[0:240:12,:,:], axis=0)
		
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)
	mam_obs = np.nanmean(season_obs[0:80:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:80:4], axis=0)
	son_obs = np.nanmean(season_obs[2:80:4], axis=0)
	annual_obs = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_annual, djf_obs, mam_obs, jja_obs, son_obs, annual_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]

	std_djf = np.std(season_rcm[3:80:4], axis=0)
	std_mam = np.std(season_rcm[0:80:4], axis=0)
	std_jja = np.std(season_rcm[1:80:4], axis=0)
	std_son = np.std(season_rcm[2:80:4], axis=0)
	std_annual = np.std(value[0:240:12,:,:], axis=0)
	
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)
	annual_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_annual, djf_rcm, mam_rcm, jja_rcm, son_rcm, annual_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]

	std_djf = np.std(season_gcm[3:80:4], axis=0)
	std_mam = np.std(season_gcm[0:80:4], axis=0)
	std_jja = np.std(season_gcm[1:80:4], axis=0)
	std_son = np.std(season_gcm[2:80:4], axis=0)
	std_annual = np.std(value[0:240:12,:,:], axis=0)
	
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(season_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(season_gcm[2:80:4], axis=0)
	annual_gcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_annual, djf_gcm, mam_gcm, jja_gcm, son_gcm, annual_gcm
	

def ttest(mean_sample1, mean_sample2, std_sample1, std_sample2):

	# Calculate t statistics	
	p1 = mean_sample1 - mean_sample2 
	p2 = (std_sample1 - std_sample2) / np.sqrt(240)

	ttest = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(ttest, df=240)

	return p_value
	
	
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
	

# Import models and obs database 
lat, lon, pre_std_djf_cru, pre_std_mam_cru, pre_std_jja_cru, pre_std_son_cru, pre_std_annual_cru, pre_djf_cru, pre_mam_cru, pre_jja_cru, pre_son_cru, pre_annual_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
lat, lon, pre_std_djf_rcm, pre_std_mam_rcm, pre_std_jja_rcm, pre_std_son_rcm, pre_std_annual_rcm, pre_djf_rcm, pre_mam_rcm, pre_jja_rcm, pre_son_rcm, pre_annual_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_gcm, pre_std_mam_gcm, pre_std_jja_gcm, pre_std_son_gcm, pre_std_annual_gcm, pre_djf_gcm, pre_mam_gcm, pre_jja_gcm, pre_son_gcm, pre_annual_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, tas_std_djf_cru, tas_std_mam_cru, tas_std_jja_cru, tas_std_son_cru, tas_std_annual_cru, tas_djf_cru, tas_mam_cru, tas_jja_cru, tas_son_cru, tas_annual_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_std_djf_rcm, tas_std_mam_rcm, tas_std_jja_rcm, tas_std_son_rcm, tas_std_annual_rcm, tas_djf_rcm, tas_mam_rcm, tas_jja_rcm, tas_son_rcm, tas_annual_rcm = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_gcm, tas_std_mam_gcm, tas_std_jja_gcm, tas_std_son_gcm, tas_std_annual_gcm, tas_djf_gcm, tas_mam_gcm, tas_jja_gcm, tas_son_gcm, tas_annual_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute bias from models and obs database 
pre_djf_rcm_bias = pre_djf_rcm - pre_djf_cru
pre_mam_rcm_bias = pre_mam_rcm - pre_mam_cru
pre_jja_rcm_bias = pre_jja_rcm - pre_jja_cru
pre_son_rcm_bias = pre_son_rcm - pre_son_cru
pre_annual_rcm_bias = pre_annual_rcm - pre_annual_cru

pre_djf_gcm_bias = pre_djf_gcm - pre_djf_cru
pre_mam_gcm_bias = pre_mam_gcm - pre_mam_cru
pre_jja_gcm_bias = pre_jja_gcm - pre_jja_cru
pre_son_gcm_bias = pre_son_gcm - pre_son_cru
pre_annual_gcm_bias = pre_annual_gcm - pre_annual_cru

tas_djf_rcm_bias = np.nanmean(tas_djf_rcm, axis=0) - tas_djf_cru
tas_mam_rcm_bias = np.nanmean(tas_mam_rcm, axis=0) - tas_mam_cru
tas_jja_rcm_bias = np.nanmean(tas_jja_rcm, axis=0) - tas_jja_cru
tas_son_rcm_bias = np.nanmean(tas_son_rcm, axis=0) - tas_son_cru
tas_annual_rcm_bias = np.nanmean(tas_annual_rcm, axis=0) - tas_annual_cru

tas_djf_gcm_bias = tas_djf_gcm - tas_djf_cru
tas_mam_gcm_bias = tas_mam_gcm - tas_mam_cru
tas_jja_gcm_bias = tas_jja_gcm - tas_jja_cru
tas_son_gcm_bias = tas_son_gcm - tas_son_cru
tas_annual_gcm_bias = tas_annual_gcm - tas_annual_cru

# Compute added value from models and obs database 
av_pre_djf = compute_av(pre_djf_gcm, pre_djf_rcm, pre_djf_cru)
av_pre_mam = compute_av(pre_mam_gcm, pre_mam_rcm, pre_mam_cru)
av_pre_jja = compute_av(pre_jja_gcm, pre_jja_rcm, pre_jja_cru)
av_pre_son = compute_av(pre_son_gcm, pre_son_rcm, pre_son_cru)
av_pre_annual = compute_av(pre_annual_gcm, pre_annual_rcm, pre_annual_cru)

av_tas_djf = compute_av(tas_djf_gcm, np.nanmean(tas_djf_rcm, axis=0), tas_djf_cru)
av_tas_mam = compute_av(tas_mam_gcm, np.nanmean(tas_mam_rcm, axis=0), tas_mam_cru)
av_tas_jja = compute_av(tas_jja_gcm, np.nanmean(tas_jja_rcm, axis=0), tas_jja_cru)
av_tas_son = compute_av(tas_son_gcm, np.nanmean(tas_son_rcm, axis=0), tas_son_cru)
av_tas_annual = compute_av(tas_annual_gcm, np.nanmean(tas_annual_rcm, axis=0), tas_annual_cru)

# Compute ttest from models and obs database 
p_value_pre_djf_rcm_cru = ttest(pre_djf_rcm, pre_djf_cru, pre_std_djf_rcm, pre_std_djf_cru)
p_value_pre_mam_rcm_cru = ttest(pre_mam_rcm, pre_mam_cru, pre_std_mam_rcm, pre_std_mam_cru)
p_value_pre_jja_rcm_cru = ttest(pre_jja_rcm, pre_jja_cru, pre_std_jja_rcm, pre_std_jja_cru)
p_value_pre_son_rcm_cru = ttest(pre_son_rcm, pre_son_cru, pre_std_son_rcm, pre_std_son_cru)
p_value_pre_annual_rcm_cru = ttest(pre_annual_rcm, pre_annual_cru, pre_std_annual_rcm, pre_std_annual_cru)

p_value_pre_djf_gcm_cru = ttest(pre_djf_gcm, pre_djf_cru, pre_std_djf_gcm, pre_std_djf_cru)
p_value_pre_mam_gcm_cru = ttest(pre_mam_gcm, pre_mam_cru, pre_std_mam_gcm, pre_std_mam_cru)
p_value_pre_jja_gcm_cru = ttest(pre_jja_gcm, pre_jja_cru, pre_std_jja_gcm, pre_std_jja_cru)
p_value_pre_son_gcm_cru = ttest(pre_son_gcm, pre_son_cru, pre_std_son_gcm, pre_std_son_cru)
p_value_pre_annual_gcm_cru = ttest(pre_annual_gcm, pre_annual_cru, pre_std_annual_gcm, pre_std_annual_cru)

p_value_tas_djf_rcm_cru = ttest(tas_djf_rcm, tas_djf_cru, tas_std_djf_rcm, tas_std_djf_cru)
p_value_tas_mam_rcm_cru = ttest(tas_mam_rcm, tas_mam_cru, tas_std_mam_rcm, tas_std_mam_cru)
p_value_tas_jja_rcm_cru = ttest(tas_jja_rcm, tas_jja_cru, tas_std_jja_rcm, tas_std_jja_cru)
p_value_tas_son_rcm_cru = ttest(tas_son_rcm, tas_son_cru, tas_std_son_rcm, tas_std_son_cru)
p_value_tas_annual_rcm_cru = ttest(tas_annual_rcm, tas_annual_cru, tas_std_annual_rcm, tas_std_annual_cru)

p_value_tas_djf_gcm_cru = ttest(tas_djf_gcm, tas_djf_cru, tas_std_djf_gcm, tas_std_djf_cru)
p_value_tas_mam_gcm_cru = ttest(tas_mam_gcm, tas_mam_cru, tas_std_mam_gcm, tas_std_mam_cru)
p_value_tas_jja_gcm_cru = ttest(tas_jja_gcm, tas_jja_cru, tas_std_jja_gcm, tas_std_jja_cru)
p_value_tas_son_gcm_cru = ttest(tas_son_gcm, tas_son_cru, tas_std_son_gcm, tas_std_son_cru)
p_value_tas_annual_gcm_cru = ttest(tas_annual_gcm, tas_annual_cru, tas_std_annual_gcm, tas_std_annual_cru)

# Plot models and obs database 
fig = plt.figure(figsize=(8, 6))
levs1 = [-6, -3, -1, 1, 3, 6]
levs2 = [-1, -0.7, -0.4, 0.4, 0.7, 1]
	
ax = fig.add_subplot(6, 5, 1)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
plt_maps_bias = map.contourf(xx, yy, pre_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_rcm_cru = ma.masked_where(p_value_pre_djf_rcm_cru >= 0.1, p_value_pre_djf_rcm_cru) 
map.contourf(xx, yy, p_value_pre_djf_rcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 2)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_rcm_cru = ma.masked_where(p_value_pre_mam_rcm_cru >= 0.1, p_value_pre_mam_rcm_cru) 
map.contourf(xx, yy, p_value_pre_mam_rcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 3)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_rcm_cru = ma.masked_where(p_value_pre_jja_rcm_cru >= 0.1, p_value_pre_jja_rcm_cru) 
map.contourf(xx, yy, p_value_pre_jja_rcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 4)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_rcm_cru = ma.masked_where(p_value_pre_son_rcm_cru >= 0.1, p_value_pre_son_rcm_cru) 
map.contourf(xx, yy, p_value_pre_son_rcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 5)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_annual_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_rcm_cru = ma.masked_where(p_value_pre_annual_rcm_cru >= 0.1, p_value_pre_annual_rcm_cru) 
map.contourf(xx, yy, p_value_pre_annual_rcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 6)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_gcm_cru = ma.masked_where(p_value_pre_djf_gcm_cru >= 0.1, p_value_pre_djf_gcm_cru) 
map.contourf(xx, yy, p_value_pre_djf_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 7)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_gcm_cru = ma.masked_where(p_value_pre_mam_gcm_cru >= 0.1, p_value_pre_mam_gcm_cru) 
map.contourf(xx, yy, p_value_pre_mam_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 8)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_gcm_cru = ma.masked_where(p_value_pre_jja_gcm_cru >= 0.1, p_value_pre_jja_gcm_cru) 
map.contourf(xx, yy, p_value_pre_jja_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 9)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_gcm_cru = ma.masked_where(p_value_pre_son_gcm_cru >= 0.1, p_value_pre_son_gcm_cru) 
map.contourf(xx, yy, p_value_pre_son_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 10)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_annual_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 11)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_djf, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 12)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_mam, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 13)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_jja, levels=levs2, latlon=True, cmap=cm.PiYG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 14)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_son, levels=levs2, latlon=True, cmap=cm.PiYG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 15)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_annual, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 16)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_cru = ma.masked_where(p_value_tas_djf_rcm_cru >= 0.1, p_value_tas_djf_rcm_cru) 
map.contourf(xx, yy, p_value_tas_djf_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 17)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_rcm_cru = ma.masked_where(p_value_tas_mam_rcm_cru >= 0.1, p_value_tas_mam_rcm_cru) 
map.contourf(xx, yy, p_value_tas_mam_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 18)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_cru = ma.masked_where(p_value_tas_jja_rcm_cru >= 0.1, p_value_tas_jja_rcm_cru) 
map.contourf(xx, yy, p_value_tas_jja_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 19)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_rcm_cru = ma.masked_where(p_value_tas_son_rcm_cru >= 0.1, p_value_tas_son_rcm_cru) 
map.contourf(xx, yy, p_value_tas_son_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 20)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_annual_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)    
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_rcm_cru = ma.masked_where(p_value_tas_annual_rcm_cru >= 0.1, p_value_tas_annual_rcm_cru) 
map.contourf(xx, yy, p_value_tas_annual_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 21)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm_cru = ma.masked_where(p_value_tas_djf_gcm_cru >= 0.1, p_value_tas_djf_gcm_cru) 
map.contourf(xx, yy, p_value_tas_djf_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 22)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_gcm_cru = ma.masked_where(p_value_tas_mam_gcm_cru >= 0.1, p_value_tas_mam_gcm_cru) 
map.contourf(xx, yy, p_value_tas_mam_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 23)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm_cru = ma.masked_where(p_value_tas_jja_gcm_cru >= 0.1, p_value_tas_jja_gcm_cru) 
map.contourf(xx, yy, p_value_tas_jja_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 24)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_gcm_cru = ma.masked_where(p_value_tas_son_gcm_cru >= 0.1, p_value_tas_son_gcm_cru) 
map.contourf(xx, yy, p_value_tas_son_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 25)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_annual_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_gcm_cru = ma.masked_where(p_value_tas_annual_gcm_cru >= 0.1, p_value_tas_annual_gcm_cru) 
map.contourf(xx, yy, p_value_tas_annual_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(6, 5, 26)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_djf, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 27)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_mam, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 28)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_jja, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 29)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_son, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 30)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_annual, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



