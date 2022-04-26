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
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	sea_obs = value[2:240:3,:,:]

	std_djf = np.std(sea_obs[3:80:4], axis=0)
	std_mam = np.std(sea_obs[0:80:4], axis=0)
	std_jja = np.std(sea_obs[1:80:4], axis=0)
	std_son = np.std(sea_obs[2:80:4], axis=0)
	std_ann = np.std(value[0:240:12,:,:], axis=0)
		
	djf_obs = np.nanmean(sea_obs[3:80:4], axis=0)
	mam_obs = np.nanmean(sea_obs[0:80:4], axis=0)
	jja_obs = np.nanmean(sea_obs[1:80:4], axis=0)
	son_obs = np.nanmean(sea_obs[2:80:4], axis=0)
	ann_obs = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_ann, djf_obs, mam_obs, jja_obs, son_obs, ann_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	sea_rcm = value[2:240:3,:,:]

	std_djf = np.std(sea_rcm[3:80:4], axis=0)
	std_mam = np.std(sea_rcm[0:80:4], axis=0)
	std_jja = np.std(sea_rcm[1:80:4], axis=0)
	std_son = np.std(sea_rcm[2:80:4], axis=0)
	std_ann = np.std(value[0:240:12,:,:], axis=0)
	
	djf_rcm = np.nanmean(sea_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(sea_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(sea_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(sea_rcm[2:80:4], axis=0)
	ann_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_ann, djf_rcm, mam_rcm, jja_rcm, son_rcm, ann_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	sea_gcm = value[2:240:3,:,:]

	std_djf = np.std(sea_gcm[3:80:4], axis=0)
	std_mam = np.std(sea_gcm[0:80:4], axis=0)
	std_jja = np.std(sea_gcm[1:80:4], axis=0)
	std_son = np.std(sea_gcm[2:80:4], axis=0)
	std_ann = np.std(value[0:240:12,:,:], axis=0)
	
	djf_gcm = np.nanmean(sea_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(sea_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(sea_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(sea_gcm[2:80:4], axis=0)
	ann_gcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf, std_mam, std_jja, std_son, std_ann, djf_gcm, mam_gcm, jja_gcm, son_gcm, ann_gcm
	

def ttest(mean_sample1, mean_sample2, std_sample1, std_sample2):

	# Calculate t statistics	
	p1 = mean_sample1 - mean_sample2 
	p2 = (std_sample1 - std_sample2) / np.sqrt(240)
	p3 = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(p3, df=240)

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
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)

	return map, xx, yy
	

# Import models and obs database 
lat, lon, pre_std_djf_cru, pre_std_mam_cru, pre_std_jja_cru, pre_std_son_cru, pre_std_ann_cru, pre_djf_cru, pre_mam_cru, pre_jja_cru, pre_son_cru, pre_ann_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
lat, lon, pre_std_djf_rcm, pre_std_mam_rcm, pre_std_jja_rcm, pre_std_son_rcm, pre_std_ann_rcm, pre_djf_rcm, pre_mam_rcm, pre_jja_rcm, pre_son_rcm, pre_ann_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_gcm, pre_std_mam_gcm, pre_std_jja_gcm, pre_std_son_gcm, pre_std_ann_gcm, pre_djf_gcm, pre_mam_gcm, pre_jja_gcm, pre_son_gcm, pre_ann_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, tas_std_djf_cru, tas_std_mam_cru, tas_std_jja_cru, tas_std_son_cru, tas_std_ann_cru, tas_djf_cru, tas_mam_cru, tas_jja_cru, tas_son_cru, tas_ann_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_std_djf_rcm, tas_std_mam_rcm, tas_std_jja_rcm, tas_std_son_rcm, tas_std_ann_rcm, tas_djf_rcm, tas_mam_rcm, tas_jja_rcm, tas_son_rcm, tas_ann_rcm = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_gcm, tas_std_mam_gcm, tas_std_jja_gcm, tas_std_son_gcm, tas_std_ann_gcm, tas_djf_gcm, tas_mam_gcm, tas_jja_gcm, tas_son_gcm, tas_ann_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute bias from models and obs database 
pre_djf_rcm_bias = pre_djf_rcm - pre_djf_cru
pre_mam_rcm_bias = pre_mam_rcm - pre_mam_cru
pre_jja_rcm_bias = pre_jja_rcm - pre_jja_cru
pre_son_rcm_bias = pre_son_rcm - pre_son_cru
pre_ann_rcm_bias = pre_ann_rcm - pre_ann_cru

pre_djf_gcm_bias = pre_djf_gcm - pre_djf_cru
pre_mam_gcm_bias = pre_mam_gcm - pre_mam_cru
pre_jja_gcm_bias = pre_jja_gcm - pre_jja_cru
pre_son_gcm_bias = pre_son_gcm - pre_son_cru
pre_ann_gcm_bias = pre_ann_gcm - pre_ann_cru

tas_djf_rcm_bias = np.nanmean(tas_djf_rcm, axis=0) - tas_djf_cru
tas_mam_rcm_bias = np.nanmean(tas_mam_rcm, axis=0) - tas_mam_cru
tas_jja_rcm_bias = np.nanmean(tas_jja_rcm, axis=0) - tas_jja_cru
tas_son_rcm_bias = np.nanmean(tas_son_rcm, axis=0) - tas_son_cru
tas_ann_rcm_bias = np.nanmean(tas_ann_rcm, axis=0) - tas_ann_cru

tas_djf_gcm_bias = tas_djf_gcm - tas_djf_cru
tas_mam_gcm_bias = tas_mam_gcm - tas_mam_cru
tas_jja_gcm_bias = tas_jja_gcm - tas_jja_cru
tas_son_gcm_bias = tas_son_gcm - tas_son_cru
tas_ann_gcm_bias = tas_ann_gcm - tas_ann_cru

# Compute added value from models and obs database 
av_pre_djf = compute_av(pre_djf_gcm, pre_djf_rcm, pre_djf_cru)
av_pre_mam = compute_av(pre_mam_gcm, pre_mam_rcm, pre_mam_cru)
av_pre_jja = compute_av(pre_jja_gcm, pre_jja_rcm, pre_jja_cru)
av_pre_son = compute_av(pre_son_gcm, pre_son_rcm, pre_son_cru)
av_pre_ann = compute_av(pre_ann_gcm, pre_ann_rcm, pre_ann_cru)

av_tas_djf = compute_av(tas_djf_gcm, np.nanmean(tas_djf_rcm, axis=0), tas_djf_cru)
av_tas_mam = compute_av(tas_mam_gcm, np.nanmean(tas_mam_rcm, axis=0), tas_mam_cru)
av_tas_jja = compute_av(tas_jja_gcm, np.nanmean(tas_jja_rcm, axis=0), tas_jja_cru)
av_tas_son = compute_av(tas_son_gcm, np.nanmean(tas_son_rcm, axis=0), tas_son_cru)
av_tas_ann = compute_av(tas_ann_gcm, np.nanmean(tas_ann_rcm, axis=0), tas_ann_cru)

# Compute ttest from models and obs database 
p_value_pre_djf_rcm_cru = ttest(pre_djf_rcm, pre_djf_cru, pre_std_djf_rcm, pre_std_djf_cru)
p_value_pre_mam_rcm_cru = ttest(pre_mam_rcm, pre_mam_cru, pre_std_mam_rcm, pre_std_mam_cru)
p_value_pre_jja_rcm_cru = ttest(pre_jja_rcm, pre_jja_cru, pre_std_jja_rcm, pre_std_jja_cru)
p_value_pre_son_rcm_cru = ttest(pre_son_rcm, pre_son_cru, pre_std_son_rcm, pre_std_son_cru)
p_value_pre_ann_rcm_cru = ttest(pre_ann_rcm, pre_ann_cru, pre_std_ann_rcm, pre_std_ann_cru)

p_value_pre_djf_gcm_cru = ttest(pre_djf_gcm, pre_djf_cru, pre_std_djf_gcm, pre_std_djf_cru)
p_value_pre_mam_gcm_cru = ttest(pre_mam_gcm, pre_mam_cru, pre_std_mam_gcm, pre_std_mam_cru)
p_value_pre_jja_gcm_cru = ttest(pre_jja_gcm, pre_jja_cru, pre_std_jja_gcm, pre_std_jja_cru)
p_value_pre_son_gcm_cru = ttest(pre_son_gcm, pre_son_cru, pre_std_son_gcm, pre_std_son_cru)
p_value_pre_ann_gcm_cru = ttest(pre_ann_gcm, pre_ann_cru, pre_std_ann_gcm, pre_std_ann_cru)

p_value_tas_djf_rcm_cru = ttest(tas_djf_rcm, tas_djf_cru, tas_std_djf_rcm, tas_std_djf_cru)
p_value_tas_mam_rcm_cru = ttest(tas_mam_rcm, tas_mam_cru, tas_std_mam_rcm, tas_std_mam_cru)
p_value_tas_jja_rcm_cru = ttest(tas_jja_rcm, tas_jja_cru, tas_std_jja_rcm, tas_std_jja_cru)
p_value_tas_son_rcm_cru = ttest(tas_son_rcm, tas_son_cru, tas_std_son_rcm, tas_std_son_cru)
p_value_tas_ann_rcm_cru = ttest(tas_ann_rcm, tas_ann_cru, tas_std_ann_rcm, tas_std_ann_cru)

p_value_tas_djf_gcm_cru = ttest(tas_djf_gcm, tas_djf_cru, tas_std_djf_gcm, tas_std_djf_cru)
p_value_tas_mam_gcm_cru = ttest(tas_mam_gcm, tas_mam_cru, tas_std_mam_gcm, tas_std_mam_cru)
p_value_tas_jja_gcm_cru = ttest(tas_jja_gcm, tas_jja_cru, tas_std_jja_gcm, tas_std_jja_cru)
p_value_tas_son_gcm_cru = ttest(tas_son_gcm, tas_son_cru, tas_std_son_gcm, tas_std_son_cru)
p_value_tas_ann_gcm_cru = ttest(tas_ann_gcm, tas_ann_cru, tas_std_ann_gcm, tas_std_ann_cru)

# Plot models and obs database 
fig = plt.figure(figsize=(7,7))
levs1 = [-6, -4, -2, 2, 4, 6]
levs2 = [-1, -0.5, -0.1, 0.1, 0.5, 1]
	
#~ ax = fig.add_subplot(5, 3, 1)
#~ plt.title(u'A) MBE RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_rcm_cru = ma.masked_where(p_value_pre_djf_rcm_cru >= 0.05, p_value_pre_djf_rcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_rcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 2)
#~ plt.title(u'B) MBE HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)  
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_gcm_cru = ma.masked_where(p_value_pre_djf_gcm_cru >= 0.05, p_value_pre_djf_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 3)
#~ plt.title(u'C) AV DJF', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_djf, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)   
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 4)
#~ plt.title(u'D) MBE RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_mam_rcm_cru = ma.masked_where(p_value_pre_mam_rcm_cru >= 0.05, p_value_pre_mam_rcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_mam_rcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 5)
#~ plt.title(u'E) MBE HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)  
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_mam_gcm_cru = ma.masked_where(p_value_pre_mam_gcm_cru >= 0.05, p_value_pre_mam_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_mam_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 6)
#~ plt.title(u'F) AV MAM', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_mam, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)  
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 7)
#~ plt.title(u'G) MBE RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_rcm_cru = ma.masked_where(p_value_pre_jja_rcm_cru >= 0.05, p_value_pre_jja_rcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_rcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 8)
#~ plt.title(u'H) MBE HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_gcm_cru = ma.masked_where(p_value_pre_jja_gcm_cru >= 0.05, p_value_pre_jja_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 9)
#~ plt.title(u'I) AV JJA', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_jja, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 10)
#~ plt.title(u'J) MBE RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_son_rcm_cru = ma.masked_where(p_value_pre_son_rcm_cru >= 0.05, p_value_pre_son_rcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_son_rcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 11)
#~ plt.title(u'K) MBE HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_son_gcm_cru = ma.masked_where(p_value_pre_son_gcm_cru >= 0.05, p_value_pre_son_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_son_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 12)
#~ plt.title(u'L) AV SON', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_son, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 13)
#~ plt.title(u'M) MBE RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_ann_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_rcm_cru = ma.masked_where(p_value_pre_ann_rcm_cru >= 0.05, p_value_pre_ann_rcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_rcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 14)
#~ plt.title(u'N) MBE HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_ann_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_gcm_cru = ma.masked_where(p_value_pre_ann_gcm_cru >= 0.05, p_value_pre_ann_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 15)
#~ plt.title(u'O) AV ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_ann, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)   
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_bias_pre_reg_had_obs_1986-2005.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()

ax = fig.add_subplot(5, 3, 1)
plt.title(u'A) MBE RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_cru = ma.masked_where(p_value_tas_djf_rcm_cru >= 0.05, p_value_tas_djf_rcm_cru) 
map.contourf(xx, yy, p_value_tas_djf_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 2)
plt.title(u'B) MBE HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm_cru = ma.masked_where(p_value_tas_djf_gcm_cru >= 0.05, p_value_tas_djf_gcm_cru) 
map.contourf(xx, yy, p_value_tas_djf_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 3)
plt.title(u'C) AV DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_djf, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 4)
plt.title(u'D) MBE RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_rcm_cru = ma.masked_where(p_value_tas_mam_rcm_cru >= 0.05, p_value_tas_mam_rcm_cru) 
map.contourf(xx, yy, p_value_tas_mam_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 5)
plt.title(u'E) MBE HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_gcm_cru = ma.masked_where(p_value_tas_mam_gcm_cru >= 0.05, p_value_tas_mam_gcm_cru) 
map.contourf(xx, yy, p_value_tas_mam_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 6)
plt.title(u'F) AV MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_mam, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 7)
plt.title(u'G) MBE RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_cru = ma.masked_where(p_value_tas_jja_rcm_cru >= 0.05, p_value_tas_jja_rcm_cru) 
map.contourf(xx, yy, p_value_tas_jja_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 8)
plt.title(u'H) MBE HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm_cru = ma.masked_where(p_value_tas_jja_gcm_cru >= 0.05, p_value_tas_jja_gcm_cru) 
map.contourf(xx, yy, p_value_tas_jja_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 9)
plt.title(u'I) AV JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_jja, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 10)
plt.title(u'J) MBE RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_rcm_cru = ma.masked_where(p_value_tas_son_rcm_cru >= 0.05, p_value_tas_son_rcm_cru) 
map.contourf(xx, yy, p_value_tas_son_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 11)
plt.title(u'K) MBE HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_gcm_cru = ma.masked_where(p_value_tas_son_gcm_cru >= 0.05, p_value_tas_son_gcm_cru) 
map.contourf(xx, yy, p_value_tas_son_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 12)
plt.title(u'L) AV SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_son, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 13)
plt.title(u'M) MBE RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_ann_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_rcm_cru = ma.masked_where(p_value_tas_ann_rcm_cru >= 0.05, p_value_tas_ann_rcm_cru) 
map.contourf(xx, yy, p_value_tas_ann_rcm_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 14)
plt.title(u'N) MBE HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_ann_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_gcm_cru = ma.masked_where(p_value_tas_ann_gcm_cru >= 0.05, p_value_tas_ann_gcm_cru) 
map.contourf(xx, yy, p_value_tas_ann_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 15)
plt.title(u'O) AV ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_ann, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



