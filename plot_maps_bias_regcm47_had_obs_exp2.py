# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
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
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	season_obs = value[2:240:3,:,:]

	std_djf_obs = np.std(season_obs[3:80:4], axis=0)
	std_jja_obs = np.std(season_obs[1:80:4], axis=0)
	std_ann_obs = np.std(value[0:240:12,:,:], axis=0)
		
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:80:4], axis=0)
	ann_obs = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf_obs, std_jja_obs, std_ann_obs, djf_obs, jja_obs, ann_obs
	
	
def import_rcm_exp1(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	season_rcm = value[2:240:3,:,:]

	std_djf_rcm = np.std(season_rcm[3:80:4], axis=0)
	std_jja_rcm = np.std(season_rcm[1:80:4], axis=0)
	std_ann_rcm = np.std(value[0:240:12,:,:], axis=0)
	
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	ann_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf_rcm, std_jja_rcm, std_ann_rcm, djf_rcm, jja_rcm, ann_rcm


def import_rcm_exp2(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	season_rcm = value[2:240:3,:,:]

	std_djf_rcm = np.std(season_rcm[3:80:4], axis=0)
	std_jja_rcm = np.std(season_rcm[1:80:4], axis=0)
	std_ann_rcm = np.std(value[0:240:12,:,:], axis=0)
	
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	ann_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf_rcm, std_jja_rcm, std_ann_rcm, djf_rcm, jja_rcm, ann_rcm
	
	
def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	season_gcm = value[2:240:3,:,:]

	std_djf_gcm = np.std(season_gcm[3:80:4], axis=0)
	std_jja_gcm = np.std(season_gcm[1:80:4], axis=0)
	std_ann_gcm = np.std(value[0:240:12,:,:], axis=0)
	
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)
	ann_gcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, std_djf_gcm, std_jja_gcm, std_ann_gcm, djf_gcm, jja_gcm, ann_gcm
	

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
lat, lon, pre_std_djf_cru, pre_std_jja_cru, pre_std_ann_cru, pre_djf_cru, pre_jja_cru, pre_ann_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
lat, lon, pre_std_djf_rcm_exp1, pre_std_jja_rcm_exp1, pre_std_ann_rcm_exp1, pre_djf_rcm_exp1, pre_jja_rcm_exp1, pre_ann_rcm_exp1 = import_rcm_exp1('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_rcm_exp2, pre_std_jja_rcm_exp2, pre_std_ann_rcm_exp2, pre_djf_rcm_exp2, pre_jja_rcm_exp2, pre_ann_rcm_exp2 = import_rcm_exp2('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_gcm, pre_std_jja_gcm, pre_std_ann_gcm, pre_djf_gcm, pre_jja_gcm, pre_ann_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, tas_std_djf_cru, tas_std_jja_cru, tas_std_ann_cru, tas_djf_cru, tas_jja_cru, tas_ann_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_std_djf_rcm_exp1, tas_std_jja_rcm_exp1, tas_std_ann_rcm_exp1, tas_djf_rcm_exp1, tas_jja_rcm_exp1, tas_ann_rcm_exp1 = import_rcm_exp1('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_rcm_exp2, tas_std_jja_rcm_exp2, tas_std_ann_rcm_exp2, tas_djf_rcm_exp2, tas_jja_rcm_exp2, tas_ann_rcm_exp2 = import_rcm_exp2('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_gcm, tas_std_jja_gcm, tas_std_ann_gcm, tas_djf_gcm, tas_jja_gcm, tas_ann_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute bias from models and obs database 
pre_djf_rcm_exp1_bias = pre_djf_rcm_exp1 - pre_djf_cru
pre_jja_rcm_exp1_bias = pre_jja_rcm_exp1 - pre_jja_cru
pre_ann_rcm_exp1_bias = pre_ann_rcm_exp1 - pre_ann_cru
pre_djf_rcm_exp2_bias = pre_djf_rcm_exp2 - pre_djf_cru
pre_jja_rcm_exp2_bias = pre_jja_rcm_exp2 - pre_jja_cru
pre_ann_rcm_exp2_bias = pre_ann_rcm_exp2 - pre_ann_cru
pre_djf_gcm_bias = pre_djf_gcm - pre_djf_cru
pre_jja_gcm_bias = pre_jja_gcm - pre_jja_cru
pre_ann_gcm_bias = pre_ann_gcm - pre_ann_cru

tas_djf_rcm_exp1_bias = np.nanmean(tas_djf_rcm_exp1, axis=0) - tas_djf_cru
tas_jja_rcm_exp1_bias = np.nanmean(tas_jja_rcm_exp1, axis=0) - tas_jja_cru
tas_ann_rcm_exp1_bias = np.nanmean(tas_ann_rcm_exp1, axis=0) - tas_ann_cru
tas_djf_rcm_exp2_bias = np.nanmean(tas_djf_rcm_exp2, axis=0) - tas_djf_cru
tas_jja_rcm_exp2_bias = np.nanmean(tas_jja_rcm_exp2, axis=0) - tas_jja_cru
tas_ann_rcm_exp2_bias = np.nanmean(tas_ann_rcm_exp2, axis=0) - tas_ann_cru
tas_djf_gcm_bias = tas_djf_gcm - tas_djf_cru
tas_jja_gcm_bias = tas_jja_gcm - tas_jja_cru
tas_ann_gcm_bias = tas_ann_gcm - tas_ann_cru

# Compute ttest from models and obs database 
p_value_pre_djf_rcm_exp1_cru = ttest(pre_djf_rcm_exp1, pre_djf_cru, pre_std_djf_rcm_exp1, pre_std_djf_cru)
p_value_pre_jja_rcm_exp1_cru = ttest(pre_jja_rcm_exp1, pre_djf_cru, pre_std_jja_rcm_exp1, pre_std_jja_cru)
p_value_pre_ann_rcm_exp1_cru = ttest(pre_ann_rcm_exp1, pre_djf_cru, pre_std_ann_rcm_exp1, pre_std_ann_cru)
p_value_pre_djf_rcm_exp2_cru = ttest(pre_djf_rcm_exp2, pre_djf_cru, pre_std_djf_rcm_exp2, pre_std_djf_cru)
p_value_pre_jja_rcm_exp2_cru = ttest(pre_jja_rcm_exp2, pre_djf_cru, pre_std_jja_rcm_exp2, pre_std_jja_cru)
p_value_pre_ann_rcm_exp2_cru = ttest(pre_ann_rcm_exp2, pre_djf_cru, pre_std_ann_rcm_exp2, pre_std_ann_cru)
p_value_pre_djf_gcm_cru = ttest(pre_djf_gcm, pre_djf_cru, pre_std_djf_gcm, pre_std_djf_cru)
p_value_pre_jja_gcm_cru = ttest(pre_jja_gcm, pre_djf_cru, pre_std_jja_gcm, pre_std_jja_cru)
p_value_pre_ann_gcm_cru = ttest(pre_ann_gcm, pre_djf_cru, pre_std_ann_gcm, pre_std_ann_cru)

p_value_tas_djf_rcm_exp1_cru = ttest(tas_djf_rcm_exp1, pre_djf_cru, pre_std_djf_rcm_exp1, pre_std_djf_cru)
p_value_tas_jja_rcm_exp1_cru = ttest(tas_jja_rcm_exp1, pre_djf_cru, pre_std_jja_rcm_exp1, pre_std_jja_cru)
p_value_tas_ann_rcm_exp1_cru = ttest(tas_ann_rcm_exp1, pre_djf_cru, pre_std_ann_rcm_exp1, pre_std_ann_cru)
p_value_tas_djf_rcm_exp2_cru = ttest(tas_djf_rcm_exp2, pre_djf_cru, pre_std_djf_rcm_exp2, pre_std_djf_cru)
p_value_tas_jja_rcm_exp2_cru = ttest(tas_jja_rcm_exp2, pre_djf_cru, pre_std_jja_rcm_exp2, pre_std_jja_cru)
p_value_tas_ann_rcm_exp2_cru = ttest(tas_ann_rcm_exp2, pre_djf_cru, pre_std_ann_rcm_exp2, pre_std_ann_cru)
p_value_tas_djf_gcm_cru = ttest(tas_djf_gcm, pre_djf_cru, tas_std_djf_gcm, tas_std_djf_cru)
p_value_tas_jja_gcm_cru = ttest(tas_jja_gcm, pre_djf_cru, tas_std_jja_gcm, tas_std_jja_cru)
p_value_tas_ann_gcm_cru = ttest(tas_ann_gcm, pre_djf_cru, tas_std_ann_gcm, tas_std_ann_cru)

# Compute added from models and obs database 
av_pre_djf_exp1 = compute_av(pre_djf_gcm, pre_djf_rcm_exp1, pre_djf_cru)
av_pre_jja_exp1 = compute_av(pre_jja_gcm, pre_jja_rcm_exp1, pre_jja_cru)
av_pre_ann_exp1 = compute_av(pre_ann_gcm, pre_ann_rcm_exp1, pre_ann_cru)
av_pre_djf_exp2 = compute_av(pre_djf_gcm, pre_djf_rcm_exp2, pre_djf_cru)
av_pre_jja_exp2 = compute_av(pre_jja_gcm, pre_jja_rcm_exp2, pre_jja_cru)
av_pre_ann_exp2 = compute_av(pre_ann_gcm, pre_ann_rcm_exp2, pre_ann_cru)

av_tas_djf_exp1 = compute_av(tas_djf_gcm, np.nanmean(tas_djf_rcm_exp1, axis=0), tas_djf_cru)
av_tas_jja_exp1 = compute_av(tas_jja_gcm, np.nanmean(tas_jja_rcm_exp1, axis=0), tas_jja_cru)
av_tas_ann_exp1 = compute_av(tas_ann_gcm, np.nanmean(tas_ann_rcm_exp1, axis=0), tas_ann_cru)
av_tas_djf_exp2 = compute_av(tas_djf_gcm, np.nanmean(tas_djf_rcm_exp2, axis=0), tas_djf_cru)
av_tas_jja_exp2 = compute_av(tas_jja_gcm, np.nanmean(tas_jja_rcm_exp2, axis=0), tas_jja_cru)
av_tas_ann_exp2 = compute_av(tas_ann_gcm, np.nanmean(tas_ann_rcm_exp2, axis=0), tas_ann_cru)

# Plot models and obs database 
fig = plt.figure(figsize=(5, 5))
levs1 = [-6, -3, -1, 1, 3, 6]
levs2 = [-1, -0.7, -0.4, 0.4, 0.7, 1]

#~ ax = fig.add_subplot(5, 3, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_djf_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_rcm_exp1_cru = ma.masked_where(p_value_pre_djf_rcm_exp1_cru >= 0.1, p_value_pre_djf_rcm_exp1_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_rcm_exp1_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_jja_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_rcm_exp1_cru = ma.masked_where(p_value_pre_jja_rcm_exp1_cru >= 0.1, p_value_pre_jja_rcm_exp1_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_rcm_exp1_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_ann_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)    
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_rcm_exp1_cru = ma.masked_where(p_value_pre_ann_rcm_exp1_cru >= 0.1, p_value_pre_ann_rcm_exp1_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_rcm_exp1_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_djf_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_rcm_exp2_cru = ma.masked_where(p_value_pre_djf_rcm_exp2_cru >= 0.1, p_value_pre_djf_rcm_exp2_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_rcm_exp2_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_jja_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_rcm_exp2_cru = ma.masked_where(p_value_pre_jja_rcm_exp2_cru >= 0.1, p_value_pre_jja_rcm_exp2_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_rcm_exp2_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 6)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_ann_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)    
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_rcm_exp2_cru = ma.masked_where(p_value_pre_ann_rcm_exp2_cru >= 0.1, p_value_pre_ann_rcm_exp2_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_rcm_exp2_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 7)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_gcm_cru = ma.masked_where(p_value_pre_djf_gcm_cru >= 0.1, p_value_pre_djf_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 8)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
#~ plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_gcm_cru = ma.masked_where(p_value_pre_jja_gcm_cru >= 0.1, p_value_pre_jja_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 9)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, pre_ann_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
#~ plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)    
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_gcm_cru = ma.masked_where(p_value_pre_ann_gcm_cru >= 0.1, p_value_pre_ann_gcm_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_gcm_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 10)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_djf_exp1, levels=levs2, latlon=True, cmap=cm.PiYG)
#~ plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 11)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_jja_exp1, levels=levs2, latlon=True, cmap=cm.PiYG) 
#~ plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 12)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_ann_exp1, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
#~ plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)   
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 13)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_djf_exp2, levels=levs2, latlon=True, cmap=cm.PiYG)
#~ plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 14)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_jja_exp2, levels=levs2, latlon=True, cmap=cm.PiYG) 
#~ plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

#~ ax = fig.add_subplot(5, 3, 15)
#~ map, xx, yy = basemap(lat, lon)
#~ plt_maps_bias = map.contourf(xx, yy, av_pre_ann_exp2, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
#~ plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6)   
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 1)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_exp1_cru = ma.masked_where(p_value_tas_djf_rcm_exp1_cru >= 0.1, p_value_tas_djf_rcm_exp1_cru) 
map.contourf(xx, yy, p_value_tas_djf_rcm_exp1_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 2)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_exp1_cru = ma.masked_where(p_value_tas_jja_rcm_exp1_cru >= 0.1, p_value_tas_jja_rcm_exp1_cru) 
map.contourf(xx, yy, p_value_tas_jja_rcm_exp1_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 3)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_ann_rcm_exp1_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)    
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_rcm_exp1_cru = ma.masked_where(p_value_tas_ann_rcm_exp1_cru >= 0.1, p_value_tas_ann_rcm_exp1_cru) 
map.contourf(xx, yy, p_value_tas_ann_rcm_exp1_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 4)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_exp2_cru = ma.masked_where(p_value_tas_djf_rcm_exp2_cru >= 0.1, p_value_tas_djf_rcm_exp2_cru) 
map.contourf(xx, yy, p_value_tas_djf_rcm_exp2_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 5)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_exp2_cru = ma.masked_where(p_value_tas_jja_rcm_exp2_cru >= 0.1, p_value_tas_jja_rcm_exp2_cru) 
map.contourf(xx, yy, p_value_tas_jja_rcm_exp2_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 6)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_ann_rcm_exp2_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)    
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_rcm_exp2_cru = ma.masked_where(p_value_tas_ann_rcm_exp2_cru >= 0.1, p_value_tas_ann_rcm_exp2_cru) 
map.contourf(xx, yy, p_value_tas_ann_rcm_exp2_cru[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 7)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm_cru = ma.masked_where(p_value_tas_djf_gcm_cru >= 0.1, p_value_tas_djf_gcm_cru) 
map.contourf(xx, yy, p_value_tas_djf_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 8)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm_cru = ma.masked_where(p_value_tas_jja_gcm_cru >= 0.1, p_value_tas_jja_gcm_cru) 
map.contourf(xx, yy, p_value_tas_jja_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 9)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_ann_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr, extend='both')
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)    
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_gcm_cru = ma.masked_where(p_value_tas_ann_gcm_cru >= 0.1, p_value_tas_ann_gcm_cru) 
map.contourf(xx, yy, p_value_tas_ann_gcm_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 10)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_djf_exp1, levels=levs2, latlon=True, cmap=cm.PiYG)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 11)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_jja_exp1, levels=levs2, latlon=True, cmap=cm.PiYG) 
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 12)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_ann_exp1, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 13)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_djf_exp2, levels=levs2, latlon=True, cmap=cm.PiYG)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 14)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_jja_exp2, levels=levs2, latlon=True, cmap=cm.PiYG) 
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 15)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_ann_exp2, levels=levs2, latlon=True, cmap=cm.PiYG, extend='both')
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)   
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_reg_exp2_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()





