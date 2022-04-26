# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal climatology maps from regcm47 and hadgem models and obs database"

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


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_obs = np.nanmean(value, axis=0)
	std_obs = np.std(value, axis=0)

	sea_obs = value[2:240:3,:,:]
	djf_obs = np.nanmean(sea_obs[3:80:4], axis=0)
	mam_obs = np.nanmean(sea_obs[0:80:4], axis=0)
	jja_obs = np.nanmean(sea_obs[1:80:4], axis=0)
	son_obs = np.nanmean(sea_obs[2:80:4], axis=0)
	ann_obs = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, mean_obs, std_obs, djf_obs, mam_obs, jja_obs, son_obs, ann_obs
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_rcm = np.nanmean(value, axis=0)
	std_rcm = np.std(value, axis=0)
	
	std = np.std(value, axis=0)
	sea_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(sea_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(sea_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(sea_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(sea_rcm[2:80:4], axis=0)
	ann_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, mean_rcm, std_rcm, djf_rcm, mam_rcm, jja_rcm, son_rcm, ann_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_gcm = np.nanmean(value, axis=0)
	std_gcm = np.std(value, axis=0)
	
	sea_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(sea_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(sea_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(sea_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(sea_gcm[2:80:4], axis=0)
	ann_gcm = np.nanmean(value[0:240:12,:,:], axis=0)
	
	return lat, lon, mean_gcm, std_gcm, djf_gcm, mam_gcm, jja_gcm, son_gcm, ann_gcm

	
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


def ttest(mean, std, sample):

	# Calculate t statistics
	p1 = mean - sample
	p2= std / np.sqrt(240)
	p3 = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(p3, df=240)
	
	return p_value


# Import models and obs database 
lat, lon, pre_cru_mean, pre_cru_std, pre_djf_cru, pre_mam_cru, pre_jja_cru, pre_son_cru, pre_ann_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
lat, lon, pre_rcm_mean, pre_rcm_std, pre_djf_rcm, pre_mam_rcm, pre_jja_rcm, pre_son_rcm, pre_ann_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_gcm_mean, pre_gcm_std, pre_djf_gcm, pre_mam_gcm, pre_jja_gcm, pre_son_gcm, pre_ann_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, tas_cru_mean, tas_cru_std, tas_djf_cru, tas_mam_cru, tas_jja_cru, tas_son_cru, tas_ann_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_rcm_mean, tas_rcm_std, tas_djf_rcm, tas_mam_rcm, tas_jja_rcm, tas_son_rcm, tas_ann_rcm = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_gcm_mean, tas_gcm_std, tas_djf_gcm, tas_mam_gcm, tas_jja_gcm, tas_son_gcm, tas_ann_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute ttest from models and obs database 
p_value_pre_djf_cru = ttest(pre_cru_mean, pre_cru_std, pre_djf_cru)
p_value_pre_mam_cru = ttest(pre_cru_mean, pre_cru_std, pre_mam_cru)
p_value_pre_jja_cru = ttest(pre_cru_mean, pre_cru_std, pre_jja_cru)
p_value_pre_son_cru = ttest(pre_cru_mean, pre_cru_std, pre_son_cru)
p_value_pre_ann_cru = ttest(pre_cru_mean, pre_cru_std, pre_ann_cru)

p_value_pre_djf_rcm = ttest(pre_rcm_mean, pre_rcm_std, pre_djf_rcm)
p_value_pre_mam_rcm = ttest(pre_rcm_mean, pre_rcm_std, pre_mam_rcm)
p_value_pre_jja_rcm = ttest(pre_rcm_mean, pre_rcm_std, pre_jja_rcm)
p_value_pre_son_rcm = ttest(pre_rcm_mean, pre_rcm_std, pre_son_rcm)
p_value_pre_ann_rcm = ttest(pre_rcm_mean, pre_rcm_std, pre_ann_rcm)

p_value_pre_djf_gcm = ttest(pre_gcm_mean, pre_gcm_std, pre_djf_gcm)
p_value_pre_mam_gcm = ttest(pre_gcm_mean, pre_gcm_std, pre_mam_gcm)
p_value_pre_jja_gcm = ttest(pre_gcm_mean, pre_gcm_std, pre_jja_gcm)
p_value_pre_son_gcm = ttest(pre_gcm_mean, pre_gcm_std, pre_son_gcm)
p_value_pre_ann_gcm = ttest(pre_gcm_mean, pre_gcm_std, pre_ann_gcm)

p_value_tas_djf_cru = ttest(tas_cru_mean, tas_cru_std, tas_djf_cru)
p_value_tas_mam_cru = ttest(tas_cru_mean, tas_cru_std, tas_mam_cru)
p_value_tas_jja_cru = ttest(tas_cru_mean, tas_cru_std, tas_jja_cru)
p_value_tas_son_cru = ttest(tas_cru_mean, tas_cru_std, tas_son_cru)
p_value_tas_ann_cru = ttest(tas_cru_mean, tas_cru_std, tas_ann_cru)

p_value_tas_djf_rcm = ttest(tas_rcm_mean, tas_rcm_std, tas_djf_rcm)
p_value_tas_mam_rcm = ttest(tas_rcm_mean, tas_rcm_std, tas_mam_rcm)
p_value_tas_jja_rcm = ttest(tas_rcm_mean, tas_rcm_std, tas_jja_rcm)
p_value_tas_son_rcm = ttest(tas_rcm_mean, tas_rcm_std, tas_son_rcm)
p_value_tas_ann_rcm = ttest(tas_rcm_mean, tas_rcm_std, tas_ann_rcm)

p_value_tas_djf_gcm = ttest(tas_gcm_mean, tas_gcm_std, tas_djf_gcm)
p_value_tas_mam_gcm = ttest(tas_gcm_mean, tas_gcm_std, tas_mam_gcm)
p_value_tas_jja_gcm = ttest(tas_gcm_mean, tas_gcm_std, tas_jja_gcm)
p_value_tas_son_gcm = ttest(tas_gcm_mean, tas_gcm_std, tas_son_gcm)
p_value_tas_ann_gcm = ttest(tas_gcm_mean, tas_gcm_std, tas_ann_gcm)

# Plot models and obs database 
fig = plt.figure(figsize=(7,7))
levs1 = [0, 3, 6, 9, 12, 15]
levs2 = [19, 21, 24, 27, 30, 33]

#~ ax = fig.add_subplot(5, 3, 1)
#~ plt.title(u'A) CRU DJF', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_djf_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
#~ map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 2)
#~ plt.title(u'B) RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_djf_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_rcm = ma.masked_where(p_value_pre_djf_rcm >= 0.05, p_value_pre_djf_rcm) 
#~ map.contourf(xx, yy, p_value_pre_djf_rcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 3)
#~ plt.title(u'C) HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_djf_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ cbar.set_label('Precipitação \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_djf_gcm = ma.masked_where(p_value_pre_djf_gcm >= 0.05, p_value_pre_djf_gcm) 
#~ map.contourf(xx, yy, p_value_pre_djf_gcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 4)
#~ plt.title(u'D) CRU MAM', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_mam_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_mam_cru = ma.masked_where(p_value_pre_mam_cru >= 0.05, p_value_pre_mam_cru) 
#~ map.contourf(xx, yy, p_value_pre_mam_cru, colors='none', hatches=['....'])
	
#~ ax = fig.add_subplot(5, 3, 5)
#~ plt.title(u'E) RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_mam_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)  
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_mam_rcm = ma.masked_where(p_value_pre_mam_rcm >= 0.05, p_value_pre_mam_rcm) 
#~ map.contourf(xx, yy, p_value_pre_mam_rcm, colors='none', hatches=['....'])
	
#~ ax = fig.add_subplot(5, 3, 6)
#~ plt.title(u'F) HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_mam_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ cbar.set_label('Precipitação \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_mam_gcm = ma.masked_where(p_value_pre_mam_gcm >= 0.05, p_value_pre_mam_gcm) 
#~ map.contourf(xx, yy, p_value_pre_mam_gcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 7)
#~ plt.title(u'G) CRU JJA', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_jja_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_cru = ma.masked_where(p_value_pre_jja_cru >= 0.05, p_value_pre_jja_cru) 
#~ map.contourf(xx, yy, p_value_pre_jja_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 8)
#~ plt.title(u'H) RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_jja_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
#~ map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_rcm = ma.masked_where(p_value_pre_jja_rcm >= 0.05, p_value_pre_jja_rcm) 
#~ map.contourf(xx, yy, p_value_pre_jja_rcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 9)
#~ plt.title(u'I) HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_jja_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max') 
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ cbar.set_label('Precipitação \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_jja_gcm = ma.masked_where(p_value_pre_jja_gcm >= 0.05, p_value_pre_jja_gcm) 
#~ map.contourf(xx, yy, p_value_pre_jja_gcm, colors='none', hatches=['....'])
	
#~ ax = fig.add_subplot(5, 3, 10)
#~ plt.title(u'J) CRU SON', loc='left', fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_son_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_son_cru = ma.masked_where(p_value_pre_son_cru >= 0.05, p_value_pre_son_cru) 
#~ map.contourf(xx, yy, p_value_pre_son_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 11)
#~ plt.title(u'K) RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_son_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)
#~ map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_son_rcm = ma.masked_where(p_value_pre_son_rcm >= 0.05, p_value_pre_son_rcm) 
#~ map.contourf(xx, yy, p_value_pre_son_rcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 12)
#~ plt.title(u'L) HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_son_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ cbar.set_label('Precipitação \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_son_gcm = ma.masked_where(p_value_pre_son_gcm >= 0.05, p_value_pre_son_gcm) 
#~ map.contourf(xx, yy, p_value_pre_son_gcm, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 13)
#~ plt.title(u'M) CRU ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_ann_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_cru = ma.masked_where(p_value_pre_ann_cru >= 0.05, p_value_pre_ann_cru) 
#~ map.contourf(xx, yy, p_value_pre_ann_cru, colors='none', hatches=['....'])

#~ ax = fig.add_subplot(5, 3, 14)
#~ plt.title(u'N) RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_ann_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_rcm = ma.masked_where(p_value_pre_ann_rcm >= 0.05, p_value_pre_ann_rcm) 
#~ map.contourf(xx, yy, p_value_pre_ann_rcm, colors='none', hatches=['....'])
	
#~ ax = fig.add_subplot(5, 3, 15)
#~ plt.title(u'O) HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
#~ map, xx, yy = basemap(lat, lon)
#~ plot_maps_mean = map.contourf(xx, yy, pre_ann_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')  
#~ cbar = map.colorbar()
#~ cbar.ax.tick_params(labelsize=6) 
#~ cbar.set_label('Precipitação \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_ann_gcm = ma.masked_where(p_value_pre_ann_gcm >= 0.05, p_value_pre_ann_gcm) 
#~ map.contourf(xx, yy, p_value_pre_ann_gcm, colors='none', hatches=['....'])

#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_clim_pre_reg_had_obs_1986-2005.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
#~ plt.show()
#~ exit()

ax = fig.add_subplot(5, 3, 1)
plt.title(u'A) CRU DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_cru = ma.masked_where(p_value_tas_djf_cru >= 0.05, p_value_tas_djf_cru) 
map.contourf(xx, yy, p_value_tas_djf_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 2)
plt.title(u'B) RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm = ma.masked_where(p_value_tas_djf_rcm >= 0.05, p_value_tas_djf_rcm) 
map.contourf(xx, yy, p_value_tas_djf_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 3)
plt.title(u'C) HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperatura \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm = ma.masked_where(p_value_tas_djf_gcm >= 0.05, p_value_tas_djf_gcm) 
map.contourf(xx, yy, p_value_tas_djf_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 4)
plt.title(u'D) CRU MAM', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_cru = ma.masked_where(p_value_tas_mam_cru >= 0.05, p_value_tas_mam_cru) 
map.contourf(xx, yy, p_value_tas_mam_cru, colors='none', hatches=['....'])
	
ax = fig.add_subplot(5, 3, 5)
plt.title(u'E) RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_rcm = ma.masked_where(p_value_tas_mam_rcm >= 0.05, p_value_tas_mam_rcm) 
map.contourf(xx, yy, p_value_tas_mam_rcm[0,:,:], colors='none', hatches=['....'])
	
ax = fig.add_subplot(5, 3, 6)
plt.title(u'F) HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperatura \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_gcm = ma.masked_where(p_value_tas_mam_gcm >= 0.05, p_value_tas_mam_gcm) 
map.contourf(xx, yy, p_value_tas_mam_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 7)
plt.title(u'G) CRU JJA', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_cru = ma.masked_where(p_value_tas_jja_cru >= 0.05, p_value_tas_jja_cru) 
map.contourf(xx, yy, p_value_tas_jja_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 8)
plt.title(u'H) RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm = ma.masked_where(p_value_tas_jja_rcm >= 0.05, p_value_tas_jja_rcm) 
map.contourf(xx, yy, p_value_tas_jja_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 9)
plt.title(u'I) HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperatura \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm = ma.masked_where(p_value_tas_jja_gcm >= 0.05, p_value_tas_jja_gcm) 
map.contourf(xx, yy, p_value_tas_jja_gcm, colors='none', hatches=['....'])
	
ax = fig.add_subplot(5, 3, 10)
plt.title(u'J) CRU SON', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_cru = ma.masked_where(p_value_tas_son_cru >= 0.05, p_value_tas_son_cru) 
map.contourf(xx, yy, p_value_tas_son_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 11)
plt.title(u'K) RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_rcm = ma.masked_where(p_value_tas_son_rcm >= 0.05, p_value_tas_son_rcm) 
map.contourf(xx, yy, p_value_tas_son_rcm[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 12)
plt.title(u'L) HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperatura \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_gcm = ma.masked_where(p_value_tas_son_gcm >= 0.05, p_value_tas_son_gcm) 
map.contourf(xx, yy, p_value_tas_son_gcm, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 13)
plt.title(u'M) CRU ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_cru = ma.masked_where(p_value_tas_ann_cru >= 0.05, p_value_tas_ann_cru) 
map.contourf(xx, yy, p_value_tas_ann_cru, colors='none', hatches=['....'])

ax = fig.add_subplot(5, 3, 14)
plt.title(u'N) RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_rcm = ma.masked_where(p_value_tas_ann_rcm >= 0.05, p_value_tas_ann_rcm) 
map.contourf(xx, yy, p_value_tas_ann_rcm[0,:,:], colors='none', hatches=['....'])
	
ax = fig.add_subplot(5, 3, 15)
plt.title(u'O) HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')  
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperatura \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_ann_gcm = ma.masked_where(p_value_tas_ann_gcm >= 0.05, p_value_tas_ann_gcm) 
map.contourf(xx, yy, p_value_tas_ann_gcm, colors='none', hatches=['....'])

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()





