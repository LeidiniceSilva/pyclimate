# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal bias maps from Reg and Had models output"

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from comp_statist_indices import compute_added_value


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_obs = value[2:240:3,:,:]
	djf_obs = np.nanmean(season_obs[3:80:4], axis=0)
	mam_obs = np.nanmean(season_obs[0:80:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:80:4], axis=0)
	son_obs = np.nanmean(season_obs[2:80:4], axis=0)
	annual_obs = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, djf_obs, mam_obs, jja_obs, son_obs, annual_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)
	annual_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm, annual_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(season_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(season_gcm[2:80:4], axis=0)
	annual_gcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm, annual_gcm
	

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
	

def ttest(data1, data2):
	
	# Calculate the mean and standard error
	x1_bar, x2_bar = np.nanmean(data1, axis=0), np.nanmean(data2, axis=0)
	var_x1, var_x2= np.var(data1, ddof=1), np.var(data2, ddof=1)

	# Pooled sample variance
	n1, n2 = 240, 240
	pool_var = ( ((n1-1)*var_x1) + ((n2-1)*var_x2) ) / (n1+n2-2)

	# Standard error
	std_error = np.sqrt(pool_var * (1.0 / n1 + 1.0 / n2))

	# Calculate t statistics
	t = abs(x1_bar - x2_bar) / std_error

	# Calculate p value
	p_value = 1 - stats.t.cdf(x=t, df=238)
	
	return p_value
	

# Import regcm exp and cru databases 	   
lat, lon, pre_djf_cru, pre_mam_cru, pre_jja_cru, pre_son_cru, pre_annual_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
lat, lon, pre_djf_rcm, pre_mam_rcm, pre_jja_rcm, pre_son_rcm, pre_annual_rcm = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_djf_gcm, pre_mam_gcm, pre_jja_gcm, pre_son_gcm, pre_annual_gcm = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')

lat, lon, tas_djf_cru, tas_mam_cru, tas_jja_cru, tas_son_cru, tas_annual_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
lat, lon, tas_djf_rcm, tas_mam_rcm, tas_jja_rcm, tas_son_rcm, tas_annual_rcm = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_djf_gcm, tas_mam_gcm, tas_jja_gcm, tas_son_gcm, tas_annual_gcm = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')

# Compute and plot bias from regcm exp and cru database
pre_djf_rcm_bias = pre_djf_rcm - pre_djf_cru
pre_mam_rcm_bias = pre_mam_rcm - pre_mam_cru
pre_jja_rcm_bias = pre_jja_rcm - pre_jja_cru
pre_son_rcm_bias = pre_son_rcm - pre_son_cru

pre_djf_gcm_bias = pre_djf_gcm - pre_djf_cru
pre_mam_gcm_bias = pre_mam_gcm - pre_mam_cru
pre_jja_gcm_bias = pre_jja_gcm - pre_jja_cru
pre_son_gcm_bias = pre_son_gcm - pre_son_cru

tas_djf_rcm_bias = np.nanmean(tas_djf_rcm, axis=0) - tas_djf_cru
tas_mam_rcm_bias = np.nanmean(tas_mam_rcm, axis=0) - tas_mam_cru
tas_jja_rcm_bias = np.nanmean(tas_jja_rcm, axis=0) - tas_jja_cru
tas_son_rcm_bias = np.nanmean(tas_son_rcm, axis=0) - tas_son_cru

tas_djf_gcm_bias = tas_djf_gcm - tas_djf_cru
tas_mam_gcm_bias = tas_mam_gcm - tas_mam_cru
tas_jja_gcm_bias = tas_jja_gcm - tas_jja_cru
tas_son_gcm_bias = tas_son_gcm - tas_son_cru

# Compute added value from regcm output
av_pre_djf = compute_added_value(pre_djf_gcm, pre_djf_rcm, pre_djf_cru)
av_pre_mam = compute_added_value(pre_mam_gcm, pre_mam_rcm, pre_mam_cru)
av_pre_jja = compute_added_value(pre_jja_gcm, pre_jja_rcm, pre_jja_cru)
av_pre_son = compute_added_value(pre_son_gcm, pre_son_rcm, pre_son_cru)

av_tas_djf = compute_added_value(tas_djf_gcm, np.nanmean(tas_djf_rcm, axis=0), tas_djf_cru)
av_tas_mam = compute_added_value(tas_mam_gcm, np.nanmean(tas_mam_rcm, axis=0), tas_mam_cru)
av_tas_jja = compute_added_value(tas_jja_gcm, np.nanmean(tas_jja_rcm, axis=0), tas_jja_cru)
av_tas_son = compute_added_value(tas_son_gcm, np.nanmean(tas_son_rcm, axis=0), tas_son_cru)

p_value_pre_djf_rcm_cru = ttest(pre_djf_rcm, pre_djf_cru)
p_value_pre_mam_rcm_cru = ttest(pre_mam_rcm, pre_mam_cru)
p_value_pre_jja_rcm_cru = ttest(pre_jja_rcm, pre_jja_cru)
p_value_pre_son_rcm_cru = ttest(pre_son_rcm, pre_son_cru)
p_value_pre_annual_rcm_cru = ttest(pre_annual_rcm, pre_annual_cru)

# Compute ttest
#~ p_value_pre_djf_gcm_cru = ttest(pre_djf_gcm, pre_djf_cru)
#~ p_value_pre_mam_gcm_cru = ttest(pre_mam_gcm, pre_mam_cru)
#~ p_value_pre_jja_gcm_cru = ttest(pre_jja_gcm, pre_jja_cru)
#~ p_value_pre_son_gcm_cru = ttest(pre_son_gcm, pre_son_cru)
#~ p_value_pre_annual_gcm_cru = ttest(pre_annual_gcm, pre_annual_cru)

#~ p_value_pre_djf_rcm_gcm = ttest(pre_djf_rcm, pre_djf_gcm)
#~ p_value_pre_mam_rcm_gcm = ttest(pre_mam_rcm, pre_mam_gcm)
#~ p_value_pre_jja_rcm_gcm = ttest(pre_jja_rcm, pre_jja_gcm)
#~ p_value_pre_son_rcm_gcm = ttest(pre_son_rcm, pre_son_gcm)
#~ p_value_pre_annual_rcm_gcm = ttest(pre_annual_rcm, pre_annual_gcm)

#~ p_value_tas_djf_rcm_cru = ttest(tas_djf_rcm, tas_djf_cru)
#~ p_value_tas_mam_rcm_cru = ttest(tas_mam_rcm, tas_mam_cru)
#~ p_value_tas_jja_rcm_cru = ttest(tas_jja_rcm, tas_jja_cru)
#~ p_value_tas_son_rcm_cru = ttest(tas_son_rcm, tas_son_cru)
#~ p_value_tas_annual_rcm_cru = ttest(tas_annual_rcm, tas_annual_cru)

#~ p_value_tas_djf_gcm_cru = ttest(tas_djf_gcm, tas_djf_cru)
#~ p_value_tas_mam_gcm_cru = ttest(tas_mam_gcm, tas_mam_cru)
#~ p_value_tas_jja_gcm_cru = ttest(tas_jja_gcm, tas_jja_cru)
#~ p_value_tas_son_gcm_cru = ttest(tas_son_gcm, tas_son_cru)
#~ p_value_tas_annual_gcm_cru = ttest(tas_annual_gcm, tas_annual_cru)

#~ p_value_tas_djf_rcm_gcm = ttest(tas_djf_rcm, tas_djf_gcm)
#~ p_value_tas_mam_rcm_gcm = ttest(tas_mam_rcm, tas_mam_gcm)
#~ p_value_tas_jja_rcm_gcm = ttest(tas_jja_rcm, tas_jja_gcm)
#~ p_value_tas_son_rcm_gcm = ttest(tas_son_rcm, tas_son_gcm)
#~ p_value_tas_annual_rcm_gcm = ttest(tas_annual_rcm, tas_annual_gcm)

# Plot bias maps 
fig = plt.figure()

levs1 = [-6, -4, -2, -1, 1, 2, 4, 6]
levs2 = [-1, -0.7, -0.4, 0, 0.4, 0.7, 1]
	
ax = fig.add_subplot(6, 5, 1)
plt.title(u'A)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
plt_maps_bias = map.contourf(xx, yy, pre_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_bias = ma.masked_where(p_value_pre_djf_bias >= 0.05, p_value_pre_djf_bias) 
map.contourf(xx, yy, p_value_pre_djf_bias, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 2)
plt.title(u'B)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_bias = ma.masked_where(p_value_pre_mam_bias >= 0.05, p_value_pre_mam_bias) 
map.contourf(xx, yy, p_value_pre_mam_bias, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 3)
plt.title(u'C)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_bias = ma.masked_where(p_value_pre_jja_bias >= 0.05, p_value_pre_jja_bias) 
map.contourf(xx, yy, p_value_pre_jja_bias, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 4)
plt.title(u'D)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_bias = ma.masked_where(p_value_pre_son_bias >= 0.05, p_value_pre_son_bias) 
map.contourf(xx, yy, p_value_pre_son_bias, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 5)
plt.title(u'E)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_annual_rcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_bias = ma.masked_where(p_value_pre_annual_bias >= 0.05, p_value_pre_annual_bias) 
map.contourf(xx, yy, p_value_pre_annual_bias, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 6)
plt.title(u'F)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 7)
plt.title(u'G)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(6, 5, 8)
plt.title(u'H)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 9)
plt.title(u'I)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 10)
plt.title(u'J)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, pre_annual_gcm_bias, levels=levs1, latlon=True, cmap=cm.BrBG)
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 11)
plt.title(u'K)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_djf, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 12)
plt.title(u'L)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_mam, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 13)
plt.title(u'M)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_jja, levels=levs2, latlon=True, cmap=cm.PiYG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 14)
plt.title(u'N)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_son, levels=levs2, latlon=True, cmap=cm.PiYG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 15)
plt.title(u'O)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_pre_annual, levels=levs2, latlon=True, cmap=cm.PiYG)
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 16)
plt.title(u'P)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 17)
plt.title(u'Q)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 18)
plt.title(u'R)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 19)
plt.title(u'S)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_rcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 20)
plt.title(u'T)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_annual_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 21)
plt.title(u'U)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_djf_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 22)
plt.title(u'V)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_mam_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 23)
plt.title(u'W)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_jja_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 24)
plt.title(u'X)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_son_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 25)
plt.title(u'Y)', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, tas_annual_gcm_bias, levels=levs1, latlon=True, cmap=cm.bwr)
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_cru = ma.masked_where(p_value_pre_djf_cru >= 0.05, p_value_pre_djf_cru) 
map.contourf(xx, yy, p_value_pre_djf_cru, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 26)
plt.title(u'Z)', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_djf, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf = ma.masked_where(p_value_tas_djf >= 0.05, p_value_tas_djf) 
map.contourf(xx, yy, p_value_tas_djf, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 27)
plt.title(u'A.1)', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_mam, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam = ma.masked_where(p_value_tas_mam >= 0.05, p_value_tas_mam) 
map.contourf(xx, yy, p_value_tas_mam, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 28)
plt.title(u'B.1)', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_jja, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja = ma.masked_where(p_value_tas_jja >= 0.05, p_value_tas_jja) 
map.contourf(xx, yy, p_value_tas_jja, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 29)
plt.title(u'C.1)', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_son, levels=levs2, latlon=True, cmap=cm.PiYG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son = ma.masked_where(p_value_tas_son >= 0.05, p_value_tas_son) 
map.contourf(xx, yy, p_value_tas_son, colors='none', hatches=["////"])

ax = fig.add_subplot(6, 5, 30)
plt.title(u'D.1)', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
map, xx, yy = basemap(lat, lon)
plt_maps_bias = map.contourf(xx, yy, av_tas_annual, levels=levs2, latlon=True, cmap=cm.PiYG) 
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual = ma.masked_where(p_value_tas_annual >= 0.05, p_value_tas_annual) 
map.contourf(xx, yy, p_value_tas_annual, colors='none', hatches=["////"])

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_seasonal_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()



