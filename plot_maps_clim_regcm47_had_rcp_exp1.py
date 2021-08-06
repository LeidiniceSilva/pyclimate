# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal climatology maps from regcm47 and hadgem models to rcp"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib as mpl 
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt
from comp_statist_indices import compute_relative_change

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
    
    
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
	

# Import models
lat, lon, pre_std_djf_rcm_hist, pre_std_mam_rcm_hist, pre_std_jja_rcm_hist, pre_std_son_rcm_hist, pre_std_annual_rcm_hist, pre_djf_rcm_hist, pre_mam_rcm_hist, pre_jja_rcm_hist, pre_son_rcm_hist, pre_annual_rcm_hist = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_gcm_hist, pre_std_mam_gcm_hist, pre_std_jja_gcm_hist, pre_std_son_gcm_hist, pre_std_annual_gcm_hist, pre_djf_gcm_hist, pre_mam_gcm_hist, pre_jja_gcm_hist, pre_son_gcm_hist, pre_annual_gcm_hist = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_std_djf_rcm_rcp26, pre_std_mam_rcm_rcp26, pre_std_jja_rcm_rcp26, pre_std_son_rcm_rcp26, pre_std_annual_rcm_rcp26, pre_djf_rcm_rcp26, pre_mam_rcm_rcp26, pre_jja_rcm_rcp26, pre_son_rcm_rcp26, pre_annual_rcm_rcp26 = import_rcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_std_djf_gcm_rcp26, pre_std_mam_gcm_rcp26, pre_std_jja_gcm_rcp26, pre_std_son_gcm_rcp26, pre_std_annual_gcm_rcp26, pre_djf_gcm_rcp26, pre_mam_gcm_rcp26, pre_jja_gcm_rcp26, pre_son_gcm_rcp26, pre_annual_gcm_rcp26 = import_gcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_std_djf_rcm_rcp85, pre_std_mam_rcm_rcp85, pre_std_jja_rcm_rcp85, pre_std_son_rcm_rcp85, pre_std_annual_rcm_rcp85, pre_djf_rcm_rcp85, pre_mam_rcm_rcp85, pre_jja_rcm_rcp85, pre_son_rcm_rcp85, pre_annual_rcm_rcp85 = import_rcm('pr', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, pre_std_djf_gcm_rcp85, pre_std_mam_gcm_rcp85, pre_std_jja_gcm_rcp85, pre_std_son_gcm_rcp85, pre_std_annual_gcm_rcp85, pre_djf_gcm_rcp85, pre_mam_gcm_rcp85, pre_jja_gcm_rcp85, pre_son_gcm_rcp85, pre_annual_gcm_rcp85 = import_gcm('pr', 'amz_neb', 'rcp85', '2080-2099')

lat, lon, tas_std_djf_rcm_hist, tas_std_mam_rcm_hist, tas_std_jja_rcm_hist, tas_std_son_rcm_hist, tas_std_annual_rcm_hist, tas_djf_rcm_hist, tas_mam_rcm_hist, tas_jja_rcm_hist, tas_son_rcm_hist, tas_annual_rcm_hist = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_gcm_hist, tas_std_mam_gcm_hist, tas_std_jja_gcm_hist, tas_std_son_gcm_hist, tas_std_annual_gcm_hist, tas_djf_gcm_hist, tas_mam_gcm_hist, tas_jja_gcm_hist, tas_son_gcm_hist, tas_annual_gcm_hist = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_std_djf_rcm_rcp26, tas_std_mam_rcm_rcp26, tas_std_jja_rcm_rcp26, tas_std_son_rcm_rcp26, tas_std_annual_rcm_rcp26, tas_djf_rcm_rcp26, tas_mam_rcm_rcp26, tas_jja_rcm_rcp26, tas_son_rcm_rcp26, tas_annual_rcm_rcp26 = import_rcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_std_djf_gcm_rcp26, tas_std_mam_gcm_rcp26, tas_std_jja_gcm_rcp26, tas_std_son_gcm_rcp26, tas_std_annual_gcm_rcp26, tas_djf_gcm_rcp26, tas_mam_gcm_rcp26, tas_jja_gcm_rcp26, tas_son_gcm_rcp26, tas_annual_gcm_rcp26 = import_gcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_std_djf_rcm_rcp85, tas_std_mam_rcm_rcp85, tas_std_jja_rcm_rcp85, tas_std_son_rcm_rcp85, tas_std_annual_rcm_rcp85, tas_djf_rcm_rcp85, tas_mam_rcm_rcp85, tas_jja_rcm_rcp85, tas_son_rcm_rcp85, tas_annual_rcm_rcp85 = import_rcm('tas', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, tas_std_djf_gcm_rcp85, tas_std_mam_gcm_rcp85, tas_std_jja_gcm_rcp85, tas_std_son_gcm_rcp85, tas_std_annual_gcm_rcp85, tas_djf_gcm_rcp85, tas_mam_gcm_rcp85, tas_jja_gcm_rcp85, tas_son_gcm_rcp85, tas_annual_gcm_rcp85 = import_gcm('tas', 'amz_neb', 'rcp85', '2080-2099')

# Compute change from rcp and historical period
pre_djf_rcm_rcp26_hist = compute_relative_change(pre_djf_rcm_rcp26, pre_djf_rcm_hist)
pre_mam_rcm_rcp26_hist = compute_relative_change(pre_mam_rcm_rcp26, pre_mam_rcm_hist)
pre_jja_rcm_rcp26_hist = compute_relative_change(pre_jja_rcm_rcp26, pre_jja_rcm_hist)
pre_son_rcm_rcp26_hist = compute_relative_change(pre_son_rcm_rcp26, pre_son_rcm_hist)
pre_annual_rcm_rcp26_hist = compute_relative_change(pre_annual_rcm_rcp26, pre_annual_rcm_hist)
pre_djf_rcm_rcp85_hist = compute_relative_change(pre_djf_rcm_rcp85, pre_djf_rcm_hist)
pre_mam_rcm_rcp85_hist = compute_relative_change(pre_mam_rcm_rcp85, pre_mam_rcm_hist)
pre_jja_rcm_rcp85_hist = compute_relative_change(pre_jja_rcm_rcp85, pre_jja_rcm_hist)
pre_son_rcm_rcp85_hist = compute_relative_change(pre_son_rcm_rcp85, pre_son_rcm_hist)
pre_annual_rcm_rcp85_hist = compute_relative_change(pre_annual_rcm_rcp85, pre_annual_rcm_hist)
pre_djf_gcm_rcp26_hist = compute_relative_change(pre_djf_gcm_rcp26, pre_djf_gcm_hist)
pre_mam_gcm_rcp26_hist = compute_relative_change(pre_mam_gcm_rcp26, pre_mam_gcm_hist)
pre_jja_gcm_rcp26_hist = compute_relative_change(pre_jja_gcm_rcp26, pre_jja_gcm_hist)
pre_son_gcm_rcp26_hist = compute_relative_change(pre_son_gcm_rcp26, pre_son_gcm_hist)
pre_annual_gcm_rcp26_hist = compute_relative_change(pre_annual_gcm_rcp26, pre_annual_gcm_hist)
pre_djf_gcm_rcp85_hist = compute_relative_change(pre_djf_gcm_rcp85, pre_djf_gcm_hist)
pre_mam_gcm_rcp85_hist = compute_relative_change(pre_mam_gcm_rcp85, pre_mam_gcm_hist)
pre_jja_gcm_rcp85_hist = compute_relative_change(pre_jja_gcm_rcp85, pre_jja_gcm_hist)
pre_son_gcm_rcp85_hist = compute_relative_change(pre_son_gcm_rcp85, pre_son_gcm_hist)
pre_annual_gcm_rcp85_hist = compute_relative_change(pre_annual_gcm_rcp85, pre_annual_gcm_hist)

tas_djf_rcm_rcp26_hist = np.nanmean(tas_djf_rcm_rcp26, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp26_hist = np.nanmean(tas_mam_rcm_rcp26, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp26_hist = np.nanmean(tas_jja_rcm_rcp26, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp26_hist = np.nanmean(tas_son_rcm_rcp26, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_annual_rcm_rcp26_hist = np.nanmean(tas_annual_rcm_rcp26, axis=0) - np.nanmean(tas_annual_rcm_hist, axis=0)
tas_djf_rcm_rcp85_hist = np.nanmean(tas_djf_rcm_rcp85, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp85_hist = np.nanmean(tas_mam_rcm_rcp85, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp85_hist = np.nanmean(tas_jja_rcm_rcp85, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp85_hist = np.nanmean(tas_son_rcm_rcp85, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_annual_rcm_rcp85_hist = np.nanmean(tas_annual_rcm_rcp85, axis=0) - np.nanmean(tas_annual_rcm_hist, axis=0)
tas_djf_gcm_rcp26_hist = tas_djf_gcm_rcp26 - tas_djf_gcm_hist
tas_mam_gcm_rcp26_hist = tas_mam_gcm_rcp26 - tas_mam_gcm_hist
tas_jja_gcm_rcp26_hist = tas_jja_gcm_rcp26 - tas_jja_gcm_hist
tas_son_gcm_rcp26_hist = tas_son_gcm_rcp26 - tas_son_gcm_hist
tas_annual_gcm_rcp26_hist = tas_annual_gcm_rcp26 - tas_annual_gcm_hist
tas_djf_gcm_rcp85_hist = tas_djf_gcm_rcp85 - tas_djf_gcm_hist
tas_mam_gcm_rcp85_hist = tas_mam_gcm_rcp85 - tas_mam_gcm_hist
tas_jja_gcm_rcp85_hist = tas_jja_gcm_rcp85 - tas_jja_gcm_hist
tas_son_gcm_rcp85_hist = tas_son_gcm_rcp85 - tas_son_gcm_hist
tas_annual_gcm_rcp85_hist = tas_annual_gcm_rcp85 - tas_annual_gcm_hist

# Compute ttest
# Precipitation
p_value_pre_djf_rcm_rcp26 = ttest(pre_djf_rcm_rcp26, pre_djf_rcm_hist, pre_std_djf_rcm_rcp26, pre_std_djf_rcm_hist)
p_value_pre_mam_rcm_rcp26 = ttest(pre_mam_rcm_rcp26, pre_mam_rcm_hist, pre_std_mam_rcm_rcp26, pre_std_mam_rcm_hist)
p_value_pre_jja_rcm_rcp26 = ttest(pre_jja_rcm_rcp26, pre_jja_rcm_hist, pre_std_jja_rcm_rcp26, pre_std_jja_rcm_hist)
p_value_pre_son_rcm_rcp26 = ttest(pre_son_rcm_rcp26, pre_son_rcm_hist, pre_std_son_rcm_rcp26, pre_std_son_rcm_hist)
p_value_pre_annual_rcm_rcp26 = ttest(pre_annual_rcm_rcp26, pre_annual_rcm_hist, pre_std_annual_rcm_rcp26, pre_std_annual_rcm_hist)

p_value_pre_djf_gcm_rcp26 = ttest(pre_djf_gcm_rcp26, pre_djf_gcm_hist, pre_std_djf_gcm_rcp26, pre_std_djf_gcm_hist)
p_value_pre_mam_gcm_rcp26 = ttest(pre_mam_gcm_rcp26, pre_mam_gcm_hist, pre_std_mam_gcm_rcp26, pre_std_mam_gcm_hist)
p_value_pre_jja_gcm_rcp26 = ttest(pre_jja_gcm_rcp26, pre_jja_gcm_hist, pre_std_jja_gcm_rcp26, pre_std_jja_gcm_hist)
p_value_pre_son_gcm_rcp26 = ttest(pre_son_gcm_rcp26, pre_son_gcm_hist, pre_std_son_gcm_rcp26, pre_std_son_gcm_hist)
p_value_pre_annual_gcm_rcp26 = ttest(pre_annual_gcm_rcp26, pre_annual_gcm_hist, pre_std_annual_gcm_rcp26, pre_std_annual_gcm_hist)

p_value_pre_djf_rcm_rcp85 = ttest(pre_djf_rcm_rcp85, pre_djf_rcm_hist, pre_std_djf_rcm_rcp85, pre_std_djf_rcm_hist)
p_value_pre_mam_rcm_rcp85 = ttest(pre_mam_rcm_rcp85, pre_mam_rcm_hist, pre_std_mam_rcm_rcp85, pre_std_mam_rcm_hist)
p_value_pre_jja_rcm_rcp85 = ttest(pre_jja_rcm_rcp85, pre_jja_rcm_hist, pre_std_jja_rcm_rcp85, pre_std_jja_rcm_hist)
p_value_pre_son_rcm_rcp85 = ttest(pre_son_rcm_rcp85, pre_son_rcm_hist, pre_std_son_rcm_rcp85, pre_std_son_rcm_hist)
p_value_pre_annual_rcm_rcp85 = ttest(pre_annual_rcm_rcp85, pre_annual_rcm_hist, pre_std_annual_rcm_rcp85, pre_std_annual_rcm_hist)

p_value_pre_djf_gcm_rcp85 = ttest(pre_djf_gcm_rcp85, pre_djf_gcm_hist, pre_std_djf_gcm_rcp85, pre_std_djf_gcm_hist)
p_value_pre_mam_gcm_rcp85 = ttest(pre_mam_gcm_rcp85, pre_mam_gcm_hist, pre_std_mam_gcm_rcp85, pre_std_mam_gcm_hist)
p_value_pre_jja_gcm_rcp85 = ttest(pre_jja_gcm_rcp85, pre_jja_gcm_hist, pre_std_jja_gcm_rcp85, pre_std_jja_gcm_hist)
p_value_pre_son_gcm_rcp85 = ttest(pre_son_gcm_rcp85, pre_son_gcm_hist, pre_std_son_gcm_rcp85, pre_std_son_gcm_hist)
p_value_pre_annual_gcm_rcp85 = ttest(pre_annual_gcm_rcp26, pre_annual_gcm_hist, pre_std_annual_gcm_rcp85, pre_std_annual_gcm_hist)

# Temperature
p_value_tas_djf_rcm_rcp26 = ttest(tas_djf_rcm_rcp26, tas_djf_rcm_hist, tas_std_djf_rcm_rcp26, tas_std_djf_rcm_hist)
p_value_tas_mam_rcm_rcp26 = ttest(tas_mam_rcm_rcp26, tas_mam_rcm_hist, tas_std_mam_rcm_rcp26, tas_std_mam_rcm_hist)
p_value_tas_jja_rcm_rcp26 = ttest(tas_jja_rcm_rcp26, tas_jja_rcm_hist, tas_std_jja_rcm_rcp26, tas_std_jja_rcm_hist)
p_value_tas_son_rcm_rcp26 = ttest(tas_son_rcm_rcp26, tas_son_rcm_hist, tas_std_son_rcm_rcp26, tas_std_son_rcm_hist)
p_value_tas_annual_rcm_rcp26 = ttest(tas_annual_rcm_rcp26, tas_annual_rcm_hist, tas_std_annual_rcm_rcp26, tas_std_annual_rcm_hist)

p_value_tas_djf_gcm_rcp26 = ttest(tas_djf_gcm_rcp26, tas_djf_gcm_hist, tas_std_djf_gcm_rcp26, tas_std_djf_gcm_hist)
p_value_tas_mam_gcm_rcp26 = ttest(tas_mam_gcm_rcp26, tas_mam_gcm_hist, tas_std_mam_gcm_rcp26, tas_std_mam_gcm_hist)
p_value_tas_jja_gcm_rcp26 = ttest(tas_jja_gcm_rcp26, tas_jja_gcm_hist, tas_std_jja_gcm_rcp26, tas_std_jja_gcm_hist)
p_value_tas_son_gcm_rcp26 = ttest(tas_son_gcm_rcp26, tas_son_gcm_hist, tas_std_son_gcm_rcp26, tas_std_son_gcm_hist)
p_value_tas_annual_gcm_rcp26 = ttest(tas_annual_gcm_rcp26, tas_annual_gcm_hist, tas_std_annual_gcm_rcp26, tas_std_annual_gcm_hist)

p_value_tas_djf_rcm_rcp85 = ttest(tas_djf_rcm_rcp85, tas_djf_rcm_hist, tas_std_djf_rcm_rcp85, tas_std_djf_rcm_hist)
p_value_tas_mam_rcm_rcp85 = ttest(tas_mam_rcm_rcp85, tas_mam_rcm_hist, tas_std_mam_rcm_rcp85, tas_std_mam_rcm_hist)
p_value_tas_jja_rcm_rcp85 = ttest(tas_jja_rcm_rcp85, tas_jja_rcm_hist, tas_std_jja_rcm_rcp85, tas_std_jja_rcm_hist)
p_value_tas_son_rcm_rcp85 = ttest(tas_son_rcm_rcp85, tas_son_rcm_hist, tas_std_son_rcm_rcp85, tas_std_son_rcm_hist)
p_value_tas_annual_rcm_rcp85 = ttest(tas_annual_rcm_rcp85, tas_annual_rcm_hist, tas_std_annual_rcm_rcp85, tas_std_annual_rcm_hist)

p_value_tas_djf_gcm_rcp85 = ttest(tas_djf_gcm_rcp85, tas_djf_gcm_hist, tas_std_djf_gcm_rcp85, tas_std_djf_gcm_hist)
p_value_tas_mam_gcm_rcp85 = ttest(tas_mam_gcm_rcp85, tas_mam_gcm_hist, tas_std_mam_gcm_rcp85, tas_std_mam_gcm_hist)
p_value_tas_jja_gcm_rcp85 = ttest(tas_jja_gcm_rcp85, tas_jja_gcm_hist, tas_std_jja_gcm_rcp85, tas_std_jja_gcm_hist)
p_value_tas_son_gcm_rcp85 = ttest(tas_son_gcm_rcp85, tas_son_gcm_hist, tas_std_son_gcm_rcp85, tas_std_son_gcm_hist)
p_value_tas_annual_gcm_rcp85 = ttest(tas_annual_gcm_rcp26, tas_annual_gcm_hist, tas_std_annual_gcm_rcp85, tas_std_annual_gcm_hist)

# Plot models 
fig = plt.figure(figsize=(8, 9))
levs1 = [-60, -40, -20, 20, 40, 60]
levs2 = [-90, -60, -30, 30, 60, 90]
levs3 = [0.5, 1, 1.5, 2, 2.5, 3]
levs4 = [1, 2, 4, 6, 8, 10]

ax = fig.add_subplot(8, 5, 1)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_rcm_rcp26 = ma.masked_where(p_value_pre_djf_rcm_rcp26 >= 0.1, p_value_pre_djf_rcm_rcp26) 
map.contourf(xx, yy, p_value_pre_djf_rcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 2)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_rcm_rcp26 = ma.masked_where(p_value_pre_mam_rcm_rcp26 >= 0.1, p_value_pre_mam_rcm_rcp26) 
map.contourf(xx, yy, p_value_pre_mam_rcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 3)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_rcm_rcp26 = ma.masked_where(p_value_pre_jja_rcm_rcp26 >= 0.1, p_value_pre_jja_rcm_rcp26) 
map.contourf(xx, yy, p_value_pre_jja_rcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 4)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_rcm_rcp26 = ma.masked_where(p_value_pre_son_rcm_rcp26 >= 0.1, p_value_pre_son_rcm_rcp26) 
map.contourf(xx, yy, p_value_pre_son_rcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 5)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_rcm_rcp26 = ma.masked_where(p_value_pre_annual_rcm_rcp26 >= 0.1, p_value_pre_annual_rcm_rcp26) 
map.contourf(xx, yy, p_value_pre_annual_rcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 6)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_gcm_rcp26 = ma.masked_where(p_value_pre_djf_gcm_rcp26 >= 0.1, p_value_pre_djf_gcm_rcp26) 
map.contourf(xx, yy, p_value_pre_djf_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 7)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_gcm_rcp26 = ma.masked_where(p_value_pre_mam_gcm_rcp26 >= 0.1, p_value_pre_mam_gcm_rcp26) 
map.contourf(xx, yy, p_value_pre_mam_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 8)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_gcm_rcp26 = ma.masked_where(p_value_pre_jja_gcm_rcp26 >= 0.1, p_value_pre_jja_gcm_rcp26) 
map.contourf(xx, yy, p_value_pre_jja_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 9)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_gcm_rcp26 = ma.masked_where(p_value_pre_son_gcm_rcp26 >= 0.1, p_value_pre_son_gcm_rcp26) 
map.contourf(xx, yy, p_value_pre_son_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 10)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_gcm_rcp26 = ma.masked_where(p_value_pre_annual_gcm_rcp26 >= 0.1, p_value_pre_annual_gcm_rcp26) 
map.contourf(xx, yy, p_value_pre_annual_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 11)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_rcm_rcp85 = ma.masked_where(p_value_pre_djf_rcm_rcp85 >= 0.1, p_value_pre_djf_rcm_rcp85) 
map.contourf(xx, yy, p_value_pre_djf_rcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 12)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_rcm_rcp85 = ma.masked_where(p_value_pre_mam_rcm_rcp85 >= 0.1, p_value_pre_mam_rcm_rcp85) 
map.contourf(xx, yy, p_value_pre_mam_rcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 13)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_rcm_rcp85 = ma.masked_where(p_value_pre_jja_rcm_rcp85  >= 0.1, p_value_pre_jja_rcm_rcp85) 
map.contourf(xx, yy, p_value_pre_jja_rcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 14)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_rcm_rcp85 = ma.masked_where(p_value_pre_son_rcm_rcp85 >= 0.1, p_value_pre_son_rcm_rcp85) 
map.contourf(xx, yy, p_value_pre_son_rcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 15)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_rcm_rcp85 = ma.masked_where(p_value_pre_annual_rcm_rcp85 >= 0.1, p_value_pre_annual_rcm_rcp85) 
map.contourf(xx, yy, p_value_pre_annual_rcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 16)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_pre_djf_gcm_rcp85 = ma.masked_where(p_value_pre_djf_gcm_rcp85 >= 0.1, p_value_pre_djf_gcm_rcp85) 
map.contourf(xx, yy, p_value_pre_djf_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 17)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_mam_gcm_rcp85 = ma.masked_where(p_value_pre_mam_gcm_rcp85 >= 0.1, p_value_pre_mam_gcm_rcp85) 
map.contourf(xx, yy, p_value_pre_mam_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 18)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_jja_gcm_rcp85 = ma.masked_where(p_value_pre_jja_gcm_rcp85 >= 0.1, p_value_pre_jja_gcm_rcp85) 
map.contourf(xx, yy, p_value_pre_jja_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 19)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_son_gcm_rcp85 = ma.masked_where(p_value_pre_son_gcm_rcp85 >= 0.1, p_value_pre_son_gcm_rcp85) 
map.contourf(xx, yy, p_value_pre_son_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 20)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_pre_annual_gcm_rcp85 = ma.masked_where(p_value_pre_annual_gcm_rcp85 >= 0.1, p_value_pre_annual_gcm_rcp85) 
map.contourf(xx, yy, p_value_pre_annual_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 21)
plt.title(u'U)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_rcp26 = ma.masked_where(p_value_tas_djf_rcm_rcp26 >= 0.1, p_value_tas_djf_rcm_rcp26) 
map.contourf(xx, yy, p_value_tas_djf_rcm_rcp26[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 22)
plt.title(u'V)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_rcm_rcp26 = ma.masked_where(p_value_tas_mam_rcm_rcp26 >= 0.1, p_value_tas_mam_rcm_rcp26) 
map.contourf(xx, yy, p_value_tas_mam_rcm_rcp26[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 23)
plt.title(u'W)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_rcp26 = ma.masked_where(p_value_tas_jja_rcm_rcp26 >= 0.1, p_value_tas_jja_rcm_rcp26) 
map.contourf(xx, yy, p_value_tas_jja_rcm_rcp26[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 24)
plt.title(u'X)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_rcm_rcp26 = ma.masked_where(p_value_tas_son_rcm_rcp26 >= 0.1, p_value_tas_son_rcm_rcp26) 
map.contourf(xx, yy, p_value_tas_son_rcm_rcp26[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 25)
plt.title(u'Y)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_rcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_rcm_rcp26 = ma.masked_where(p_value_tas_annual_rcm_rcp26 >= 0.1, p_value_tas_annual_rcm_rcp26) 
map.contourf(xx, yy, p_value_tas_annual_rcm_rcp26[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 26)
plt.title(u'Z)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm_rcp26 = ma.masked_where(p_value_tas_djf_gcm_rcp26 >= 0.1, p_value_tas_djf_gcm_rcp26) 
map.contourf(xx, yy, p_value_tas_djf_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 27)
plt.title(u'A.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_gcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_gcm_rcp26 = ma.masked_where(p_value_tas_mam_gcm_rcp26 >= 0.1, p_value_tas_mam_gcm_rcp26) 
map.contourf(xx, yy, p_value_tas_mam_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 28)
plt.title(u'B.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm_rcp26 = ma.masked_where(p_value_tas_jja_gcm_rcp26 >= 0.1, p_value_tas_jja_gcm_rcp26) 
map.contourf(xx, yy, p_value_tas_jja_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 29)
plt.title(u'C.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_gcm_rcp26 = ma.masked_where(p_value_tas_son_gcm_rcp26 >= 0.1, p_value_tas_son_gcm_rcp26) 
map.contourf(xx, yy, p_value_tas_son_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 30)
plt.title(u'D.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_gcm_rcp26_hist, levels=levs3, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_gcm_rcp26 = ma.masked_where(p_value_tas_annual_gcm_rcp26 >= 0.1, p_value_tas_annual_gcm_rcp26) 
map.contourf(xx, yy, p_value_tas_annual_gcm_rcp26, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 31)
plt.title(u'E.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_rcm_rcp85 = ma.masked_where(p_value_tas_djf_rcm_rcp85 >= 0.1, p_value_tas_djf_rcm_rcp85) 
map.contourf(xx, yy, p_value_tas_djf_rcm_rcp85[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 32)
plt.title(u'F.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_rcm_rcp85 = ma.masked_where(p_value_tas_mam_rcm_rcp85 >= 0.1, p_value_tas_mam_rcm_rcp85) 
map.contourf(xx, yy, p_value_tas_mam_rcm_rcp85[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 33)
plt.title(u'G.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_rcm_rcp85 = ma.masked_where(p_value_tas_jja_rcm_rcp85 >= 0.1, p_value_tas_jja_rcm_rcp85) 
map.contourf(xx, yy, p_value_tas_jja_rcm_rcp85[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 34)
plt.title(u'H.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_rcm_rcp85 = ma.masked_where(p_value_tas_son_rcm_rcp85 >= 0.1, p_value_tas_son_rcm_rcp85) 
map.contourf(xx, yy, p_value_tas_son_rcm_rcp85[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 35)
plt.title(u'I.1)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_rcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_rcm_rcp85 = ma.masked_where(p_value_tas_annual_rcm_rcp85 >= 0.1, p_value_tas_annual_rcm_rcp85) 
map.contourf(xx, yy, p_value_tas_annual_rcm_rcp85[0,:,:], colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 36)
plt.title(u'J.1)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value_tas_djf_gcm_rcp85 = ma.masked_where(p_value_tas_djf_gcm_rcp85 >= 0.1, p_value_tas_djf_gcm_rcp85) 
map.contourf(xx, yy, p_value_tas_djf_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 37)
plt.title(u'K.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_gcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_mam_gcm_rcp85 = ma.masked_where(p_value_tas_mam_gcm_rcp85 >= 0.1, p_value_tas_mam_gcm_rcp85) 
map.contourf(xx, yy, p_value_tas_mam_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 38)
plt.title(u'L.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_jja_gcm_rcp85 = ma.masked_where(p_value_tas_jja_gcm_rcp85 >= 0.1, p_value_tas_jja_gcm_rcp85) 
map.contourf(xx, yy, p_value_tas_jja_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 39)
plt.title(u'M.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_son_gcm_rcp85 = ma.masked_where(p_value_tas_son_gcm_rcp85 >= 0.1, p_value_tas_son_gcm_rcp85) 
map.contourf(xx, yy, p_value_tas_son_gcm_rcp85, colors='none', hatches=['....'])

ax = fig.add_subplot(8, 5, 40)
plt.title(u'N.1)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_gcm_rcp85_hist, levels=levs4, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
p_value_tas_annual_gcm_rcp85 = ma.masked_where(p_value_tas_annual_gcm_rcp85 >= 0.1, p_value_tas_annual_gcm_rcp85) 
map.contourf(xx, yy, p_value_tas_annual_gcm_rcp85, colors='none', hatches=['....'])

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.10, hspace=0.10)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()









