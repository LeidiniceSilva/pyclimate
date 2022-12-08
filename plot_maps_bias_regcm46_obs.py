# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot bias maps from regcm46 and obs database"

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

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_obs(area, obs, time):
	
	param = 'precip' 
	date  = '2001-2010'

	path  = '/home/nice/Documentos/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(value, axis=0)
	std_obs = np.std(value, axis=0)

	return lat, lon, mean_obs, std_obs
	
	
def import_sim(area, exp, time):
	
	param = 'pr' 
	date  = '2001-2010'

	path  = '/home/nice/Documentos/dataset/rcm/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, time, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_sim = np.nanmean(value, axis=0)
	std_sim = np.std(value, axis=0)

	return lat, lon, mean_sim, std_sim
	

def ttest(mean_sample1, mean_sample2, std_sample1):

	# Calculate t statistics	
	p1 = mean_sample1 - mean_sample2 
	p2 = std_sample1 / np.sqrt(10)

	ttest = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(ttest, df=10)

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
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)

	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
def plot_maps_bias(bias_exp1_djf, bias_exp2_djf, bias_exp1_jja, bias_exp2_jja, bias_exp1_ann, bias_exp2_ann, p_djf_exp1, p_djf_exp2, p_jja_exp1, p_jja_exp2, p_ann_exp1, p_ann_exp2, diff_djf, diff_jja, diff_ann):
		
	fig = plt.figure(figsize=(8,4))
	levs = [-8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8]

	ax = fig.add_subplot(331)
	plt.title(u'A) Reg_Holtslag - GPCP DJF', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_djf, levels=levs, latlon=True, cmap=cm.BrBG) 
	p_djf_exp1 = ma.masked_where(p_djf_exp1 >= 0.05, p_djf_exp1) 
	map.contourf(xx, yy, p_djf_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	ax = fig.add_subplot(332)
	plt.title(u'B) Reg_UW-PBL - GPCP DJF', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_djf, levels=levs, latlon=True, cmap=cm.BrBG)
	p_djf_exp2 = ma.masked_where(p_djf_exp2 >= 0.05, p_djf_exp2) 
	map.contourf(xx, yy, p_djf_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(333)
	plt.title(u'C) Reg_Holtslag - Reg_UW-PBL', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_djf, levels=levs, latlon=True, cmap=cm.PiYG)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(334)
	plt.title(u'D) Reg_Holtslag - GPCP JJA', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_jja, levels=levs, latlon=True, cmap=cm.BrBG)
	p_jja_exp1 = ma.masked_where(p_jja_exp1 >= 0.05, p_jja_exp1) 
	map.contourf(xx, yy, p_jja_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	ax = fig.add_subplot(335)
	plt.title(u'E) Reg_UW-PBL - GPCP JJA', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_jja, levels=levs, latlon=True, cmap=cm.BrBG)
	p_jja_exp2 = ma.masked_where(p_jja_exp2 >= 0.05, p_jja_exp2) 
	map.contourf(xx, yy, p_jja_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(336)
	plt.title(u'F) Reg_Holtslag - Reg_UW-PBL', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_jja, levels=levs, latlon=True, cmap=cm.PiYG)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(337)
	plt.title(u'G) Reg_Holtslag - GPCP ANN', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_ann, levels=levs, latlon=True, cmap=cm.BrBG)
	p_ann_exp1 = ma.masked_where(p_ann_exp1 >= 0.05, p_ann_exp1) 
	map.contourf(xx, yy, p_ann_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(338)
	plt.title(u'H) Reg_UW-PBL - GPCP ANN', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_ann, levels=levs, latlon=True, cmap=cm.BrBG, extend='both')
	p_ann_exp2 = ma.masked_where(p_ann_exp2 >= 0.05, p_ann_exp2) 
	map.contourf(xx, yy, p_ann_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
	bounds = [-8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8]
	cbar = fig.colorbar(plt_maps_bias, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
	cbar.set_label(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
	cbar.ax.tick_params(labelsize=8)  

	ax = fig.add_subplot(339)
	plt.title(u'I) Reg_Holtslag - Reg_UW-PBL', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_ann, levels=levs, latlon=True, cmap=cm.PiYG, extend='both')
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	cb_ax = fig.add_axes([1, 0.2, 0.016, 0.6])
	bounds = [-8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8]
	cbar = fig.colorbar(plt_maps_bias, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
	cbar.set_label(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
	cbar.ax.tick_params(labelsize=8)  

	return plt_maps_bias


# Import regcm exps and obs database 
lat, lon, djf_obs, djf_obs_std = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'djf')
lat, lon, djf_exp1, djf_exp1_std = import_sim(u'amz_neb', u'regcm_exp1', 'djf')
lat, lon, djf_exp2, djf_exp2_std = import_sim(u'amz_neb', u'regcm_exp2', 'djf')

lat, lon, jja_obs, jja_obs_std = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'jja')
lat, lon, jja_exp1, jja_exp1_std = import_sim(u'amz_neb', u'regcm_exp1', 'jja')
lat, lon, jja_exp2, jja_exp2_std = import_sim(u'amz_neb', u'regcm_exp2', 'jja')

lat, lon, ann_obs, ann_obs_std = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'ann')
lat, lon, ann_exp1, ann_exp1_std = import_sim(u'amz_neb', u'regcm_exp1', 'ann')
lat, lon, ann_exp2, ann_exp2_std = import_sim(u'amz_neb', u'regcm_exp2', 'ann')

# Compute and plot bias from regcm exps and obs database 
bias_exp1_djf = djf_exp1 - djf_obs
bias_exp1_jja = jja_exp1 - jja_obs
bias_exp1_ann = ann_exp1 - ann_obs

bias_exp2_djf = djf_exp2 - djf_obs
bias_exp2_jja = jja_exp2 - jja_obs
bias_exp2_ann = ann_exp2 - ann_obs

diff_djf = djf_exp1 - djf_exp2
diff_jja = jja_exp1 - jja_exp2
diff_ann = ann_exp1 - ann_exp2

# Compute ttest from models and obs database 
p_djf_exp1 = ttest(djf_obs, djf_exp1, djf_obs_std)
p_djf_exp2 = ttest(djf_obs, djf_exp2, djf_obs_std)
p_jja_exp1 = ttest(jja_obs, jja_exp1, jja_obs_std)
p_jja_exp2 = ttest(jja_obs, jja_exp2, jja_obs_std)
p_ann_exp1 = ttest(djf_obs, ann_exp1, ann_obs_std)
p_ann_exp2 = ttest(djf_obs, ann_exp2, ann_obs_std)

# Plot maps with the function
plt_map = plot_maps_bias(bias_exp1_djf, bias_exp2_djf, bias_exp1_jja, bias_exp2_jja, bias_exp1_ann, bias_exp2_ann, p_djf_exp1, p_djf_exp2, p_jja_exp1, p_jja_exp2, p_ann_exp1, p_ann_exp2, diff_djf, diff_jja, diff_ann)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_bias_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



