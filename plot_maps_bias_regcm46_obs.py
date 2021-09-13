# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot bias maps from from regcm46 and obs database"

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

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_obs(area, obs):
	
	param = 'precip' # precip, pre or tmp
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	season_obs = value[2:120:3,:,:]
	djf_obs = np.nanmean(season_obs[3:40:4], axis=0)
	mam_obs = np.nanmean(season_obs[0:40:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:40:4], axis=0)

	djf_obs_std = np.std(season_obs[3:40:4], axis=0)
	mam_obs_std = np.std(season_obs[0:40:4], axis=0)
	jja_obs_std = np.std(season_obs[1:40:4], axis=0)

	return lat, lon, djf_obs, mam_obs, jja_obs, djf_obs_std, mam_obs_std, jja_obs_std
	

def import_sim(area, exp):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	season_sim = value[2:120:3,:,:]
	
	djf_sim = np.nanmean(season_sim[3:40:4], axis=0)
	mam_sim = np.nanmean(season_sim[0:40:4], axis=0)
	jja_sim = np.nanmean(season_sim[1:40:4], axis=0)

	djf_sim_std = np.std(season_sim[3:40:4], axis=0)
	mam_sim_std = np.std(season_sim[0:40:4], axis=0)
	jja_sim_std = np.std(season_sim[1:40:4], axis=0)
	
	return lat, lon, djf_sim, mam_sim, jja_sim, djf_sim_std, mam_sim_std, jja_sim_std
	

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
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)
	
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
def plot_maps_bias(bias_exp1_djf, bias_exp2_djf, bias_exp1_mam, bias_exp2_mam, bias_exp1_jja, bias_exp2_jja, p_djf_exp1, p_djf_exp2, p_mam_exp1, p_mam_exp2, p_jja_exp1, p_jja_exp2, diff_djf, diff_mam, diff_jja, p_djf, p_mam, p_jja):
		
	fig = plt.figure(figsize=(8,4))
	levs = [-6, -4, -2, 2, 4, 6]
	levs1 = [-3, -2, -1, 1, 2, 3]

	ax = fig.add_subplot(331)
	plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_djf, levels=levs, latlon=True, cmap=cm.RdBu) 
	p_djf_exp1 = ma.masked_where(p_djf_exp1 >= 0.05, p_djf_exp1) 
	map.contourf(xx, yy, p_djf_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	ax = fig.add_subplot(332)
	plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_djf, levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_djf_exp2 = ma.masked_where(p_djf_exp2 >= 0.05, p_djf_exp2) 
	map.contourf(xx, yy, p_djf_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(333)
	plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_djf, levels=levs1, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_djf = ma.masked_where(p_djf >= 0.05, p_djf) 
	map.contourf(xx, yy, p_djf, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(334)
	plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_mam, levels=levs, latlon=True, cmap=cm.RdBu)
	p_mam_exp1 = ma.masked_where(p_mam_exp1 >= 0.05, p_mam_exp1) 
	map.contourf(xx, yy, p_mam_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	ax = fig.add_subplot(335)
	plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_mam, levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_mam_exp2 = ma.masked_where(p_mam_exp2 >= 0.05, p_mam_exp2) 
	map.contourf(xx, yy, p_mam_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(336)
	plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_mam, levels=levs1, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_mam = ma.masked_where(p_mam >= 0.05, p_mam) 
	map.contourf(xx, yy, p_mam, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(337)
	plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_jja, levels=levs, latlon=True, cmap=cm.RdBu)
	p_jja_exp1 = ma.masked_where(p_jja_exp1 >= 0.05, p_jja_exp1) 
	map.contourf(xx, yy, p_jja_exp1, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(338)
	plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_jja, levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	p_jja_exp2 = ma.masked_where(p_jja_exp2 >= 0.05, p_jja_exp2) 
	map.contourf(xx, yy, p_jja_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(339)
	plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, diff_jja, levels=levs1, latlon=True, cmap=cm.PiYG)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_jja_exp2 = ma.masked_where(p_jja_exp2 >= 0.05, p_jja_exp2) 
	map.contourf(xx, yy, p_jja_exp2, colors='none', hatches=['....'])
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
			
	return plt_maps_bias


# Import regcm exps and obs database 
lat, lon, djf_obs, mam_obs, jja_obs, djf_obs_std, mam_obs_std, jja_obs_std = import_obs(u'amz_neb', u'gpcp_v2.3_obs')
lat, lon, djf_exp1, mam_exp1, jja_exp1, djf_exp1_std, mam_exp1_std, jja_exp1_std = import_sim(u'amz_neb', u'regcm_exp1')
lat, lon, djf_exp2, mam_exp2, jja_exp2, djf_exp2_std, mam_exp2_std, jja_exp2_std = import_sim(u'amz_neb', u'regcm_exp2')

# Compute and plot bias from regcm exps and obs database 
bias_exp1_djf = djf_exp1 - djf_obs
bias_exp1_mam = mam_exp1 - mam_obs
bias_exp1_jja = jja_exp1 - jja_obs

bias_exp2_djf = djf_exp2 - djf_obs
bias_exp2_mam = mam_exp2 - mam_obs
bias_exp2_jja = jja_exp2 - jja_obs

diff_djf = djf_exp1 - djf_exp2
diff_mam = mam_exp1 - mam_exp2
diff_jja = jja_exp1 - jja_exp2

# Compute ttest from models and obs database 
p_djf_exp1 = ttest(djf_exp1, djf_obs, djf_exp1_std, djf_obs_std)
p_djf_exp2 = ttest(djf_exp2, djf_obs, djf_exp2_std, djf_obs_std)
p_mam_exp1 = ttest(mam_exp1, mam_obs, mam_exp1_std, mam_obs_std)
p_mam_exp2 = ttest(mam_exp2, mam_obs, mam_exp2_std, mam_obs_std)
p_jja_exp1 = ttest(jja_exp1, mam_obs, jja_exp1_std, jja_obs_std)
p_jja_exp2 = ttest(jja_exp2, mam_obs, jja_exp2_std, jja_obs_std)
p_djf = ttest(jja_exp2, jja_exp1, jja_exp2_std, jja_exp1_std)
p_mam = ttest(jja_exp2, jja_exp1, jja_exp2_std, jja_exp1_std)
p_jja = ttest(jja_exp2, jja_exp1, jja_exp2_std, jja_exp1_std)

# Plot maps with the function
plt_map = plot_maps_bias(bias_exp1_djf, bias_exp2_djf, bias_exp1_mam, bias_exp2_mam, bias_exp1_jja, bias_exp2_jja, p_djf_exp1, p_djf_exp2, p_mam_exp1, p_mam_exp2, p_jja_exp1, p_jja_exp2, diff_djf, diff_mam, diff_jja, p_djf, p_mam, p_jja)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/figs'
name_out = 'pyplt_maps_bias_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



