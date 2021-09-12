# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology maps from from regcm46 and obs database"

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

	mean_obs = np.nanmean(value, axis=0)
	std_obs = np.std(value, axis=0)
	
	season_obs = value[2:120:3,:,:]
	djf_obs = np.nanmean(season_obs[3:40:4], axis=0)
	mam_obs = np.nanmean(season_obs[0:40:4], axis=0)
	jja_obs = np.nanmean(season_obs[1:40:4], axis=0)

	return lat, lon, mean_obs, std_obs, djf_obs, mam_obs, jja_obs
	
	
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

	mean_sim = np.nanmean(value, axis=0)
	std_sim = np.std(value, axis=0)
	
	season_sim = value[2:120:3,:,:]
	djf_sim = np.nanmean(season_sim[3:40:4], axis=0)
	mam_sim = np.nanmean(season_sim[0:40:4], axis=0)
	jja_sim = np.nanmean(season_sim[1:40:4], axis=0)

	return lat, lon, mean_sim, std_sim, djf_sim, mam_sim, jja_sim
	

def ttest(mean, std, sample):

	# Calculate t statistics
	p2= std / np.sqrt(40)
	p3 = sample / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(p3, df=40)
	
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
	
	
def plot_maps_clim(djf_obs, djf_exp1, djf_exp2, mam_obs, mam_exp1, mam_exp2, jja_obs, jja_exp1, jja_exp2, p_djf_obs, p_mam_obs, p_jja_obs, p_djf_exp1, p_mam_exp1, p_jja_exp1, p_djf_exp2, p_mam_exp2, p_jja_exp2):
		
	fig = plt.figure(figsize=(8,4))
	levs = [1, 3, 6, 9, 12, 15]

	ax = fig.add_subplot(331)
	plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, djf_obs, levels=levs, latlon=True, cmap=cm.YlGnBu) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	p_djf_obs = ma.masked_where(p_djf_obs >= 0.05, p_djf_obs) 
	map.contourf(xx, yy, p_djf_obs, colors='none', hatches=['....'])
	
	ax = fig.add_subplot(332)
	plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, djf_exp1, levels=levs, latlon=True, cmap=cm.YlGnBu)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_djf_exp1 = ma.masked_where(p_djf_exp1 >= 0.05, p_djf_exp1) 
	map.contourf(xx, yy, p_djf_exp1, colors='none', hatches=['....'])
	
	ax = fig.add_subplot(333)
	plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, djf_exp2, levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_djf_exp2 = ma.masked_where(p_djf_exp2 >= 0.05, p_djf_exp2) 
	map.contourf(xx, yy, p_djf_exp2, colors='none', hatches=['....'])
			
	ax = fig.add_subplot(334)
	plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, mam_obs, levels=levs, latlon=True, cmap=cm.YlGnBu) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	p_mam_obs = ma.masked_where(p_mam_obs >= 0.05, p_mam_obs) 
	map.contourf(xx, yy, p_mam_obs, colors='none', hatches=['....'])
			
	ax = fig.add_subplot(335)
	plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mam_exp1, levels=levs, latlon=True, cmap=cm.YlGnBu)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_mam_exp1 = ma.masked_where(p_mam_exp1 >= 0.05, p_mam_exp1) 
	map.contourf(xx, yy, p_mam_exp1, colors='none', hatches=['....'])
		
	ax = fig.add_subplot(336)
	plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mam_exp2, levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_mam_exp2 = ma.masked_where(p_mam_exp2 >= 0.05, p_mam_exp2) 
	map.contourf(xx, yy, p_mam_exp2, colors='none', hatches=['....'])
				
	ax = fig.add_subplot(337)
	plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, jja_obs, levels=levs, latlon=True, cmap=cm.YlGnBu)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[1,0,0,0], linewidth=0.5, color='black')
	p_jja_obs = ma.masked_where(p_jja_obs >= 0.05, p_jja_obs) 
	map.contourf(xx, yy, p_jja_obs, colors='none', hatches=['....'])
	
	ax = fig.add_subplot(338)
	plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, jja_exp1, levels=levs, latlon=True, cmap=cm.YlGnBu)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_jja_exp1 = ma.masked_where(p_jja_exp1 >= 0.05, p_jja_exp1) 
	map.contourf(xx, yy, p_jja_exp1, colors='none', hatches=['....'])
				
	ax = fig.add_subplot(339)
	plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, jja_exp2, levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=7, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=7, labels=[0,0,0,0], linewidth=0.5, color='black')
	p_jja_exp2 = ma.masked_where(p_jja_exp2 >= 0.05, p_jja_exp2) 
	map.contourf(xx, yy, p_jja_exp2, colors='none', hatches=['....'])
	
	return plt_maps_clim


# Import regcm exps and obs database 
lat, lon, mean_obs, std_obs, djf_obs, mam_obs, jja_obs = import_obs(u'amz_neb', u'gpcp_v2.3_obs')
lat, lon, mean_exp1, std_exp1, djf_exp1, mam_exp1, jja_exp1 = import_sim(u'amz_neb', u'regcm_exp1')
lat, lon, mean_exp2, std_exp2, djf_exp2, mam_exp2, jja_exp2 = import_sim(u'amz_neb', u'regcm_exp2')

# Compute and plot bias from regcm exps and obs database 
bias_exp1_djf = djf_exp1 - djf_obs
bias_exp1_mam = mam_exp1 - mam_obs
bias_exp1_jja = jja_exp1 - jja_obs

bias_exp2_djf = djf_exp2 - djf_obs
bias_exp2_mam = mam_exp2 - mam_obs
bias_exp2_jja = jja_exp2 - jja_obs

# Compute ttest from models and obs database 
p_djf_obs = ttest(mean_obs, std_obs, djf_obs)
p_mam_obs = ttest(mean_obs, std_obs, mam_obs)
p_jja_obs = ttest(mean_obs, std_obs, jja_obs)

p_djf_exp1 = ttest(mean_exp1, std_exp1, djf_exp1)
p_mam_exp1 = ttest(mean_exp1, std_exp1, mam_exp1)
p_jja_exp1 = ttest(mean_exp1, std_exp1, jja_exp1)

p_djf_exp2 = ttest(mean_exp1, std_exp2, djf_exp2)
p_mam_exp2 = ttest(mean_exp1, std_exp2, mam_exp2)
p_jja_exp2 = ttest(mean_exp1, std_exp2, jja_exp2)

# Plot maps with the function
plt_map = plot_maps_clim(djf_obs, djf_exp1, djf_exp2, mam_obs, mam_exp1, mam_exp2, jja_obs, jja_exp1, jja_exp2, p_djf_obs, p_mam_obs, p_jja_obs, p_djf_exp1, p_mam_exp1, p_jja_exp1, p_djf_exp2, p_mam_exp2, p_jja_exp2)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/figs'
name_out = 'pyplt_maps_clim_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



