# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot portrait diagram from regcm47 and hadgem models and obs database"

import os
import netCDF4
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties
from comp_statist_indices import compute_ioa


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_obs = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_obs = season_obs[3:80:4]
	mam_obs = season_obs[0:80:4]
	jja_obs = season_obs[1:80:4]
	son_obs = season_obs[2:80:4]

	return annual_obs, djf_obs, mam_obs, jja_obs, son_obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	mam_rcm = season_rcm[0:80:4]
	jja_rcm = season_rcm[1:80:4]
	son_rcm = season_rcm[2:80:4]

	return annual_rcm, djf_rcm, mam_rcm, jja_rcm, son_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_gcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_gcm = season_gcm[3:80:4]
	mam_gcm = season_gcm[0:80:4]
	jja_gcm = season_gcm[1:80:4]
	son_gcm = season_gcm[2:80:4]

	return annual_gcm, djf_gcm, mam_gcm, jja_gcm, son_gcm

	              
# Import models and obs database 
# Precipitation
annual_cru_pre_samz, djf_cru_pre_samz, mam_cru_pre_samz, jja_cru_pre_samz, son_cru_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
annual_rcm_pre_samz, djf_rcm_pre_samz, mam_rcm_pre_samz, jja_rcm_pre_samz, son_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
annual_gcm_pre_samz, djf_gcm_pre_samz, mam_gcm_pre_samz, jja_gcm_pre_samz, son_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

annual_cru_pre_eneb, djf_cru_pre_eneb, mam_cru_pre_eneb, jja_cru_pre_eneb, son_cru_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
annual_rcm_pre_eneb, djf_rcm_pre_eneb, mam_rcm_pre_eneb, jja_rcm_pre_eneb, son_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
annual_gcm_pre_eneb, djf_gcm_pre_eneb, mam_gcm_pre_eneb, jja_gcm_pre_eneb, son_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

annual_cru_pre_matopiba, djf_cru_pre_matopiba, mam_cru_pre_matopiba, jja_cru_pre_matopiba, son_cru_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
annual_rcm_pre_matopiba, djf_rcm_pre_matopiba, mam_rcm_pre_matopiba, jja_rcm_pre_matopiba, son_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
annual_gcm_pre_matopiba, djf_gcm_pre_matopiba, mam_gcm_pre_matopiba, jja_gcm_pre_matopiba, son_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
annual_cru_tas_samz, djf_cru_tas_samz, mam_cru_tas_samz, jja_cru_tas_samz, son_cru_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
annual_rcm_tas_samz, djf_rcm_tas_samz, mam_rcm_tas_samz, jja_rcm_tas_samz, son_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
annual_gcm_tas_samz, djf_gcm_tas_samz, mam_gcm_tas_samz, jja_gcm_tas_samz, son_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

annual_cru_tas_eneb, djf_cru_tas_eneb, mam_cru_tas_eneb, jja_cru_tas_eneb, son_cru_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
annual_rcm_tas_eneb, djf_rcm_tas_eneb, mam_rcm_tas_eneb, jja_rcm_tas_eneb, son_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
annual_gcm_tas_eneb, djf_gcm_tas_eneb, mam_gcm_tas_eneb, jja_gcm_tas_eneb, son_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

annual_cru_tas_matopiba, djf_cru_tas_matopiba, mam_cru_tas_matopiba, jja_cru_tas_matopiba, son_cru_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
annual_rcm_tas_matopiba, djf_rcm_tas_matopiba, mam_rcm_tas_matopiba, jja_rcm_tas_matopiba, son_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
annual_gcm_tas_matopiba, djf_gcm_tas_matopiba, mam_gcm_tas_matopiba, jja_gcm_tas_matopiba, son_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Compute Nash–Sutcliffe Efficient Coefficient
# Precipitation
# RegCM
nse_djf_rcm_pre_samz = compute_ioa(djf_rcm_pre_samz, djf_cru_pre_samz)
nse_mam_rcm_pre_samz = compute_ioa(mam_rcm_pre_samz, mam_cru_pre_samz)
nse_jja_rcm_pre_samz = compute_ioa(jja_rcm_pre_samz, jja_cru_pre_samz)
nse_son_rcm_pre_samz = compute_ioa(son_rcm_pre_samz, son_cru_pre_samz)
nse_annual_rcm_pre_samz = compute_ioa(annual_rcm_pre_samz, annual_cru_pre_samz)
nse_djf_rcm_pre_eneb = compute_ioa(djf_rcm_pre_eneb, djf_cru_pre_eneb)
nse_mam_rcm_pre_eneb = compute_ioa(mam_rcm_pre_eneb, mam_cru_pre_eneb)
nse_jja_rcm_pre_eneb = compute_ioa(jja_rcm_pre_eneb, jja_cru_pre_eneb)
nse_son_rcm_pre_eneb = compute_ioa(son_rcm_pre_eneb, son_cru_pre_eneb)
nse_annual_rcm_pre_eneb = compute_ioa(annual_rcm_pre_eneb, annual_cru_pre_eneb)
nse_djf_rcm_pre_matopiba = compute_ioa(djf_rcm_pre_matopiba, djf_cru_pre_matopiba)
nse_mam_rcm_pre_matopiba = compute_ioa(mam_rcm_pre_matopiba, mam_cru_pre_matopiba)
nse_jja_rcm_pre_matopiba = compute_ioa(jja_rcm_pre_matopiba, jja_cru_pre_matopiba)
nse_son_rcm_pre_matopiba = compute_ioa(son_rcm_pre_matopiba, son_cru_pre_matopiba)
nse_annual_rcm_pre_matopiba = compute_ioa(annual_rcm_pre_matopiba, annual_cru_pre_matopiba)

# HadGEM2-ES
nse_djf_gcm_pre_samz = compute_ioa(djf_gcm_pre_samz, djf_cru_pre_samz)
nse_mam_gcm_pre_samz = compute_ioa(mam_gcm_pre_samz, mam_cru_pre_samz)
nse_jja_gcm_pre_samz = compute_ioa(jja_gcm_pre_samz, jja_cru_pre_samz)
nse_son_gcm_pre_samz = compute_ioa(son_gcm_pre_samz, son_cru_pre_samz)
nse_annual_gcm_pre_samz = compute_ioa(annual_gcm_pre_samz, annual_cru_pre_samz)
nse_djf_gcm_pre_eneb = compute_ioa(djf_gcm_pre_eneb, djf_cru_pre_eneb)
nse_mam_gcm_pre_eneb = compute_ioa(mam_gcm_pre_eneb, mam_cru_pre_eneb)
nse_jja_gcm_pre_eneb = compute_ioa(jja_gcm_pre_eneb, jja_cru_pre_eneb)
nse_son_gcm_pre_eneb = compute_ioa(son_gcm_pre_eneb, son_cru_pre_eneb)
nse_annual_gcm_pre_eneb = compute_ioa(annual_gcm_pre_eneb, annual_cru_pre_eneb)
nse_djf_gcm_pre_matopiba = compute_ioa(djf_gcm_pre_matopiba, djf_cru_pre_matopiba)
nse_mam_gcm_pre_matopiba = compute_ioa(mam_gcm_pre_matopiba, mam_cru_pre_matopiba)
nse_jja_gcm_pre_matopiba = compute_ioa(jja_gcm_pre_matopiba, jja_cru_pre_matopiba)
nse_son_gcm_pre_matopiba = compute_ioa(son_gcm_pre_matopiba, son_cru_pre_matopiba)
nse_annual_gcm_pre_matopiba = compute_ioa(annual_gcm_pre_matopiba, annual_cru_pre_matopiba)

# Temperature
# RegCM
nse_djf_rcm_tas_samz = compute_ioa(np.nanmean(djf_rcm_tas_samz, axis=1), djf_cru_tas_samz)
nse_mam_rcm_tas_samz = compute_ioa(np.nanmean(mam_rcm_tas_samz, axis=1), mam_cru_tas_samz)
nse_jja_rcm_tas_samz = compute_ioa(np.nanmean(jja_rcm_tas_samz, axis=1), jja_cru_tas_samz)
nse_son_rcm_tas_samz = compute_ioa(np.nanmean(son_rcm_tas_samz, axis=1), son_cru_tas_samz)
nse_annual_rcm_tas_samz = compute_ioa(np.nanmean(annual_rcm_tas_samz, axis=1), annual_cru_tas_samz)
nse_djf_rcm_tas_eneb = compute_ioa(np.nanmean(djf_rcm_tas_eneb, axis=1), djf_cru_tas_eneb)
nse_mam_rcm_tas_eneb = compute_ioa(np.nanmean(mam_rcm_tas_eneb, axis=1), mam_cru_tas_eneb)
nse_jja_rcm_tas_eneb = compute_ioa(np.nanmean(jja_rcm_tas_eneb, axis=1), jja_cru_tas_eneb)
nse_son_rcm_tas_eneb = compute_ioa(np.nanmean(son_rcm_tas_eneb, axis=1), son_cru_tas_eneb)
nse_annual_rcm_tas_eneb = compute_ioa(np.nanmean(annual_rcm_tas_eneb, axis=1), annual_cru_tas_eneb)
nse_djf_rcm_tas_matopiba = compute_ioa(np.nanmean(djf_rcm_tas_matopiba, axis=1), djf_cru_tas_matopiba)
nse_mam_rcm_tas_matopiba = compute_ioa(np.nanmean(mam_rcm_tas_matopiba, axis=1), mam_cru_tas_matopiba)
nse_jja_rcm_tas_matopiba = compute_ioa(np.nanmean(jja_rcm_tas_matopiba, axis=1), jja_cru_tas_matopiba)
nse_son_rcm_tas_matopiba = compute_ioa(np.nanmean(son_rcm_tas_matopiba, axis=1), son_cru_tas_matopiba)
nse_annual_rcm_tas_matopiba = compute_ioa(np.nanmean(annual_rcm_tas_matopiba, axis=1), annual_cru_tas_matopiba)

# HadGEM2-ES
nse_djf_gcm_tas_samz = compute_ioa(djf_gcm_tas_samz, djf_cru_tas_samz)
nse_mam_gcm_tas_samz = compute_ioa(mam_gcm_tas_samz, mam_cru_tas_samz)
nse_jja_gcm_tas_samz = compute_ioa(jja_gcm_tas_samz, jja_cru_tas_samz)
nse_son_gcm_tas_samz = compute_ioa(son_gcm_tas_samz, son_cru_tas_samz)
nse_annual_gcm_tas_samz = compute_ioa(annual_gcm_tas_samz, annual_cru_tas_samz)
nse_djf_gcm_tas_eneb = compute_ioa(djf_gcm_tas_eneb, djf_cru_tas_eneb)
nse_mam_gcm_tas_eneb = compute_ioa(mam_gcm_tas_eneb, mam_cru_tas_eneb)
nse_jja_gcm_tas_eneb = compute_ioa(jja_gcm_tas_eneb, jja_cru_tas_eneb)
nse_son_gcm_tas_eneb = compute_ioa(son_gcm_tas_eneb, son_cru_tas_eneb)
nse_annual_gcm_tas_eneb = compute_ioa(annual_gcm_tas_eneb, annual_cru_tas_eneb)
nse_djf_gcm_tas_matopiba = compute_ioa(djf_gcm_tas_matopiba, djf_cru_tas_matopiba)
nse_mam_gcm_tas_matopiba = compute_ioa(mam_gcm_tas_matopiba, mam_cru_tas_matopiba)
nse_jja_gcm_tas_matopiba = compute_ioa(jja_gcm_tas_matopiba, jja_cru_tas_matopiba)
nse_son_gcm_tas_matopiba = compute_ioa(son_gcm_tas_matopiba, son_cru_tas_matopiba)
nse_annual_gcm_tas_matopiba = compute_ioa(annual_gcm_tas_matopiba, annual_cru_tas_matopiba)

nse_rcm_pre = np.array([[nse_djf_rcm_pre_samz, nse_mam_rcm_pre_samz, nse_jja_rcm_pre_samz, nse_son_rcm_pre_samz, nse_annual_rcm_pre_samz],
[nse_djf_rcm_pre_eneb, nse_mam_rcm_pre_eneb, nse_jja_rcm_pre_eneb, nse_son_rcm_pre_eneb, nse_annual_rcm_pre_eneb], 
[nse_djf_rcm_pre_matopiba, nse_mam_rcm_pre_matopiba, nse_jja_rcm_pre_matopiba, nse_son_rcm_pre_matopiba, nse_annual_rcm_pre_matopiba]])

nse_gcm_pre = np.array([[nse_djf_gcm_pre_samz, nse_mam_gcm_pre_samz, nse_jja_gcm_pre_samz, nse_son_gcm_pre_samz, nse_annual_gcm_pre_samz],
[nse_djf_gcm_pre_eneb, nse_mam_gcm_pre_eneb, nse_jja_gcm_pre_eneb, nse_son_gcm_pre_eneb, nse_annual_gcm_pre_eneb], 
[nse_djf_gcm_pre_matopiba, nse_mam_gcm_pre_matopiba, nse_jja_gcm_pre_matopiba, nse_son_gcm_pre_matopiba, nse_annual_gcm_pre_matopiba]])

nse_rcm_tas = np.array([[nse_djf_rcm_tas_samz, nse_mam_rcm_tas_samz, nse_jja_rcm_tas_samz, nse_son_rcm_tas_samz, nse_annual_rcm_tas_samz],
[nse_djf_rcm_tas_eneb, nse_mam_rcm_tas_eneb, nse_jja_rcm_tas_eneb, nse_son_rcm_tas_eneb, nse_annual_rcm_tas_eneb], 
[nse_djf_rcm_tas_matopiba, nse_mam_rcm_tas_matopiba, nse_jja_rcm_tas_matopiba, nse_son_rcm_tas_matopiba, nse_annual_rcm_tas_matopiba]])

nse_gcm_tas = np.array([[nse_djf_gcm_tas_samz, nse_mam_gcm_tas_samz, nse_jja_gcm_tas_samz, nse_son_gcm_tas_samz, nse_annual_gcm_tas_samz],
[nse_djf_gcm_tas_eneb, nse_mam_gcm_tas_eneb, nse_jja_gcm_tas_eneb, nse_son_gcm_tas_eneb, nse_annual_gcm_tas_eneb], 
[nse_djf_gcm_tas_matopiba, nse_mam_gcm_tas_matopiba, nse_jja_gcm_tas_matopiba, nse_son_gcm_tas_matopiba, nse_annual_gcm_tas_matopiba]])

# Plot models and obs database 
fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
norm = colors.BoundaryNorm(boundaries=np.arange(0, 1, 0.1), ncolors=256)

xlabels = [u'DJF', u'MAM', u'JJA', u'SON', u'ANN']
ylabels = [u'MATOPIBA', u'ENEB', u'SAMZ']

# First column heatmaps with same colormap
pcm1 = axes[0, 0].pcolormesh(nse_rcm_pre, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='YlGnBu')
axes[0, 0].set_title(u'A) RegCM4.7', loc='left', fontweight='bold', fontsize=8)
axes[0, 0].set_xticks(np.arange(nse_rcm_pre.shape[1]) + 0.5)
axes[0, 0].set_yticks(np.arange(nse_rcm_pre.shape[0]) + 0.5)
axes[0, 0].set_yticklabels(ylabels, fontsize=8)
plt.setp(axes[0, 0].get_xticklabels(), visible=False)
axes[0, 0].tick_params(bottom = False)
for y in range(nse_rcm_pre.shape[0]):
    for x in range(nse_rcm_pre.shape[1]):
        axes[0, 0].text(x + 0.5, y + 0.5, '%.1f' % nse_rcm_pre[y, x],
                 ha="center", va="center", color='k',
                 )
                 
pcm2 = axes[1, 0].pcolormesh(nse_gcm_pre, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='YlGnBu')
axes[1, 0].set_title(u'B) HadGEM2-ES', loc='left', fontweight='bold', fontsize=8)
axes[1, 0].set_xticks(np.arange(nse_gcm_pre.shape[1]) + 0.5)
axes[1, 0].set_yticks(np.arange(nse_gcm_pre.shape[0]) + 0.5)
axes[1, 0].set_xticklabels(xlabels, fontsize=8)
axes[1, 0].set_yticklabels(ylabels, fontsize=8)
clb=fig.colorbar(pcm2, ax=axes[:, 0], extend='max', shrink=0.7)
clb.set_label('IOA de precipitação', rotation=90)
clb.ax.yaxis.set_label_position('left')
for y in range(nse_gcm_pre.shape[0]):
    for x in range(nse_gcm_pre.shape[1]):
        axes[1, 0].text(x + 0.5, y + 0.5, '%.1f' % nse_gcm_pre[y, x],
                 ha="center", va="center", color='k',
                 )
                 
# Second column heatmaps with same colormap
pcm3 = axes[0, 1].pcolormesh(nse_rcm_tas, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='YlOrRd')
axes[0, 1].set_title('C) RegCM4.7', loc='left', fontweight='bold', fontsize=8)
axes[0, 1].set_xticks(np.arange(nse_rcm_tas.shape[1]) + 0.5)
axes[0, 1].set_yticks(np.arange(nse_rcm_tas.shape[0]) + 0.5)
axes[0, 1].set_yticklabels(ylabels, fontsize=8)
axes[0, 1].tick_params(left = False)
axes[0, 1].tick_params(bottom = False)
plt.setp(axes[0, 1].get_xticklabels(), visible=False)
plt.setp(axes[0, 1].get_yticklabels(), visible=False)
for y in range(nse_rcm_tas.shape[0]):
    for x in range(nse_rcm_tas.shape[1]):
        axes[0, 1].text(x + 0.5, y + 0.5, '%.1f' % nse_rcm_tas[y, x],
                 ha="center", va="center", color='k',
                 )
                 
pcm4 = axes[1, 1].pcolormesh(nse_gcm_tas, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='YlOrRd')
axes[1, 1].set_title('D) HadGEM2-ES', loc='left', fontweight='bold', fontsize=8)
axes[1, 1].set_xticks(np.arange(nse_gcm_tas.shape[1]) + 0.5)
axes[1, 1].set_yticks(np.arange(nse_gcm_tas.shape[0]) + 0.5)
axes[1, 1].set_xticklabels(xlabels, fontsize=8)
axes[1, 1].set_yticklabels(ylabels, fontsize=8)
axes[1, 1].tick_params(left = False)
plt.setp(axes[1, 1].get_yticklabels(), visible=False)
clb=fig.colorbar(pcm4, ax=axes[:, 1], extend='max', shrink=0.7)
clb.set_label('IOA de temperatura', rotation=90)
clb.ax.yaxis.set_label_position('left')
for y in range(nse_gcm_tas.shape[0]):
    for x in range(nse_gcm_tas.shape[1]):
        axes[1, 1].text(x + 0.5, y + 0.5, '%.1f' % nse_gcm_tas[y, x],
                 ha="center", va="center", color='k',
                 )
                       
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_portrait_diagram_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


