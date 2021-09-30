# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot annual climatology from cmip6"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import gridspec

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from sklearn import metrics
from scipy.stats import norm
from matplotlib.font_manager import FontProperties


def import_obs(param, area):
	
	path  = '/home/nice/Documents/dataset/obs/cmip6'
	arq   = '{0}/{1}_{2}_cru_ts4.05_obs_mon_1961-2014_lonlat.nc'.format(path, param, area)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	obs_data1 = np.nanmean(np.nanmean(value[0:44,:,:], axis=1), axis=1)
	obs_clim1 = []
	for mon in range(1, 12 + 1):
		obs1 = np.nanmean(obs_data1[mon::12], axis=0)
		obs_clim1.append(obs1)

	obs_data2 = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
	obs_clim2 = []
	for mon in range(1, 12 + 1):
		obs2 = np.nanmean(obs_data2[mon::12], axis=0)
		obs_clim2.append(obs2)
			
	return obs_clim1, obs_clim2
	

def import_cmip(param, area, cmip, exp, date):
	
	path  = '/home/nice/Documents/dataset/gcm/cmip6/{0}'.format(cmip)
	arq   = '{0}/{1}_{2}_Amon_ensmean_{3}_historical_{4}_{5}_lonlat_mask.nc'.format(path, param, area, cmip, exp, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	sim_data = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
	sim_clim = []
	for mon in range(1, 12 + 1):
		sim = np.nanmean(sim_data[mon::12], axis=0)
		sim_clim.append(sim)
	
	return sim_clim	
	              
               
# Import regcm exps and obs database 
pre_clim1_samz_cru, pre_clim2_samz_cru = import_obs(u'pre', 'SAMZ')
pre_clim1_slpb_cru, pre_clim2_slpb_cru = import_obs(u'pre', 'SLPB')
pre_clim1_neb_cru,  pre_clim2_neb_cru  = import_obs(u'pre', 'NEB')

pre_clim_samz_cmip5 = import_cmip(u'pr', 'SAMZ', u'cmip5', 'r1i1p1', '1961-2005')
pre_clim_slpb_cmip5 = import_cmip(u'pr', 'SLPB', u'cmip5', 'r1i1p1', '1961-2005')
pre_clim_neb_cmip5  = import_cmip(u'pr', 'NEB', u'cmip5', 'r1i1p1', '1961-2005')

pre_clim_samz_cmip6 = import_cmip(u'pr', 'SAMZ', u'cmip6', 'r1i1p1f1_gn', '1961-2014')
pre_clim_slpb_cmip6 = import_cmip(u'pr', 'SLPB', u'cmip6', 'r1i1p1f1_gn', '1961-2014')
pre_clim_neb_cmip6  = import_cmip(u'pr', 'NEB', u'cmip6', 'r1i1p1f1_gn', '1961-2014')

pre_clim1_samz_cru, clim2_samz_cru = import_obs(u'tmp', 'SAMZ')
pre_clim1_slpb_cru, clim2_slpb_cru = import_obs(u'tmp', 'SLPB')
pre_clim1_neb_cru,  clim2_neb_cru  = import_obs(u'tmp', 'NEB')

pre_clim_samz_cmip5 = import_cmip(u'tas', 'SAMZ', u'cmip5', 'r1i1p1', '1961-2005')
pre_clim_slpb_cmip5 = import_cmip(u'tas', 'SLPB', u'cmip5', 'r1i1p1', '1961-2005')
pre_clim_neb_cmip5  = import_cmip(u'tas', 'NEB', u'cmip5', 'r1i1p1', '1961-2005')

pre_clim_samz_cmip6 = import_cmip(u'tas', 'SAMZ', u'cmip6', 'r1i1p1f1_gn', '1961-2014')
pre_clim_slpb_cmip6 = import_cmip(u'tas', 'SLPB', u'cmip6', 'r1i1p1f1_gn', '1961-2014')
pre_clim_neb_cmip6  = import_cmip(u'tas', 'NEB', u'cmip6', 'r1i1p1f1_gn', '1961-2014')

# Plot regcm exps and obs database 
fig = plt.figure(figsize=(10, 8)) 
gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1]) 
time = np.arange(0.5, 12 + 0.5)
time1 = np.arange(0, 45)
time2 = np.arange(0, 54)

ax = plt.subplot(gs[0])
annual_cycle1 = ax.plot(time2, ann_samz_cru,   linewidth=1.5, color='black', linestyle='-', label='CRU 61-14')
annual_cycle2 = ax.plot(time1, ann_samz_cmip5, linewidth=1.5, color='red',   linestyle='-', label='CMIP5 MME')
annual_cycle3 = ax.plot(time2, ann_samz_cmip6, linewidth=1.5, color='blue',  linestyle='-', label='CMIP6 MME')
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time2, ('1961', '', '', '', '1965', '', '', '', '1969', '', '', '', '1973', '', '', '', '1977', '', '', '', '1981', '', '', '', '1985', '', '', '', '1989', '', '', '', '1993', '', '', '', '1997', '', '', '', '2001', '', '', '', '2005', '', '', '', '2009', '', '', '', '2013', ''), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
plt.axvline(44, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(fontsize=6, loc=9, ncol=3, frameon=False)

ax = plt.subplot(gs[1])
annual_cycle1 = ax.plot(time, clim1_samz_cru,  linewidth=1.5, color='black', marker='.', linestyle='--', label='CRU 61-05')
annual_cycle2 = ax.plot(time, clim2_samz_cru,  linewidth=1.5, color='black', marker='.', linestyle='-',  label='CRU 61-14')
annual_cycle3 = ax.plot(time, clim_samz_cmip5, linewidth=1.5, color='red',   marker='.', linestyle='-',  label='CMIP5 MME')
annual_cycle3 = ax.plot(time, clim_samz_cmip6, linewidth=1.5, color='blue',  marker='.', linestyle='-',  label='CMP6 MME')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(fontsize=6, loc=9, ncol=2, frameon=False)

ax = plt.subplot(gs[2])
annual_cycle1 = ax.plot(time2, ann_slpb_cru,   linewidth=1.5, color='black', linestyle='-', label='CRU 61-14')
annual_cycle2 = ax.plot(time1, ann_slpb_cmip5, linewidth=1.5, color='red',   linestyle='-', label='CMIP5 MME')
annual_cycle3 = ax.plot(time2, ann_slpb_cmip6, linewidth=1.5, color='blue',  linestyle='-', label='CMIP6 MME')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time2, ('1961', '', '', '', '1965', '', '', '', '1969', '', '', '', '1973', '', '', '', '1977', '', '', '', '1981', '', '', '', '1985', '', '', '', '1989', '', '', '', '1993', '', '', '', '1997', '', '', '', '2001', '', '', '', '2005', '', '', '', '2009', '', '', '', '2013', ''), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
plt.axvline(44, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)

ax = plt.subplot(gs[3])
annual_cycle1 = ax.plot(time, clim1_slpb_cru,  linewidth=1.5, color='black', marker='.', linestyle='--', label='CRU 61-05')
annual_cycle2 = ax.plot(time, clim2_slpb_cru,  linewidth=1.5, color='black', marker='.', linestyle='-',  label='CRU 61-14')
annual_cycle3 = ax.plot(time, clim_slpb_cmip5, linewidth=1.5, color='red',   marker='.', linestyle='-',  label='CMIP5 MME')
annual_cycle3 = ax.plot(time, clim_slpb_cmip6, linewidth=1.5, color='blue',  marker='.', linestyle='-',  label='CMP6 MME')
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

ax = plt.subplot(gs[4])
annual_cycle1 = ax.plot(time2, ann_neb_cru,   linewidth=1.5, color='black', linestyle='-', label='CRU 61-14')
annual_cycle2 = ax.plot(time1, ann_neb_cmip5, linewidth=1.5, color='red',   linestyle='-', label='CMIP5 MME')
annual_cycle3 = ax.plot(time2, ann_neb_cmip6, linewidth=1.5, color='blue',  linestyle='-', label='CMIP6 MME')
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time2, ('1961', '', '', '', '1965', '', '', '', '1969', '', '', '', '1973', '', '', '', '1977', '', '', '', '1981', '', '', '', '1985', '', '', '', '1989', '', '', '', '1993', '', '', '', '1997', '', '', '', '2001', '', '', '', '2005', '', '', '', '2009', '', '', '', '2013', ''), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)
plt.axvline(44, linewidth=1., linestyle='--', color='black')

ax = plt.subplot(gs[5])
annual_cycle1 = ax.plot(time, clim1_neb_cru,  linewidth=1.5, color='black', marker='.', linestyle='--', label='CRU 61-05')
annual_cycle2 = ax.plot(time, clim2_neb_cru,  linewidth=1.5, color='black', marker='.', linestyle='-',  label='CRU 61-14')
annual_cycle3 = ax.plot(time, clim_neb_cmip5, linewidth=1.5, color='red',   marker='.', linestyle='-',  label='CMIP5 MME')
annual_cycle3 = ax.plot(time, clim_neb_cmip6, linewidth=1.5, color='blue',  marker='.', linestyle='-',  label='CMP6 MME')
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2), fontsize=8)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/paper_cmip6/figs'
name_out = 'pyplt_clim_cmip6_1961-2014_pre.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






