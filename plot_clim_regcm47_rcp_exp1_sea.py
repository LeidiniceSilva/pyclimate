# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual climatology from regcm47 and hadgem models to rcp"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from scipy.stats import t
from scipy.stats import norm
from comp_statist_indices import compute_relative_change
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1))
	sea_rcm = np.nanmean(np.nanmean(value[0:240:3,:,:], axis=1), axis=1)
	djf_rcm = np.nanmean(sea_rcm[3:80:4])
	mam_rcm = np.nanmean(sea_rcm[0:80:4])
	jja_rcm = np.nanmean(sea_rcm[1:80:4])
	son_rcm = np.nanmean(sea_rcm[2:80:4])

	return djf_rcm, mam_rcm, jja_rcm, son_rcm, ann_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = var[:][:,:,:]

	ann_gcm = np.nanmean(np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1))
	sea_gcm = np.nanmean(np.nanmean(value[0:239:3,:,:], axis=1), axis=1)
	djf_gcm = np.nanmean(sea_gcm[3:80:4])
	mam_gcm = np.nanmean(sea_gcm[0:80:4])
	jja_gcm = np.nanmean(sea_gcm[1:80:4])
	son_gcm = np.nanmean(sea_gcm[2:80:4])

	return djf_gcm, mam_gcm, jja_gcm, son_gcm, ann_gcm
	

def comp_diff_rcp_hist_pre(rcp, hist):
	
	p1 = (rcp - hist) / hist
	p2 = p1 * 100
	
	return p2
	
	
# Import models 
# Precipitation 
djf_rcm_pre_hist_samz, mam_rcm_pre_hist_samz, jja_rcm_pre_hist_samz, son_rcm_pre_hist_samz, ann_rcm_pre_hist_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
djf_gcm_pre_hist_samz, mam_gcm_pre_hist_samz, jja_gcm_pre_hist_samz, son_gcm_pre_hist_samz, ann_gcm_pre_hist_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
djf_rcm_pre_rcp26_samz, mam_rcm_pre_rcp26_samz, jja_rcm_pre_rcp26_samz, son_rcm_pre_rcp26_samz, ann_rcm_pre_rcp26_samz = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
djf_gcm_pre_rcp26_samz, mam_gcm_pre_rcp26_samz, jja_gcm_pre_rcp26_samz, son_gcm_pre_rcp26_samz, ann_gcm_pre_rcp26_samz = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
djf_rcm_pre_rcp85_samz, mam_rcm_pre_rcp85_samz, jja_rcm_pre_rcp85_samz, son_rcm_pre_rcp85_samz, ann_rcm_pre_rcp85_samz = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
djf_gcm_pre_rcp85_samz, mam_gcm_pre_rcp85_samz, jja_gcm_pre_rcp85_samz, son_gcm_pre_rcp85_samz, ann_gcm_pre_rcp85_samz = import_gcm('pr', 'samz', 'rcp85', '2080-2099')

djf_rcm_pre_hist_eneb, mam_rcm_pre_hist_eneb, jja_rcm_pre_hist_eneb, son_rcm_pre_hist_eneb, ann_rcm_pre_hist_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
djf_gcm_pre_hist_eneb, mam_gcm_pre_hist_eneb, jja_gcm_pre_hist_eneb, son_gcm_pre_hist_eneb, ann_gcm_pre_hist_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
djf_rcm_pre_rcp26_eneb, mam_rcm_pre_rcp26_eneb, jja_rcm_pre_rcp26_eneb, son_rcm_pre_rcp26_eneb, ann_rcm_pre_rcp26_eneb = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
djf_gcm_pre_rcp26_eneb, mam_gcm_pre_rcp26_eneb, jja_gcm_pre_rcp26_eneb, son_gcm_pre_rcp26_eneb, ann_gcm_pre_rcp26_eneb = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
djf_rcm_pre_rcp85_eneb, mam_rcm_pre_rcp85_eneb, jja_rcm_pre_rcp85_eneb, son_rcm_pre_rcp85_eneb, ann_rcm_pre_rcp85_eneb = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
djf_gcm_pre_rcp85_eneb, mam_gcm_pre_rcp85_eneb, jja_gcm_pre_rcp85_eneb, son_gcm_pre_rcp85_eneb, ann_gcm_pre_rcp85_eneb = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')

djf_rcm_pre_hist_matopiba, mam_rcm_pre_hist_matopiba, jja_rcm_pre_hist_matopiba, son_rcm_pre_hist_matopiba, ann_rcm_pre_hist_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
djf_gcm_pre_hist_matopiba, mam_gcm_pre_hist_matopiba, jja_gcm_pre_hist_matopiba, son_gcm_pre_hist_matopiba, ann_gcm_pre_hist_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
djf_rcm_pre_rcp26_matopiba, mam_rcm_pre_rcp26_matopiba, jja_rcm_pre_rcp26_matopiba, son_rcm_pre_rcp26_matopiba, ann_rcm_pre_rcp26_matopiba = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
djf_gcm_pre_rcp26_matopiba, mam_gcm_pre_rcp26_matopiba, jja_gcm_pre_rcp26_matopiba, son_gcm_pre_rcp26_matopiba, ann_gcm_pre_rcp26_matopiba = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')
djf_rcm_pre_rcp85_matopiba, mam_rcm_pre_rcp85_matopiba, jja_rcm_pre_rcp85_matopiba, son_rcm_pre_rcp85_matopiba, ann_rcm_pre_rcp85_matopiba = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
djf_gcm_pre_rcp85_matopiba, mam_gcm_pre_rcp85_matopiba, jja_gcm_pre_rcp85_matopiba, son_gcm_pre_rcp85_matopiba, ann_gcm_pre_rcp85_matopiba = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')

# Temperature
djf_rcm_tas_hist_samz, mam_rcm_tas_hist_samz, jja_rcm_tas_hist_samz, son_rcm_tas_hist_samz, ann_rcm_tas_hist_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
djf_gcm_tas_hist_samz, mam_gcm_tas_hist_samz, jja_gcm_tas_hist_samz, son_gcm_tas_hist_samz, ann_gcm_tas_hist_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
djf_rcm_tas_rcp26_samz, mam_rcm_tas_rcp26_samz, jja_rcm_tas_rcp26_samz, son_rcm_tas_rcp26_samz, ann_rcm_tas_rcp26_samz = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
djf_gcm_tas_rcp26_samz, mam_gcm_tas_rcp26_samz, jja_gcm_tas_rcp26_samz, son_gcm_tas_rcp26_samz, ann_gcm_tas_rcp26_samz = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
djf_rcm_tas_rcp85_samz, mam_rcm_tas_rcp85_samz, jja_rcm_tas_rcp85_samz, son_rcm_tas_rcp85_samz, ann_rcm_tas_rcp85_samz = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
djf_gcm_tas_rcp85_samz, mam_gcm_tas_rcp85_samz, jja_gcm_tas_rcp85_samz, son_gcm_tas_rcp85_samz, ann_gcm_tas_rcp85_samz = import_gcm('tas', 'samz', 'rcp85', '2080-2099')

djf_rcm_tas_hist_eneb, mam_rcm_tas_hist_eneb, jja_rcm_tas_hist_eneb, son_rcm_tas_hist_eneb, ann_rcm_tas_hist_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
djf_gcm_tas_hist_eneb, mam_gcm_tas_hist_eneb, jja_gcm_tas_hist_eneb, son_gcm_tas_hist_eneb, ann_gcm_tas_hist_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
djf_rcm_tas_rcp26_eneb, mam_rcm_tas_rcp26_eneb, jja_rcm_tas_rcp26_eneb, son_rcm_tas_rcp26_eneb, ann_rcm_tas_rcp26_eneb = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
djf_gcm_tas_rcp26_eneb, mam_gcm_tas_rcp26_eneb, jja_gcm_tas_rcp26_eneb, son_gcm_tas_rcp26_eneb, ann_gcm_tas_rcp26_eneb = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
djf_rcm_tas_rcp85_eneb, mam_rcm_tas_rcp85_eneb, jja_rcm_tas_rcp85_eneb, son_rcm_tas_rcp85_eneb, ann_rcm_tas_rcp85_eneb = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
djf_gcm_tas_rcp85_eneb, mam_gcm_tas_rcp85_eneb, jja_gcm_tas_rcp85_eneb, son_gcm_tas_rcp85_eneb, ann_gcm_tas_rcp85_eneb = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')

djf_rcm_tas_hist_matopiba, mam_rcm_tas_hist_matopiba, jja_rcm_tas_hist_matopiba, son_rcm_tas_hist_matopiba, ann_rcm_tas_hist_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
djf_gcm_tas_hist_matopiba, mam_gcm_tas_hist_matopiba, jja_gcm_tas_hist_matopiba, son_gcm_tas_hist_matopiba, ann_gcm_tas_hist_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
djf_rcm_tas_rcp26_matopiba, mam_rcm_tas_rcp26_matopiba, jja_rcm_tas_rcp26_matopiba, son_rcm_tas_rcp26_matopiba, ann_rcm_tas_rcp26_matopiba = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
djf_gcm_tas_rcp26_matopiba, mam_gcm_tas_rcp26_matopiba, jja_gcm_tas_rcp26_matopiba, son_gcm_tas_rcp26_matopiba, ann_gcm_tas_rcp26_matopiba = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')
djf_rcm_tas_rcp85_matopiba, mam_rcm_tas_rcp85_matopiba, jja_rcm_tas_rcp85_matopiba, son_rcm_tas_rcp85_matopiba, ann_rcm_tas_rcp85_matopiba = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
djf_gcm_tas_rcp85_matopiba, mam_gcm_tas_rcp85_matopiba, jja_gcm_tas_rcp85_matopiba, son_gcm_tas_rcp85_matopiba, ann_gcm_tas_rcp85_matopiba = import_gcm('tas', 'matopiba', 'rcp85', '2080-2099')

# Compute skill metrics
# Precipitation - samz
diff_djf_rcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp26_samz, djf_rcm_pre_hist_samz)
diff_djf_gcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp26_samz, djf_gcm_pre_hist_samz)
diff_djf_rcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp85_samz, djf_rcm_pre_hist_samz)
diff_djf_gcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp85_samz, djf_gcm_pre_hist_samz)

diff_mam_rcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp26_samz, mam_rcm_pre_hist_samz)
diff_mam_gcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp26_samz, mam_gcm_pre_hist_samz)
diff_mam_rcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp85_samz, mam_rcm_pre_hist_samz)
diff_mam_gcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp85_samz, mam_gcm_pre_hist_samz)

diff_jja_rcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp26_samz, jja_rcm_pre_hist_samz)
diff_jja_gcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp26_samz, jja_gcm_pre_hist_samz)
diff_jja_rcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp85_samz, jja_rcm_pre_hist_samz)
diff_jja_gcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp85_samz, jja_gcm_pre_hist_samz)

diff_son_rcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(son_rcm_pre_rcp26_samz, son_rcm_pre_hist_samz)
diff_son_gcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(son_gcm_pre_rcp26_samz, son_gcm_pre_hist_samz)
diff_son_rcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(son_rcm_pre_rcp85_samz, son_rcm_pre_hist_samz)
diff_son_gcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(son_gcm_pre_rcp85_samz, son_gcm_pre_hist_samz)

diff_ann_rcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp26_samz, ann_rcm_pre_hist_samz)
diff_ann_gcm_pre_rcp26_samz = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp26_samz, ann_gcm_pre_hist_samz)
diff_ann_rcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp85_samz, ann_rcm_pre_hist_samz)
diff_ann_gcm_pre_rcp85_samz = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp85_samz, ann_gcm_pre_hist_samz)

# eneb
diff_djf_rcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp26_eneb, djf_rcm_pre_hist_eneb)
diff_djf_gcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp26_eneb, djf_gcm_pre_hist_eneb)
diff_djf_rcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp85_eneb, djf_rcm_pre_hist_eneb)
diff_djf_gcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp85_eneb, djf_gcm_pre_hist_eneb)

diff_mam_rcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp26_eneb, mam_rcm_pre_hist_eneb)
diff_mam_gcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp26_eneb, mam_gcm_pre_hist_eneb)
diff_mam_rcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp85_eneb, mam_rcm_pre_hist_eneb)
diff_mam_gcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp85_eneb, mam_gcm_pre_hist_eneb)

diff_jja_rcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp26_eneb, jja_rcm_pre_hist_eneb)
diff_jja_gcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp26_eneb, jja_gcm_pre_hist_eneb)
diff_jja_rcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp85_eneb, jja_rcm_pre_hist_eneb)
diff_jja_gcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp85_eneb, jja_gcm_pre_hist_eneb)

diff_son_rcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(son_rcm_pre_rcp26_eneb, son_rcm_pre_hist_eneb)
diff_son_gcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(son_gcm_pre_rcp26_eneb, son_gcm_pre_hist_eneb)
diff_son_rcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(son_rcm_pre_rcp85_eneb, son_rcm_pre_hist_eneb)
diff_son_gcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(son_gcm_pre_rcp85_eneb, son_gcm_pre_hist_eneb)

diff_ann_rcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp26_eneb, ann_rcm_pre_hist_eneb)
diff_ann_gcm_pre_rcp26_eneb = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp26_eneb, ann_gcm_pre_hist_eneb)
diff_ann_rcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp85_eneb, ann_rcm_pre_hist_eneb)
diff_ann_gcm_pre_rcp85_eneb = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp85_eneb, ann_gcm_pre_hist_eneb)

# matopiba
diff_djf_rcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp26_matopiba, djf_rcm_pre_hist_matopiba)
diff_djf_gcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp26_matopiba, djf_gcm_pre_hist_matopiba)
diff_djf_rcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(djf_rcm_pre_rcp85_matopiba, djf_rcm_pre_hist_matopiba)
diff_djf_gcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(djf_gcm_pre_rcp85_matopiba, djf_gcm_pre_hist_matopiba)

diff_mam_rcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp26_matopiba, mam_rcm_pre_hist_matopiba)
diff_mam_gcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp26_matopiba, mam_gcm_pre_hist_matopiba)
diff_mam_rcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(mam_rcm_pre_rcp85_matopiba, mam_rcm_pre_hist_matopiba)
diff_mam_gcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(mam_gcm_pre_rcp85_matopiba, mam_gcm_pre_hist_matopiba)

diff_jja_rcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp26_matopiba, jja_rcm_pre_hist_matopiba)
diff_jja_gcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp26_matopiba, jja_gcm_pre_hist_matopiba)
diff_jja_rcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(jja_rcm_pre_rcp85_matopiba, jja_rcm_pre_hist_matopiba)
diff_jja_gcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(jja_gcm_pre_rcp85_matopiba, jja_gcm_pre_hist_matopiba)

diff_son_rcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(son_rcm_pre_rcp26_matopiba, son_rcm_pre_hist_matopiba)
diff_son_gcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(son_gcm_pre_rcp26_matopiba, son_gcm_pre_hist_matopiba)
diff_son_rcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(son_rcm_pre_rcp85_matopiba, son_rcm_pre_hist_matopiba)
diff_son_gcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(son_gcm_pre_rcp85_matopiba, son_gcm_pre_hist_matopiba)

diff_ann_rcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp26_matopiba, ann_rcm_pre_hist_matopiba)
diff_ann_gcm_pre_rcp26_matopiba = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp26_matopiba, ann_gcm_pre_hist_matopiba)
diff_ann_rcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(ann_rcm_pre_rcp85_matopiba, ann_rcm_pre_hist_matopiba)
diff_ann_gcm_pre_rcp85_matopiba = comp_diff_rcp_hist_pre(ann_gcm_pre_rcp85_matopiba, ann_gcm_pre_hist_matopiba)

# Temperature - samz
diff_djf_rcm_tas_rcp26_samz = djf_rcm_tas_rcp26_samz - djf_rcm_tas_hist_samz
diff_djf_gcm_tas_rcp26_samz = djf_gcm_tas_rcp26_samz - djf_gcm_tas_hist_samz
diff_djf_rcm_tas_rcp85_samz = djf_rcm_tas_rcp85_samz - djf_rcm_tas_hist_samz
diff_djf_gcm_tas_rcp85_samz = djf_gcm_tas_rcp85_samz - djf_gcm_tas_hist_samz

diff_mam_rcm_tas_rcp26_samz = mam_rcm_tas_rcp26_samz - mam_rcm_tas_hist_samz
diff_mam_gcm_tas_rcp26_samz = mam_gcm_tas_rcp26_samz - mam_gcm_tas_hist_samz
diff_mam_rcm_tas_rcp85_samz = mam_rcm_tas_rcp85_samz - mam_rcm_tas_hist_samz
diff_mam_gcm_tas_rcp85_samz = mam_gcm_tas_rcp85_samz - mam_gcm_tas_hist_samz

diff_jja_rcm_tas_rcp26_samz = jja_rcm_tas_rcp26_samz - jja_rcm_tas_hist_samz
diff_jja_gcm_tas_rcp26_samz = jja_gcm_tas_rcp26_samz - jja_gcm_tas_hist_samz
diff_jja_rcm_tas_rcp85_samz = jja_rcm_tas_rcp85_samz - jja_rcm_tas_hist_samz
diff_jja_gcm_tas_rcp85_samz = jja_gcm_tas_rcp85_samz - jja_gcm_tas_hist_samz

diff_son_rcm_tas_rcp26_samz = son_rcm_tas_rcp26_samz - son_rcm_tas_hist_samz
diff_son_gcm_tas_rcp26_samz = son_gcm_tas_rcp26_samz - son_gcm_tas_hist_samz
diff_son_rcm_tas_rcp85_samz = son_rcm_tas_rcp85_samz - son_rcm_tas_hist_samz
diff_son_gcm_tas_rcp85_samz = son_gcm_tas_rcp85_samz - son_gcm_tas_hist_samz

diff_ann_rcm_tas_rcp26_samz = ann_rcm_tas_rcp26_samz - ann_rcm_tas_hist_samz
diff_ann_gcm_tas_rcp26_samz = ann_gcm_tas_rcp26_samz - ann_gcm_tas_hist_samz
diff_ann_rcm_tas_rcp85_samz = ann_rcm_tas_rcp85_samz - ann_rcm_tas_hist_samz
diff_ann_gcm_tas_rcp85_samz = ann_gcm_tas_rcp85_samz - ann_gcm_tas_hist_samz

# eneb
diff_djf_rcm_tas_rcp26_eneb = djf_rcm_tas_rcp26_eneb - djf_rcm_tas_hist_eneb
diff_djf_gcm_tas_rcp26_eneb = djf_gcm_tas_rcp26_eneb - djf_gcm_tas_hist_eneb
diff_djf_rcm_tas_rcp85_eneb = djf_rcm_tas_rcp85_eneb - djf_rcm_tas_hist_eneb
diff_djf_gcm_tas_rcp85_eneb = djf_gcm_tas_rcp85_eneb - djf_gcm_tas_hist_eneb

diff_mam_rcm_tas_rcp26_eneb = mam_rcm_tas_rcp26_eneb - mam_rcm_tas_hist_eneb
diff_mam_gcm_tas_rcp26_eneb = mam_gcm_tas_rcp26_eneb - mam_gcm_tas_hist_eneb
diff_mam_rcm_tas_rcp85_eneb = mam_rcm_tas_rcp85_eneb - mam_rcm_tas_hist_eneb
diff_mam_gcm_tas_rcp85_eneb = mam_gcm_tas_rcp85_eneb - mam_gcm_tas_hist_eneb

diff_jja_rcm_tas_rcp26_eneb = jja_rcm_tas_rcp26_eneb - jja_rcm_tas_hist_eneb
diff_jja_gcm_tas_rcp26_eneb = jja_gcm_tas_rcp26_eneb - jja_gcm_tas_hist_eneb
diff_jja_rcm_tas_rcp85_eneb = jja_rcm_tas_rcp85_eneb - jja_rcm_tas_hist_eneb
diff_jja_gcm_tas_rcp85_eneb = jja_gcm_tas_rcp85_eneb - jja_gcm_tas_hist_eneb

diff_son_rcm_tas_rcp26_eneb = son_rcm_tas_rcp26_eneb - son_rcm_tas_hist_eneb
diff_son_gcm_tas_rcp26_eneb = son_gcm_tas_rcp26_eneb - son_gcm_tas_hist_eneb
diff_son_rcm_tas_rcp85_eneb = son_rcm_tas_rcp85_eneb - son_rcm_tas_hist_eneb
diff_son_gcm_tas_rcp85_eneb = son_gcm_tas_rcp85_eneb - son_gcm_tas_hist_eneb

diff_ann_rcm_tas_rcp26_eneb = ann_rcm_tas_rcp26_eneb - ann_rcm_tas_hist_eneb
diff_ann_gcm_tas_rcp26_eneb = ann_gcm_tas_rcp26_eneb - ann_gcm_tas_hist_eneb
diff_ann_rcm_tas_rcp85_eneb = ann_rcm_tas_rcp85_eneb - ann_rcm_tas_hist_eneb
diff_ann_gcm_tas_rcp85_eneb = ann_gcm_tas_rcp85_eneb - ann_gcm_tas_hist_eneb

# matopiba
diff_djf_rcm_tas_rcp26_matopiba = djf_rcm_tas_rcp26_matopiba - djf_rcm_tas_hist_matopiba
diff_djf_gcm_tas_rcp26_matopiba = djf_gcm_tas_rcp26_matopiba - djf_gcm_tas_hist_matopiba
diff_djf_rcm_tas_rcp85_matopiba = djf_rcm_tas_rcp85_matopiba - djf_rcm_tas_hist_matopiba
diff_djf_gcm_tas_rcp85_matopiba = djf_gcm_tas_rcp85_matopiba - djf_gcm_tas_hist_matopiba

diff_mam_rcm_tas_rcp26_matopiba = mam_rcm_tas_rcp26_matopiba - mam_rcm_tas_hist_matopiba
diff_mam_gcm_tas_rcp26_matopiba = mam_gcm_tas_rcp26_matopiba - mam_gcm_tas_hist_matopiba
diff_mam_rcm_tas_rcp85_matopiba = mam_rcm_tas_rcp85_matopiba - mam_rcm_tas_hist_matopiba
diff_mam_gcm_tas_rcp85_matopiba = mam_gcm_tas_rcp85_matopiba - mam_gcm_tas_hist_matopiba

diff_jja_rcm_tas_rcp26_matopiba = jja_rcm_tas_rcp26_matopiba - jja_rcm_tas_hist_matopiba
diff_jja_gcm_tas_rcp26_matopiba = jja_gcm_tas_rcp26_matopiba - jja_gcm_tas_hist_matopiba
diff_jja_rcm_tas_rcp85_matopiba = jja_rcm_tas_rcp85_matopiba - jja_rcm_tas_hist_matopiba
diff_jja_gcm_tas_rcp85_matopiba = jja_gcm_tas_rcp85_matopiba - jja_gcm_tas_hist_matopiba

diff_son_rcm_tas_rcp26_matopiba = son_rcm_tas_rcp26_matopiba - son_rcm_tas_hist_matopiba
diff_son_gcm_tas_rcp26_matopiba = son_gcm_tas_rcp26_matopiba - son_gcm_tas_hist_matopiba
diff_son_rcm_tas_rcp85_matopiba = son_rcm_tas_rcp85_matopiba - son_rcm_tas_hist_matopiba
diff_son_gcm_tas_rcp85_matopiba = son_gcm_tas_rcp85_matopiba - son_gcm_tas_hist_matopiba

diff_ann_rcm_tas_rcp26_matopiba = ann_rcm_tas_rcp26_matopiba - ann_rcm_tas_hist_matopiba
diff_ann_gcm_tas_rcp26_matopiba = ann_gcm_tas_rcp26_matopiba - ann_gcm_tas_hist_matopiba
diff_ann_rcm_tas_rcp85_matopiba = ann_rcm_tas_rcp85_matopiba - ann_rcm_tas_hist_matopiba
diff_ann_gcm_tas_rcp85_matopiba = ann_gcm_tas_rcp85_matopiba - ann_gcm_tas_hist_matopiba

# Plot models
fig = plt.figure()
time = np.arange(1, 6)

ax = fig.add_subplot(3, 2, 1)
ax.plot(1, diff_djf_rcm_pre_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_samz, mfc='red', color='black', marker='s')
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 2)
ax.plot(1, diff_djf_rcm_tas_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_samz, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_samz, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_samz, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_samz, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_samz, mfc='red', color='black', marker='s')
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 3)
ax.plot(1, diff_djf_rcm_pre_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_eneb, mfc='red', color='black', marker='s')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation change (%)', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 4)
ax.plot(1, diff_djf_rcm_tas_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_eneb, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_eneb, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_eneb, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_eneb, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_eneb, mfc='red', color='black', marker='s')
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature change (°C)', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 5)
ax.plot(1, diff_djf_rcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_matopiba, mfc='red', color='black', marker='s')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 6)
ax.plot(1, diff_djf_rcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_matopiba, mfc='red', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_matopiba, mfc='blue', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_matopiba, mfc='red', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_matopiba, mfc='red', color='black', marker='s')
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_diff_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
