# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script scatter plot from Reg and Had models"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties


def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_rcm = np.nanmean(np.nanmean(value[0:239:12,:,:], axis=1), axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[0:239:3,:,:], axis=1), axis=1)
	djf_rcm = np.nanmean(season_rcm[3:80:4])
	mam_rcm = np.nanmean(season_rcm[0:80:4])
	jja_rcm = np.nanmean(season_rcm[1:80:4])
	son_rcm = np.nanmean(season_rcm[2:80:4])

	return annual_rcm, djf_rcm, mam_rcm, jja_rcm, son_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = var[:][:,:,:]

	annual_gcm = np.nanmean(np.nanmean(value[0:239:12,:,:], axis=1), axis=1) 
	season_gcm = np.nanmean(np.nanmean(value[0:239:3,:,:], axis=1), axis=1)
	djf_gcm = np.nanmean(season_gcm[3:80:4])
	mam_gcm = np.nanmean(season_gcm[0:80:4])
	jja_gcm = np.nanmean(season_gcm[1:80:4])
	son_gcm = np.nanmean(season_gcm[2:80:4])

	return annual_gcm, djf_gcm, mam_gcm, jja_gcm, son_gcm
	

# Import regcm and hadgem models
# Precipitation 
annual_rcm_pre_hist_samz, djf_rcm_pre_hist_samz, mam_rcm_pre_hist_samz, jja_rcm_pre_hist_samz, son_rcm_pre_hist_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
annual_gcm_pre_hist_samz, djf_gcm_pre_hist_samz, mam_gcm_pre_hist_samz, jja_gcm_pre_hist_samz, son_gcm_pre_hist_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
annual_rcm_pre_rcp26_samz, djf_rcm_pre_rcp26_samz, mam_rcm_pre_rcp26_samz, jja_rcm_pre_rcp26_samz, son_rcm_pre_rcp26_samz = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
annual_gcm_pre_rcp26_samz, djf_gcm_pre_rcp26_samz, mam_gcm_pre_rcp26_samz, jja_gcm_pre_rcp26_samz, son_gcm_pre_rcp26_samz = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
annual_rcm_pre_rcp85_samz, djf_rcm_pre_rcp85_samz, mam_rcm_pre_rcp85_samz, jja_rcm_pre_rcp85_samz, son_rcm_pre_rcp85_samz = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
annual_gcm_pre_rcp85_samz, djf_gcm_pre_rcp85_samz, mam_gcm_pre_rcp85_samz, jja_gcm_pre_rcp85_samz, son_gcm_pre_rcp85_samz = import_gcm('pr', 'samz', 'rcp85', '2080-2099')

annual_rcm_pre_hist_eneb, djf_rcm_pre_hist_eneb, mam_rcm_pre_hist_eneb, jja_rcm_pre_hist_eneb, son_rcm_pre_hist_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
annual_gcm_pre_hist_eneb, djf_gcm_pre_hist_eneb, mam_gcm_pre_hist_eneb, jja_gcm_pre_hist_eneb, son_gcm_pre_hist_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
annual_rcm_pre_rcp26_eneb, djf_rcm_pre_rcp26_eneb, mam_rcm_pre_rcp26_eneb, jja_rcm_pre_rcp26_eneb, son_rcm_pre_rcp26_eneb = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
annual_gcm_pre_rcp26_eneb, djf_gcm_pre_rcp26_eneb, mam_gcm_pre_rcp26_eneb, jja_gcm_pre_rcp26_eneb, son_gcm_pre_rcp26_eneb = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
annual_rcm_pre_rcp85_eneb, djf_rcm_pre_rcp85_eneb, mam_rcm_pre_rcp85_eneb, jja_rcm_pre_rcp85_eneb, son_rcm_pre_rcp85_eneb = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
annual_gcm_pre_rcp85_eneb, djf_gcm_pre_rcp85_eneb, mam_gcm_pre_rcp85_eneb, jja_gcm_pre_rcp85_eneb, son_gcm_pre_rcp85_eneb = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')

annual_rcm_pre_hist_matopiba, djf_rcm_pre_hist_matopiba, mam_rcm_pre_hist_matopiba, jja_rcm_pre_hist_matopiba, son_rcm_pre_hist_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
annual_gcm_pre_hist_matopiba, djf_gcm_pre_hist_matopiba, mam_gcm_pre_hist_matopiba, jja_gcm_pre_hist_matopiba, son_gcm_pre_hist_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
annual_rcm_pre_rcp26_matopiba, djf_rcm_pre_rcp26_matopiba, mam_rcm_pre_rcp26_matopiba, jja_rcm_pre_rcp26_matopiba, son_rcm_pre_rcp26_matopiba = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
annual_gcm_pre_rcp26_matopiba, djf_gcm_pre_rcp26_matopiba, mam_gcm_pre_rcp26_matopiba, jja_gcm_pre_rcp26_matopiba, son_gcm_pre_rcp26_matopiba = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')
annual_rcm_pre_rcp85_matopiba, djf_rcm_pre_rcp85_matopiba, mam_rcm_pre_rcp85_matopiba, jja_rcm_pre_rcp85_matopiba, son_rcm_pre_rcp85_matopiba = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
annual_gcm_pre_rcp85_matopiba, djf_gcm_pre_rcp85_matopiba, mam_gcm_pre_rcp85_matopiba, jja_gcm_pre_rcp85_matopiba, son_gcm_pre_rcp85_matopiba = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')

# Temperature
annual_rcm_tas_hist_samz, djf_rcm_tas_hist_samz, mam_rcm_tas_hist_samz, jja_rcm_tas_hist_samz, son_rcm_tas_hist_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
annual_gcm_tas_hist_samz, djf_gcm_tas_hist_samz, mam_gcm_tas_hist_samz, jja_gcm_tas_hist_samz, son_gcm_tas_hist_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
annual_rcm_tas_rcp26_samz, djf_rcm_tas_rcp26_samz, mam_rcm_tas_rcp26_samz, jja_rcm_tas_rcp26_samz, son_rcm_tas_rcp26_samz = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
annual_gcm_tas_rcp26_samz, djf_gcm_tas_rcp26_samz, mam_gcm_tas_rcp26_samz, jja_gcm_tas_rcp26_samz, son_gcm_tas_rcp26_samz = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
annual_rcm_tas_rcp85_samz, djf_rcm_tas_rcp85_samz, mam_rcm_tas_rcp85_samz, jja_rcm_tas_rcp85_samz, son_rcm_tas_rcp85_samz = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
annual_gcm_tas_rcp85_samz, djf_gcm_tas_rcp85_samz, mam_gcm_tas_rcp85_samz, jja_gcm_tas_rcp85_samz, son_gcm_tas_rcp85_samz = import_gcm('tas', 'samz', 'rcp85', '2080-2099')

annual_rcm_tas_hist_eneb, djf_rcm_tas_hist_eneb, mam_rcm_tas_hist_eneb, jja_rcm_tas_hist_eneb, son_rcm_tas_hist_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
annual_gcm_tas_hist_eneb, djf_gcm_tas_hist_eneb, mam_gcm_tas_hist_eneb, jja_gcm_tas_hist_eneb, son_gcm_tas_hist_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
annual_rcm_tas_rcp26_eneb, djf_rcm_tas_rcp26_eneb, mam_rcm_tas_rcp26_eneb, jja_rcm_tas_rcp26_eneb, son_rcm_tas_rcp26_eneb = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
annual_gcm_tas_rcp26_eneb, djf_gcm_tas_rcp26_eneb, mam_gcm_tas_rcp26_eneb, jja_gcm_tas_rcp26_eneb, son_gcm_tas_rcp26_eneb = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
annual_rcm_tas_rcp85_eneb, djf_rcm_tas_rcp85_eneb, mam_rcm_tas_rcp85_eneb, jja_rcm_tas_rcp85_eneb, son_rcm_tas_rcp85_eneb = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
annual_gcm_tas_rcp85_eneb, djf_gcm_tas_rcp85_eneb, mam_gcm_tas_rcp85_eneb, jja_gcm_tas_rcp85_eneb, son_gcm_tas_rcp85_eneb = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')

annual_rcm_tas_hist_matopiba, djf_rcm_tas_hist_matopiba, mam_rcm_tas_hist_matopiba, jja_rcm_tas_hist_matopiba, son_rcm_tas_hist_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
annual_gcm_tas_hist_matopiba, djf_gcm_tas_hist_matopiba, mam_gcm_tas_hist_matopiba, jja_gcm_tas_hist_matopiba, son_gcm_tas_hist_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
annual_rcm_tas_rcp26_matopiba, djf_rcm_tas_rcp26_matopiba, mam_rcm_tas_rcp26_matopiba, jja_rcm_tas_rcp26_matopiba, son_rcm_tas_rcp26_matopiba = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
annual_gcm_tas_rcp26_matopiba, djf_gcm_tas_rcp26_matopiba, mam_gcm_tas_rcp26_matopiba, jja_gcm_tas_rcp26_matopiba, son_gcm_tas_rcp26_matopiba = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')
annual_rcm_tas_rcp85_matopiba, djf_rcm_tas_rcp85_matopiba, mam_rcm_tas_rcp85_matopiba, jja_rcm_tas_rcp85_matopiba, son_rcm_tas_rcp85_matopiba = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
annual_gcm_tas_rcp85_matopiba, djf_gcm_tas_rcp85_matopiba, mam_gcm_tas_rcp85_matopiba, jja_gcm_tas_rcp85_matopiba, son_gcm_tas_rcp85_matopiba = import_gcm('tas', 'matopiba', 'rcp85', '2080-2099')

# Compute skill metrics
# Precipitation - samz
diff_annual_rcm_pre_rcp26_hist_samz = annual_rcm_pre_rcp26_samz - annual_rcm_pre_hist_samz
diff_annual_gcm_pre_rcp26_hist_samz = annual_gcm_pre_rcp26_samz - annual_gcm_pre_hist_samz
diff_annual_rcm_pre_rcp85_hist_samz = annual_rcm_pre_rcp85_samz - annual_rcm_pre_hist_samz
diff_annual_gcm_pre_rcp85_hist_samz = annual_gcm_pre_rcp85_samz - annual_gcm_pre_hist_samz

diff_djf_rcm_pre_rcp26_hist_samz = djf_rcm_pre_rcp26_samz - djf_rcm_pre_hist_samz
diff_djf_gcm_pre_rcp26_hist_samz = djf_gcm_pre_rcp26_samz - djf_gcm_pre_hist_samz
diff_djf_rcm_pre_rcp85_hist_samz = djf_rcm_pre_rcp85_samz - djf_rcm_pre_hist_samz
diff_djf_gcm_pre_rcp85_hist_samz = djf_gcm_pre_rcp85_samz - djf_gcm_pre_hist_samz

diff_mam_rcm_pre_rcp26_hist_samz = mam_rcm_pre_rcp26_samz - mam_rcm_pre_hist_samz
diff_mam_gcm_pre_rcp26_hist_samz = mam_gcm_pre_rcp26_samz - mam_gcm_pre_hist_samz
diff_mam_rcm_pre_rcp85_hist_samz = mam_rcm_pre_rcp85_samz - mam_rcm_pre_hist_samz
diff_mam_gcm_pre_rcp85_hist_samz = mam_gcm_pre_rcp85_samz - mam_gcm_pre_hist_samz

diff_jja_rcm_pre_rcp26_hist_samz = jja_rcm_pre_rcp26_samz - jja_rcm_pre_hist_samz
diff_jja_gcm_pre_rcp26_hist_samz = jja_gcm_pre_rcp26_samz - jja_gcm_pre_hist_samz
diff_jja_rcm_pre_rcp85_hist_samz = jja_rcm_pre_rcp85_samz - jja_rcm_pre_hist_samz
diff_jja_gcm_pre_rcp85_hist_samz = jja_gcm_pre_rcp85_samz - jja_gcm_pre_hist_samz

diff_son_rcm_pre_rcp26_hist_samz = son_rcm_pre_rcp26_samz - son_rcm_pre_hist_samz
diff_son_gcm_pre_rcp26_hist_samz = son_gcm_pre_rcp26_samz - son_gcm_pre_hist_samz
diff_son_rcm_pre_rcp85_hist_samz = son_rcm_pre_rcp85_samz - son_rcm_pre_hist_samz
diff_son_gcm_pre_rcp85_hist_samz = son_gcm_pre_rcp85_samz - son_gcm_pre_hist_samz

# eneb
diff_annual_rcm_pre_rcp26_hist_eneb = annual_rcm_pre_rcp26_eneb - annual_rcm_pre_hist_eneb
diff_annual_gcm_pre_rcp26_hist_eneb = annual_gcm_pre_rcp26_eneb - annual_gcm_pre_hist_eneb
diff_annual_rcm_pre_rcp85_hist_eneb = annual_rcm_pre_rcp85_eneb - annual_rcm_pre_hist_eneb
diff_annual_gcm_pre_rcp85_hist_eneb = annual_gcm_pre_rcp85_eneb - annual_gcm_pre_hist_eneb

diff_djf_rcm_pre_rcp26_hist_eneb = djf_rcm_pre_rcp26_eneb - djf_rcm_pre_hist_eneb
diff_djf_gcm_pre_rcp26_hist_eneb = djf_gcm_pre_rcp26_eneb - djf_gcm_pre_hist_eneb
diff_djf_rcm_pre_rcp85_hist_eneb = djf_rcm_pre_rcp85_eneb - djf_rcm_pre_hist_eneb
diff_djf_gcm_pre_rcp85_hist_eneb = djf_gcm_pre_rcp85_eneb - djf_gcm_pre_hist_eneb

diff_mam_rcm_pre_rcp26_hist_eneb = mam_rcm_pre_rcp26_eneb - mam_rcm_pre_hist_eneb
diff_mam_gcm_pre_rcp26_hist_eneb = mam_gcm_pre_rcp26_eneb - mam_gcm_pre_hist_eneb
diff_mam_rcm_pre_rcp85_hist_eneb = mam_rcm_pre_rcp85_eneb - mam_rcm_pre_hist_eneb
diff_mam_gcm_pre_rcp85_hist_eneb = mam_gcm_pre_rcp85_eneb - mam_gcm_pre_hist_eneb

diff_jja_rcm_pre_rcp26_hist_eneb = jja_rcm_pre_rcp26_eneb - jja_rcm_pre_hist_eneb
diff_jja_gcm_pre_rcp26_hist_eneb = jja_gcm_pre_rcp26_eneb - jja_gcm_pre_hist_eneb
diff_jja_rcm_pre_rcp85_hist_eneb = jja_rcm_pre_rcp85_eneb - jja_rcm_pre_hist_eneb
diff_jja_gcm_pre_rcp85_hist_eneb = jja_gcm_pre_rcp85_eneb - jja_gcm_pre_hist_eneb

diff_son_rcm_pre_rcp26_hist_eneb = son_rcm_pre_rcp26_eneb - son_rcm_pre_hist_eneb
diff_son_gcm_pre_rcp26_hist_eneb = son_gcm_pre_rcp26_eneb - son_gcm_pre_hist_eneb
diff_son_rcm_pre_rcp85_hist_eneb = son_rcm_pre_rcp85_eneb - son_rcm_pre_hist_eneb
diff_son_gcm_pre_rcp85_hist_eneb = son_gcm_pre_rcp85_eneb - son_gcm_pre_hist_eneb

# matopiba
diff_annual_rcm_pre_rcp26_hist_matopiba = annual_rcm_pre_rcp26_matopiba - annual_rcm_pre_hist_matopiba
diff_annual_gcm_pre_rcp26_hist_matopiba = annual_gcm_pre_rcp26_matopiba - annual_gcm_pre_hist_matopiba
diff_annual_rcm_pre_rcp85_hist_matopiba = annual_rcm_pre_rcp85_matopiba - annual_rcm_pre_hist_matopiba
diff_annual_gcm_pre_rcp85_hist_matopiba = annual_gcm_pre_rcp85_matopiba - annual_gcm_pre_hist_matopiba

diff_djf_rcm_pre_rcp26_hist_matopiba = djf_rcm_pre_rcp26_matopiba - djf_rcm_pre_hist_matopiba
diff_djf_gcm_pre_rcp26_hist_matopiba = djf_gcm_pre_rcp26_matopiba - djf_gcm_pre_hist_matopiba
diff_djf_rcm_pre_rcp85_hist_matopiba = djf_rcm_pre_rcp85_matopiba - djf_rcm_pre_hist_matopiba
diff_djf_gcm_pre_rcp85_hist_matopiba = djf_gcm_pre_rcp85_matopiba - djf_gcm_pre_hist_matopiba

diff_mam_rcm_pre_rcp26_hist_matopiba = mam_rcm_pre_rcp26_matopiba - mam_rcm_pre_hist_matopiba
diff_mam_gcm_pre_rcp26_hist_matopiba = mam_gcm_pre_rcp26_matopiba - mam_gcm_pre_hist_matopiba
diff_mam_rcm_pre_rcp85_hist_matopiba = mam_rcm_pre_rcp85_matopiba - mam_rcm_pre_hist_matopiba
diff_mam_gcm_pre_rcp85_hist_matopiba = mam_gcm_pre_rcp85_matopiba - mam_gcm_pre_hist_matopiba

diff_jja_rcm_pre_rcp26_hist_matopiba = jja_rcm_pre_rcp26_matopiba - jja_rcm_pre_hist_matopiba
diff_jja_gcm_pre_rcp26_hist_matopiba = jja_gcm_pre_rcp26_matopiba - jja_gcm_pre_hist_matopiba
diff_jja_rcm_pre_rcp85_hist_matopiba = jja_rcm_pre_rcp85_matopiba - jja_rcm_pre_hist_matopiba
diff_jja_gcm_pre_rcp85_hist_matopiba = jja_gcm_pre_rcp85_matopiba - jja_gcm_pre_hist_matopiba

diff_son_rcm_pre_rcp26_hist_matopiba = son_rcm_pre_rcp26_matopiba - son_rcm_pre_hist_matopiba
diff_son_gcm_pre_rcp26_hist_matopiba = son_gcm_pre_rcp26_matopiba - son_gcm_pre_hist_matopiba
diff_son_rcm_pre_rcp85_hist_matopiba = son_rcm_pre_rcp85_matopiba - son_rcm_pre_hist_matopiba
diff_son_gcm_pre_rcp85_hist_matopiba = son_gcm_pre_rcp85_matopiba - son_gcm_pre_hist_matopiba

# Temperature - samz
diff_annual_rcm_tas_rcp26_hist_samz = annual_rcm_tas_rcp26_samz - annual_rcm_tas_hist_samz
diff_annual_gcm_tas_rcp26_hist_samz = annual_gcm_tas_rcp26_samz - annual_gcm_tas_hist_samz
diff_annual_rcm_tas_rcp85_hist_samz = annual_rcm_tas_rcp85_samz - annual_rcm_tas_hist_samz
diff_annual_gcm_tas_rcp85_hist_samz = annual_gcm_tas_rcp85_samz - annual_gcm_tas_hist_samz

diff_djf_rcm_tas_rcp26_hist_samz = djf_rcm_tas_rcp26_samz - djf_rcm_tas_hist_samz
diff_djf_gcm_tas_rcp26_hist_samz = djf_gcm_tas_rcp26_samz - djf_gcm_tas_hist_samz
diff_djf_rcm_tas_rcp85_hist_samz = djf_rcm_tas_rcp85_samz - djf_rcm_tas_hist_samz
diff_djf_gcm_tas_rcp85_hist_samz = djf_gcm_tas_rcp85_samz - djf_gcm_tas_hist_samz

diff_mam_rcm_tas_rcp26_hist_samz = mam_rcm_tas_rcp26_samz - mam_rcm_tas_hist_samz
diff_mam_gcm_tas_rcp26_hist_samz = mam_gcm_tas_rcp26_samz - mam_gcm_tas_hist_samz
diff_mam_rcm_tas_rcp85_hist_samz = mam_rcm_tas_rcp85_samz - mam_rcm_tas_hist_samz
diff_mam_gcm_tas_rcp85_hist_samz = mam_gcm_tas_rcp85_samz - mam_gcm_tas_hist_samz

diff_jja_rcm_tas_rcp26_hist_samz = jja_rcm_tas_rcp26_samz - jja_rcm_tas_hist_samz
diff_jja_gcm_tas_rcp26_hist_samz = jja_gcm_tas_rcp26_samz - jja_gcm_tas_hist_samz
diff_jja_rcm_tas_rcp85_hist_samz = jja_rcm_tas_rcp85_samz - jja_rcm_tas_hist_samz
diff_jja_gcm_tas_rcp85_hist_samz = jja_gcm_tas_rcp85_samz - jja_gcm_tas_hist_samz

diff_son_rcm_tas_rcp26_hist_samz = son_rcm_tas_rcp26_samz - son_rcm_tas_hist_samz
diff_son_gcm_tas_rcp26_hist_samz = son_gcm_tas_rcp26_samz - son_gcm_tas_hist_samz
diff_son_rcm_tas_rcp85_hist_samz = son_rcm_tas_rcp85_samz - son_rcm_tas_hist_samz
diff_son_gcm_tas_rcp85_hist_samz = son_gcm_tas_rcp85_samz - son_gcm_tas_hist_samz

# eneb
diff_annual_rcm_tas_rcp26_hist_eneb = annual_rcm_tas_rcp26_eneb - annual_rcm_tas_hist_eneb
diff_annual_gcm_tas_rcp26_hist_eneb = annual_gcm_tas_rcp26_eneb - annual_gcm_tas_hist_eneb
diff_annual_rcm_tas_rcp85_hist_eneb = annual_rcm_tas_rcp85_eneb - annual_rcm_tas_hist_eneb
diff_annual_gcm_tas_rcp85_hist_eneb = annual_gcm_tas_rcp85_eneb - annual_gcm_tas_hist_eneb

diff_djf_rcm_tas_rcp26_hist_eneb = djf_rcm_tas_rcp26_eneb - djf_rcm_tas_hist_eneb
diff_djf_gcm_tas_rcp26_hist_eneb = djf_gcm_tas_rcp26_eneb - djf_gcm_tas_hist_eneb
diff_djf_rcm_tas_rcp85_hist_eneb = djf_rcm_tas_rcp85_eneb - djf_rcm_tas_hist_eneb
diff_djf_gcm_tas_rcp85_hist_eneb = djf_gcm_tas_rcp85_eneb - djf_gcm_tas_hist_eneb

diff_mam_rcm_tas_rcp26_hist_eneb = mam_rcm_tas_rcp26_eneb - mam_rcm_tas_hist_eneb
diff_mam_gcm_tas_rcp26_hist_eneb = mam_gcm_tas_rcp26_eneb - mam_gcm_tas_hist_eneb
diff_mam_rcm_tas_rcp85_hist_eneb = mam_rcm_tas_rcp85_eneb - mam_rcm_tas_hist_eneb
diff_mam_gcm_tas_rcp85_hist_eneb = mam_gcm_tas_rcp85_eneb - mam_gcm_tas_hist_eneb

diff_jja_rcm_tas_rcp26_hist_eneb = jja_rcm_tas_rcp26_eneb - jja_rcm_tas_hist_eneb
diff_jja_gcm_tas_rcp26_hist_eneb = jja_gcm_tas_rcp26_eneb - jja_gcm_tas_hist_eneb
diff_jja_rcm_tas_rcp85_hist_eneb = jja_rcm_tas_rcp85_eneb - jja_rcm_tas_hist_eneb
diff_jja_gcm_tas_rcp85_hist_eneb = jja_gcm_tas_rcp85_eneb - jja_gcm_tas_hist_eneb

diff_son_rcm_tas_rcp26_hist_eneb = son_rcm_tas_rcp26_eneb - son_rcm_tas_hist_eneb
diff_son_gcm_tas_rcp26_hist_eneb = son_gcm_tas_rcp26_eneb - son_gcm_tas_hist_eneb
diff_son_rcm_tas_rcp85_hist_eneb = son_rcm_tas_rcp85_eneb - son_rcm_tas_hist_eneb
diff_son_gcm_tas_rcp85_hist_eneb = son_gcm_tas_rcp85_eneb - son_gcm_tas_hist_eneb

# matopiba
diff_annual_rcm_tas_rcp26_hist_matopiba = annual_rcm_tas_rcp26_matopiba - annual_rcm_tas_hist_matopiba
diff_annual_gcm_tas_rcp26_hist_matopiba = annual_gcm_tas_rcp26_matopiba - annual_gcm_tas_hist_matopiba
diff_annual_rcm_tas_rcp85_hist_matopiba = annual_rcm_tas_rcp85_matopiba - annual_rcm_tas_hist_matopiba
diff_annual_gcm_tas_rcp85_hist_matopiba = annual_gcm_tas_rcp85_matopiba - annual_gcm_tas_hist_matopiba

diff_djf_rcm_tas_rcp26_hist_matopiba = djf_rcm_tas_rcp26_matopiba - djf_rcm_tas_hist_matopiba
diff_djf_gcm_tas_rcp26_hist_matopiba = djf_gcm_tas_rcp26_matopiba - djf_gcm_tas_hist_matopiba
diff_djf_rcm_tas_rcp85_hist_matopiba = djf_rcm_tas_rcp85_matopiba - djf_rcm_tas_hist_matopiba
diff_djf_gcm_tas_rcp85_hist_matopiba = djf_gcm_tas_rcp85_matopiba - djf_gcm_tas_hist_matopiba

diff_mam_rcm_tas_rcp26_hist_matopiba = mam_rcm_tas_rcp26_matopiba - mam_rcm_tas_hist_matopiba
diff_mam_gcm_tas_rcp26_hist_matopiba = mam_gcm_tas_rcp26_matopiba - mam_gcm_tas_hist_matopiba
diff_mam_rcm_tas_rcp85_hist_matopiba = mam_rcm_tas_rcp85_matopiba - mam_rcm_tas_hist_matopiba
diff_mam_gcm_tas_rcp85_hist_matopiba = mam_gcm_tas_rcp85_matopiba - mam_gcm_tas_hist_matopiba

diff_jja_rcm_tas_rcp26_hist_matopiba = jja_rcm_tas_rcp26_matopiba - jja_rcm_tas_hist_matopiba
diff_jja_gcm_tas_rcp26_hist_matopiba = jja_gcm_tas_rcp26_matopiba - jja_gcm_tas_hist_matopiba
diff_jja_rcm_tas_rcp85_hist_matopiba = jja_rcm_tas_rcp85_matopiba - jja_rcm_tas_hist_matopiba
diff_jja_gcm_tas_rcp85_hist_matopiba = jja_gcm_tas_rcp85_matopiba - jja_gcm_tas_hist_matopiba

diff_son_rcm_tas_rcp26_hist_matopiba = son_rcm_tas_rcp26_matopiba - son_rcm_tas_hist_matopiba
diff_son_gcm_tas_rcp26_hist_matopiba = son_gcm_tas_rcp26_matopiba - son_gcm_tas_hist_matopiba
diff_son_rcm_tas_rcp85_hist_matopiba = son_rcm_tas_rcp85_matopiba - son_rcm_tas_hist_matopiba
diff_son_gcm_tas_rcp85_hist_matopiba = son_gcm_tas_rcp85_matopiba - son_gcm_tas_hist_matopiba

print('Precipitation Annual DJF MAM JJA SON: RCM HAD: SAMZ')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_samz, axis=0))

print('Precipitation Annual DJF MAM JJA SON: RCM HAD: ENEB')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_eneb, axis=0))

print('Precipitation Annual DJF MAM JJA SON: RCM HAD: MATOPIBA')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_matopiba, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM HAD: SAMZ')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_samz, axis=0), axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_samz, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM HAD: ENEB')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_eneb, axis=0), axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_eneb, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM HAD: MATOPIBA')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_matopiba, axis=0), axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_matopiba, axis=0))
exit()

# Compute skill metrics
fig=plt.figure()

ax1 = plt.subplot(3, 5, 1)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax1.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax1.text(-2.0,17.3, '   Annual  ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
             
ax2 = plt.subplot(3, 5, 2)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax2.set_xticklabels([])
ax2.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax2.text(-2.0,17.3, '     DJF     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
                
ax3 = plt.subplot(3, 5, 3)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax3.set_xticklabels([])
ax3.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax3.text(-2.0,17.3, '    MAM    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})

ax4 = plt.subplot(3, 5, 4)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax4.set_xticklabels([])
ax4.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax4.text(-2.0,17.3, '      JJA     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
   
ax5 = plt.subplot(3, 5, 5)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax5.set_xticklabels([])
ax5.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax5.text(-2.0,17.3, '     SON    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
title = ax5.text(3.75,11., '     SAMZ      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax6 = plt.subplot(3, 5, 6)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax6.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

ax7 = plt.subplot(3, 5, 7)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax7.set_xticklabels([])
ax7.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax8 = plt.subplot(3, 5, 8)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax8.set_xticklabels([])
ax8.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax9 = plt.subplot(3, 5, 9)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-1, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax9.set_xticklabels([])
ax9.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax10 = plt.subplot(3, 5, 10)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax10.set_xticklabels([])
ax10.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax10.text(3.75,11., '     ENEB      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax11 = plt.subplot(3, 5, 11)
plt.scatter(3, 10, s=80, color='blue', marker='o', label='Reg (RCP26)')
plt.scatter(3, 14, s=80, color='blue', marker='D', label='Had (RCP26)')
plt.scatter(-2, -10, s=80, color='red', marker='o', label='Reg (RCP85)')
plt.scatter(1, -5, s=80, color='red', marker='D', label='Had (RCP85)')
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax11.legend(loc='lower left', bbox_to_anchor=(-0.7, -0.7), shadow=True, ncol=4)

ax12 = plt.subplot(3, 5, 12)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax12.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax13 = plt.subplot(3, 5, 13)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax13.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.xlabel(u'Temperature (°C)', fontweight='bold')

ax14 = plt.subplot(3, 5, 14)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax14.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax15 = plt.subplot(3, 5, 15)
plt.scatter(3, 10, s=80, color='blue', marker='o')
plt.scatter(3, 14, s=80, color='blue', marker='D')
plt.scatter(-2, -10, s=80, color='red', marker='o')
plt.scatter(1, -5, s=80, color='red', marker='D')
ax15.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax15.text(3.75,11., ' MATOPIBA  ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.6, 'pad':4})
                
# Path out to save bias figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_scatter_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


