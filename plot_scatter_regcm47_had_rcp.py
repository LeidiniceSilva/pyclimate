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

print('Precipitation Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: SAMZ')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_samz, axis=0))

print('Precipitation Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: ENEB')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_eneb, axis=0))

print('Precipitation Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: MATOPIBA')
print(np.nanmean(diff_annual_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_annual_gcm_pre_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_djf_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_pre_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_mam_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_pre_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_jja_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_pre_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_son_rcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_pre_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_pre_rcp85_hist_matopiba, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: SAMZ')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_samz, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_samz, axis=0), axis=0))
print()
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_samz, axis=0))
print()
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_samz, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_samz, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_samz, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: ENEB')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_eneb, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_eneb, axis=0), axis=0))
print()
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_eneb, axis=0))
print()
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_eneb, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_eneb, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_eneb, axis=0))

print('Temperature Annual DJF MAM JJA SON: RCM and GCM: RCP26 and RCP85: MATOPIBA')
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp26_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp26_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_rcm_tas_rcp85_hist_matopiba, axis=0), axis=0))
print(np.nanmean(np.nanmean(diff_annual_gcm_tas_rcp85_hist_matopiba, axis=0), axis=0))
print()
print(np.nanmean(diff_djf_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_djf_gcm_tas_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_mam_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_mam_gcm_tas_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_jja_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_jja_gcm_tas_rcp85_hist_matopiba, axis=0))
print()
print(np.nanmean(diff_son_rcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp26_hist_matopiba, axis=0))
print(np.nanmean(diff_son_rcm_tas_rcp85_hist_matopiba, axis=0))
print(np.nanmean(diff_son_gcm_tas_rcp85_hist_matopiba, axis=0))

# Plot scatter plot with diff between rcp and historical
fig=plt.figure()

ax1 = plt.subplot(3, 5, 1)
plt.scatter(-0.05, 1.56, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.47, 1.19, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-1.51, 7.44, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-1.88, 6.30, s=80, color='red', marker='s', edgecolor='black')
ax1.set_xlim(-2, 2)
ax1.set_ylim(-10, 10)
ax1.set_xticks(np.arange(-2, 3, 1))
ax1.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax1.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax1.text(-1.5, -9, 'A)', fontweight='bold')
title = ax1.text(-2, 11.5, '   Annual   ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
             
ax2 = plt.subplot(3, 5, 2)
plt.scatter(-0.57, 1.78, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(0.50, 1.89, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-1.21, 6.18, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(2.04, 9.41, s=80, color='red', marker='s', edgecolor='black')
ax2.set_xlim(-2, 2)
ax2.set_ylim(-10, 10)
ax2.set_xticks(np.arange(-2, 3, 1))
ax2.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax2.text(-1.5, -9, 'B)', fontweight='bold')
title = ax2.text(-2, 11.5, '      DJF      ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
                
ax3 = plt.subplot(3, 5, 3)
plt.scatter(-0.05, 1.56, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.47, 1.19, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-1.51, 7.44, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-1.88, 6.30, s=80, color='red', marker='s', edgecolor='black')
ax3.set_xlim(-2, 2)
ax3.set_ylim(-10, 10)
ax3.set_xticks(np.arange(-2, 3, 1))
ax3.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax3.text(-1.5, -9, 'C)', fontweight='bold')
title = ax3.text(-2, 11.5, '     MAM    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})

ax4 = plt.subplot(3, 5, 4)
plt.scatter(0.28, 1.15, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.25, 1.15, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(0.01, 4.98, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-0.48, 6.05, s=80, color='red', marker='s', edgecolor='black')
ax4.set_xlim(-2, 2)
ax4.set_ylim(-10, 10)
ax4.set_xticks(np.arange(-2, 3, 1))
ax4.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax4.text(-1.5, -9, 'D)', fontweight='bold')
title = ax4.text(-2, 11.5, '       JJA      ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
   
ax5 = plt.subplot(3, 5, 5)
plt.scatter(0.04, 1.70, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.26, 1.61, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(0.16, 5.79, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(0.32, 7.02, s=80, c='red', marker='s', edgecolor='black')
ax5.set_xlim(-2, 2)
ax5.set_ylim(-10, 10)
ax5.set_xticks(np.arange(-2, 3, 1))
ax5.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax5.text(-1.5, -9, 'E)', fontweight='bold')
title = ax5.text(-2, 11.5, '     SON     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
title = ax5.text(2.3, 6.7, '     SAMZ      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax6 = plt.subplot(3, 5, 6)
plt.scatter(-0.32, 1.78, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.29, 1.47, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(0.23, 4.76, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(1.05, 4.30, s=80, color='red', marker='s', edgecolor='black')
ax6.set_xlim(-2, 2)
ax6.set_ylim(-10, 10)
ax6.set_xticks(np.arange(-2, 3, 1))
ax6.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax6.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax6.text(-1.5, -9, 'F)', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')

ax7 = plt.subplot(3, 5, 7)
plt.scatter(-0.02, 1.32, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.06, 1.15, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-0.06, 5.18, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-0.10, 4.28, s=80, color='red', marker='s', edgecolor='black')
ax7.set_xlim(-2, 2)
ax7.set_ylim(-10, 10)
ax7.set_xticks(np.arange(-2, 3, 1))
ax7.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax7.get_xticklabels(), visible=False)
plt.setp(ax7.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax7.text(-1.5, -9, 'G)', fontweight='bold')

ax8 = plt.subplot(3, 5, 8)
plt.scatter(-0.32, 1.78, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.29, 1.47, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(0.23, 4.76, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(1.05, 4.30, s=80, color='red', marker='s', edgecolor='black')
ax8.set_xlim(-2, 2)
ax8.set_ylim(-10, 10)
ax8.set_xticks(np.arange(-2, 3, 1))
ax8.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax8.get_xticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax8.text(-1.5, -9, 'H)', fontweight='bold')

ax9 = plt.subplot(3, 5, 9)
plt.scatter(-0.42, 1.35, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.67, 1.30, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-0.38, 4.53, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-0.15, 4.04, s=80, color='red', marker='s', edgecolor='black')
ax9.set_xlim(-2, 2)
ax9.set_ylim(-10, 10)
ax9.set_xticks(np.arange(-2, 3, 1))
ax9.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax9.get_xticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax9.text(-1.5, -9, 'I)', fontweight='bold')

ax10 = plt.subplot(3, 5, 10)
plt.scatter(-0.14, 1.68, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.07, 1.51, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-0.22, 5.27, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-0.14, 4.51, s=80, color='red', marker='s', edgecolor='black')
ax10.set_xlim(-2, 2)
ax10.set_ylim(-10, 10)
ax10.set_xticks(np.arange(-2, 3, 1))
ax10.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax10.get_xticklabels(), visible=False)
plt.setp(ax10.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax10.text(-1.5, -9, 'J)', fontweight='bold')
title = ax10.text(2.3, 6.7, '      ENEB      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax11 = plt.subplot(3, 5, 11)
plt.scatter(-0.32, 1.70, s=80, color='blue', marker='o', label='Reg (RCP26)', edgecolor='black')
plt.scatter(0.28, 1.27, s=80, color='blue', marker='s', label='Had (RCP26)', edgecolor='black')
plt.scatter(-0.69, 5.44, s=80, color='red', marker='o', label='Reg (RCP85)', edgecolor='black')
plt.scatter(0.76, 4.78, s=80, color='red', marker='s', label='Had (RCP85)', edgecolor='black')
ax11.set_xlim(-2, 2)
ax11.set_ylim(-10, 10)
ax11.set_xticks(np.arange(-2, 3, 1))
ax11.set_yticks(np.arange(-10, 15, 5))
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax11.text(-1.5, -9, 'K)', fontweight='bold')
ax11.legend(loc='lower left', bbox_to_anchor=(-0.7, -0.7), shadow=True, ncol=4)

ax12 = plt.subplot(3, 5, 12)
plt.scatter(-0.12, 1.71, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.59, 1.99, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-1.15, 6.31, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-1.57, 7.30, s=80, color='red', marker='s', edgecolor='black')
ax12.set_xlim(-2, 2)
ax12.set_ylim(-10, 10)
ax12.set_xticks(np.arange(-2, 3, 1))
ax12.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax12.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax12.text(-1.5, -9, 'L)', fontweight='bold')

ax13 = plt.subplot(3, 5, 13)
plt.scatter(-.32, 1.70, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(0.28, 1.27, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-0.69, 5.44, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(0.76, 4.78, s=80, color='red', marker='s', edgecolor='black')
ax13.set_xlim(-2, 2)
ax13.set_ylim(-10, 10)
ax13.set_xticks(np.arange(-2, 3, 1))
ax13.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax13.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax13.text(-1.5, -9, 'M)', fontweight='bold')
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

ax14 = plt.subplot(3, 5, 14)
plt.scatter(-0.19, 1.36, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(-0.82, 1.16, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-1.15, 5.56, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-1.31, 5.06, s=80, color='red', marker='s', edgecolor='black')
ax14.set_xlim(-2, 2)
ax14.set_ylim(-10, 10)
ax14.set_xticks(np.arange(-2, 3, 1))
ax14.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax14.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax14.text(-1.5, -9, 'N)', fontweight='bold')

ax15 = plt.subplot(3, 5, 15)
plt.scatter(-0.02, 1.82, s=80, color='blue', marker='o', edgecolor='black')
plt.scatter(0.02, 1.79, s=80, color='blue', marker='s', edgecolor='black')
plt.scatter(-0.05, 6.30, s=80, color='red', marker='o', edgecolor='black')
plt.scatter(-0.23, 6.56, s=80, color='red', marker='o', edgecolor='black')
ax15.set_xlim(-2, 2)
ax15.set_ylim(-10, 10)
ax15.set_xticks(np.arange(-2, 3, 1))
ax15.set_yticks(np.arange(-10, 15, 5))
plt.setp(ax15.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax15.text(-1.5, -9, 'O)', fontweight='bold')
title = ax15.text(2.3, 6.7, ' MATOPIBA  ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})

fig.tight_layout()
#~ plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)
            
# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_scatter_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


