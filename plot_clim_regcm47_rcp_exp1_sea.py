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

	jan = np.nanmean(np.nanmean(np.nanmean(value[0::12,:,:], axis=1), axis=1))
	feb = np.nanmean(np.nanmean(np.nanmean(value[1::12,:,:], axis=1), axis=1))
	mar = np.nanmean(np.nanmean(np.nanmean(value[2::12,:,:], axis=1), axis=1))
	apr = np.nanmean(np.nanmean(np.nanmean(value[3::12,:,:], axis=1), axis=1))
	may = np.nanmean(np.nanmean(np.nanmean(value[4::12,:,:], axis=1), axis=1))
	jun = np.nanmean(np.nanmean(np.nanmean(value[5::12,:,:], axis=1), axis=1))
	jul = np.nanmean(np.nanmean(np.nanmean(value[6::12,:,:], axis=1), axis=1))
	aug = np.nanmean(np.nanmean(np.nanmean(value[7::12,:,:], axis=1), axis=1))
	sep = np.nanmean(np.nanmean(np.nanmean(value[8::12,:,:], axis=1), axis=1))
	out = np.nanmean(np.nanmean(np.nanmean(value[9::12,:,:], axis=1), axis=1))
	nov = np.nanmean(np.nanmean(np.nanmean(value[10::12,:,:], axis=1), axis=1))
	dec = np.nanmean(np.nanmean(np.nanmean(value[11::12,:,:], axis=1), axis=1))

	djf = (dec+jan+feb)/3
	mam = (mar+apr+may)/3
	jja = (jun+jul+aug)/3
	son = (sep+out+nov)/3

	a1 = np.nanmean(np.nanmean(np.nanmean(value[0:12,:,:], axis=1), axis=1))
	a2 = np.nanmean(np.nanmean(np.nanmean(value[12:24,:,:], axis=1), axis=1))
	a3 = np.nanmean(np.nanmean(np.nanmean(value[24:36,:,:], axis=1), axis=1))
	a4 = np.nanmean(np.nanmean(np.nanmean(value[36:48,:,:], axis=1), axis=1))
	a5 = np.nanmean(np.nanmean(np.nanmean(value[48:60,:,:], axis=1), axis=1))
	a6 = np.nanmean(np.nanmean(np.nanmean(value[60:72,:,:], axis=1), axis=1))
	a7 = np.nanmean(np.nanmean(np.nanmean(value[72:84,:,:], axis=1), axis=1))
	a8 = np.nanmean(np.nanmean(np.nanmean(value[84:106,:,:], axis=1), axis=1))
	a9 = np.nanmean(np.nanmean(np.nanmean(value[106:118,:,:], axis=1), axis=1))
	a10 = np.nanmean(np.nanmean(np.nanmean(value[118:120,:,:], axis=1), axis=1))
	a11 = np.nanmean(np.nanmean(np.nanmean(value[120:132,:,:], axis=1), axis=1))
	a12 = np.nanmean(np.nanmean(np.nanmean(value[132:144,:,:], axis=1), axis=1))
	a13 = np.nanmean(np.nanmean(np.nanmean(value[144:156,:,:], axis=1), axis=1))
	a14 = np.nanmean(np.nanmean(np.nanmean(value[156:168,:,:], axis=1), axis=1))
	a15 = np.nanmean(np.nanmean(np.nanmean(value[168:180,:,:], axis=1), axis=1))
	a16 = np.nanmean(np.nanmean(np.nanmean(value[180:192,:,:], axis=1), axis=1))
	a17 = np.nanmean(np.nanmean(np.nanmean(value[192:204,:,:], axis=1), axis=1))
	a18 = np.nanmean(np.nanmean(np.nanmean(value[204:216,:,:], axis=1), axis=1))
	a19 = np.nanmean(np.nanmean(np.nanmean(value[216:228,:,:], axis=1), axis=1))
	a20 = np.nanmean(np.nanmean(np.nanmean(value[228:240,:,:], axis=1), axis=1))

	ann = (a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20)/20
	
	return djf, mam, jja, son, ann


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = var[:][:,:,:]

	jan = np.nanmean(np.nanmean(np.nanmean(value[0::12,:,:], axis=1), axis=1))
	feb = np.nanmean(np.nanmean(np.nanmean(value[1::12,:,:], axis=1), axis=1))
	mar = np.nanmean(np.nanmean(np.nanmean(value[2::12,:,:], axis=1), axis=1))
	apr = np.nanmean(np.nanmean(np.nanmean(value[3::12,:,:], axis=1), axis=1))
	may = np.nanmean(np.nanmean(np.nanmean(value[4::12,:,:], axis=1), axis=1))
	jun = np.nanmean(np.nanmean(np.nanmean(value[5::12,:,:], axis=1), axis=1))
	jul = np.nanmean(np.nanmean(np.nanmean(value[6::12,:,:], axis=1), axis=1))
	aug = np.nanmean(np.nanmean(np.nanmean(value[7::12,:,:], axis=1), axis=1))
	sep = np.nanmean(np.nanmean(np.nanmean(value[8::12,:,:], axis=1), axis=1))
	out = np.nanmean(np.nanmean(np.nanmean(value[9::12,:,:], axis=1), axis=1))
	nov = np.nanmean(np.nanmean(np.nanmean(value[10::12,:,:], axis=1), axis=1))
	dec = np.nanmean(np.nanmean(np.nanmean(value[11::12,:,:], axis=1), axis=1))

	djf = (dec+jan+feb)/3
	mam = (mar+apr+may)/3
	jja = (jun+jul+aug)/3
	son = (sep+out+nov)/3

	a1 = np.nanmean(np.nanmean(np.nanmean(value[0:12,:,:], axis=1), axis=1))
	a2 = np.nanmean(np.nanmean(np.nanmean(value[12:24,:,:], axis=1), axis=1))
	a3 = np.nanmean(np.nanmean(np.nanmean(value[24:36,:,:], axis=1), axis=1))
	a4 = np.nanmean(np.nanmean(np.nanmean(value[36:48,:,:], axis=1), axis=1))
	a5 = np.nanmean(np.nanmean(np.nanmean(value[48:60,:,:], axis=1), axis=1))
	a6 = np.nanmean(np.nanmean(np.nanmean(value[60:72,:,:], axis=1), axis=1))
	a7 = np.nanmean(np.nanmean(np.nanmean(value[72:84,:,:], axis=1), axis=1))
	a8 = np.nanmean(np.nanmean(np.nanmean(value[84:106,:,:], axis=1), axis=1))
	a9 = np.nanmean(np.nanmean(np.nanmean(value[106:118,:,:], axis=1), axis=1))
	a10 = np.nanmean(np.nanmean(np.nanmean(value[118:120,:,:], axis=1), axis=1))
	a11 = np.nanmean(np.nanmean(np.nanmean(value[120:132,:,:], axis=1), axis=1))
	a12 = np.nanmean(np.nanmean(np.nanmean(value[132:144,:,:], axis=1), axis=1))
	a13 = np.nanmean(np.nanmean(np.nanmean(value[144:156,:,:], axis=1), axis=1))
	a14 = np.nanmean(np.nanmean(np.nanmean(value[156:168,:,:], axis=1), axis=1))
	a15 = np.nanmean(np.nanmean(np.nanmean(value[168:180,:,:], axis=1), axis=1))
	a16 = np.nanmean(np.nanmean(np.nanmean(value[180:192,:,:], axis=1), axis=1))
	a17 = np.nanmean(np.nanmean(np.nanmean(value[192:204,:,:], axis=1), axis=1))
	a18 = np.nanmean(np.nanmean(np.nanmean(value[204:216,:,:], axis=1), axis=1))
	a19 = np.nanmean(np.nanmean(np.nanmean(value[216:228,:,:], axis=1), axis=1))
	a20 = np.nanmean(np.nanmean(np.nanmean(value[228:240,:,:], axis=1), axis=1))

	ann = (a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17+a18+a19+a20)/20
	
	return djf, mam, jja, son, ann
	

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
ax.plot(1, diff_djf_rcm_pre_rcp26_samz, label='RegCM4.7 RCP2.6', mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_samz, label='RegCM4.7 RCP8.5', mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_samz, label='HadGEM2-ES RCP2.6', mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_samz, label='HadGEM2-ES RCP8.5',mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_samz, mfc='dimgray', color='black', marker='s')
plt.title(u'A) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
plt.legend(loc=9, handlelength=0.50, handleheight=0.50, shadow=True, ncol=2, fontsize=5.5)

ax = fig.add_subplot(3, 2, 2)
ax.plot(1, diff_djf_rcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_samz, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_samz, mfc='dimgray', color='black', marker='s')
plt.title(u'D) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 3)
ax.plot(1, diff_djf_rcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_eneb, mfc='dimgray', color='black', marker='s')
plt.title(u'B) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation change (%)', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 4)
ax.plot(1, diff_djf_rcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_eneb, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_eneb, mfc='dimgray', color='black', marker='s')
plt.title(u'E) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature change (°C)', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 5)
ax.plot(1, diff_djf_rcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_pre_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
plt.title(u'C) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Period', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.ylim(-80, 80)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 6)
ax.plot(1, diff_djf_rcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='o')
ax.plot(1, diff_djf_rcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(2, diff_mam_rcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(3, diff_jja_rcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(4, diff_son_rcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(5, diff_ann_rcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='o')
ax.plot(1, diff_djf_gcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp26_matopiba, mfc='lightgray', color='black', marker='s')
ax.plot(1, diff_djf_gcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(2, diff_mam_gcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(3, diff_jja_gcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(4, diff_son_gcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
ax.plot(5, diff_ann_gcm_tas_rcp85_matopiba, mfc='dimgray', color='black', marker='s')
plt.title(u'F) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Period', fontsize=8)
plt.xticks(time, ('DJF', 'MAM', 'JJA', 'SON', 'ANN'), fontsize=8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.ylim(-8, 8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_sea_diff_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
