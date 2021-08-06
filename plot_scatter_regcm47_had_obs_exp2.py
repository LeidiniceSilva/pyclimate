# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot scatter plot from Reg and Had models and obs database"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value
	
	
def import_rcm_exp1(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value


def import_rcm_exp2(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value
	
	
def compute_linear_regression(obs, sim):

	m, c, r, p, se1 = stats.linregress(obs, sim)
	eq = 'y = %2.2fx+%2.2f'%(m, c)
	r2 = 'R² = %1.2f'%(r)
   
	return m, c, eq, r2

       
# Import regcm exps model end obs database climatology
# Precipitation
mon_pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
mon_pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
mon_pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_pre_reg_exp1_samz = import_rcm_exp1('pr', 'samz', 'hist', '1986-2005')
mon_pre_reg_exp1_eneb = import_rcm_exp1('pr', 'eneb', 'hist', '1986-2005')
mon_pre_reg_exp1_matopiba = import_rcm_exp1('pr', 'matopiba', 'hist', '1986-2005')

mon_pre_reg_exp2_samz = import_rcm_exp2('pr', 'samz', 'hist', '1986-2005')
mon_pre_reg_exp2_eneb = import_rcm_exp2('pr', 'eneb', 'hist', '1986-2005')
mon_pre_reg_exp2_matopiba = import_rcm_exp2('pr', 'matopiba', 'hist', '1986-2005')

mon_pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
mon_tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
mon_tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
mon_tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_tas_reg_exp1_samz = import_rcm_exp1('tas', 'samz', 'hist', '1986-2005')
mon_tas_reg_exp1_eneb = import_rcm_exp1('tas', 'eneb', 'hist', '1986-2005')
mon_tas_reg_exp1_matopiba = import_rcm_exp1('tas', 'matopiba', 'hist', '1986-2005')

mon_tas_reg_exp2_samz = import_rcm_exp2('tas', 'samz', 'hist', '1986-2005')
mon_tas_reg_exp2_eneb = import_rcm_exp2('tas', 'eneb', 'hist', '1986-2005')
mon_tas_reg_exp2_matopiba = import_rcm_exp2('tas', 'matopiba', 'hist', '1986-2005')

mon_tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Imp linear regression
# Precipitation
m_pre_reg1_samz, c_pre_reg1_samz, eq_pre_reg1_samz, r2_pre_reg1_samz = compute_linear_regression(mon_pre_cru_samz, mon_pre_reg_exp1_samz)
m_pre_reg2_samz, c_pre_reg2_samz, eq_pre_reg2_samz, r2_pre_reg2_samz = compute_linear_regression(mon_pre_cru_samz, mon_pre_reg_exp2_samz)
m_pre_had_samz, c_pre_had_samz, eq_pre_had_samz, r2_pre_had_samz = compute_linear_regression(mon_pre_cru_samz, mon_pre_had_samz)

m_pre_reg1_eneb, c_pre_reg1_eneb, eq_pre_reg1_eneb, r2_pre_reg1_eneb = compute_linear_regression(mon_pre_cru_eneb, mon_pre_reg_exp1_eneb)
m_pre_reg2_eneb, c_pre_reg2_eneb, eq_pre_reg2_eneb, r2_pre_reg2_eneb = compute_linear_regression(mon_pre_cru_eneb, mon_pre_reg_exp2_eneb)
m_pre_had_eneb, c_pre_had_eneb, eq_pre_had_eneb, r2_pre_had_eneb = compute_linear_regression(mon_pre_cru_eneb, mon_pre_had_eneb)

m_pre_reg1_matopiba, c_pre_reg1_matopiba, eq_pre_reg1_matopiba, r2_pre_reg1_matopiba = compute_linear_regression(mon_pre_cru_matopiba, mon_pre_reg_exp1_matopiba)
m_pre_reg2_matopiba, c_pre_reg2_matopiba, eq_pre_reg2_matopiba, r2_pre_reg2_matopiba = compute_linear_regression(mon_pre_cru_matopiba, mon_pre_reg_exp2_matopiba)
m_pre_had_matopiba, c_pre_had_matopiba, eq_pre_had_matopiba, r2_pre_had_matopiba = compute_linear_regression(mon_pre_cru_matopiba, mon_pre_had_matopiba)

# Temperature
m_tas_reg1_samz, c_tas_reg1_samz, eq_tas_reg1_samz, r2_tas_reg1_samz = compute_linear_regression(mon_tas_cru_samz, np.nanmean(mon_tas_reg_exp1_samz, axis=1))
m_tas_reg2_samz, c_tas_reg2_samz, eq_tas_reg2_samz, r2_tas_reg2_samz = compute_linear_regression(mon_tas_cru_samz, np.nanmean(mon_tas_reg_exp2_samz, axis=1))
m_tas_had_samz, c_tas_had_samz, eq_tas_had_samz, r2_tas_had_samz = compute_linear_regression(mon_tas_cru_samz, mon_tas_had_samz)

m_tas_reg1_eneb, c_tas_reg1_eneb, eq_tas_reg1_eneb, r2_tas_reg1_eneb = compute_linear_regression(mon_tas_cru_eneb, np.nanmean(mon_tas_reg_exp1_eneb, axis=1))
m_tas_reg2_eneb, c_tas_reg2_eneb, eq_tas_reg2_eneb, r2_tas_reg2_eneb = compute_linear_regression(mon_tas_cru_eneb, np.nanmean(mon_tas_reg_exp2_eneb, axis=1))
m_tas_had_eneb, c_tas_had_eneb, eq_tas_had_eneb, r2_tas_had_eneb = compute_linear_regression(mon_tas_cru_eneb, mon_tas_had_eneb)

m_tas_reg1_matopiba, c_tas_reg1_matopiba, eq_tas_reg1_matopiba, r2_tas_reg1_matopiba = compute_linear_regression(mon_tas_cru_matopiba, np.nanmean(mon_tas_reg_exp1_matopiba, axis=1))
m_tas_reg2_matopiba, c_tas_reg2_matopiba, eq_tas_reg2_matopiba, r2_tas_reg2_matopiba = compute_linear_regression(mon_tas_cru_matopiba, np.nanmean(mon_tas_reg_exp2_matopiba, axis=1))
m_tas_had_matopiba, c_tas_had_matopiba, eq_tas_had_matopiba, r2_tas_had_matopiba = compute_linear_regression(mon_tas_cru_matopiba, mon_tas_had_matopiba)

# Plot model end obs data climatology
fig = plt.figure()
    
ax1 = fig.add_subplot(3, 2, 1)
plt.plot([0,15],[0,15], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_pre_cru_samz, mon_pre_reg_exp1_samz, edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_pre_cru_samz, mon_pre_reg_exp2_samz, edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_pre_cru_samz, mon_pre_had_samz, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_pre_cru_samz, m_pre_reg1_samz*mon_pre_cru_samz+c_pre_reg1_samz, color='blue', linewidth=1.)
#~ plt.plot(mon_pre_cru_samz, m_pre_reg2_samz*mon_pre_cru_samz+c_pre_reg2_samz, color='red', linewidth=1.)
#~ plt.plot(mon_pre_cru_samz, m_pre_had_samz*mon_pre_cru_samz+c_pre_had_samz, color='gray', linewidth=1.)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(0, 15)
plt.xticks(np.arange(0, 18, 3), fontsize=8)
plt.ylim(0, 15)
plt.yticks(np.arange(0, 18, 3), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.text(1, 16, r2_pre_reg1_samz, fontsize=6, color='blue')
plt.text(5.5, 16, r2_pre_reg2_samz, fontsize=6, color='red')
plt.text(10.5, 16, r2_pre_had_samz, fontsize=6, color='gray')
#~ plt.text(1, 16, eq_pre_reg1_samz, fontsize=6, color='blue')
#~ plt.text(5.5, 16, eq_pre_reg2_samz, fontsize=6, color='red')
#~ plt.text(10.5, 16, eq_pre_had_samz, fontsize=6, color='gray')

plt.legend(fontsize=6, loc=2, ncol=1, frameon=False)

ax2 = fig.add_subplot(3, 2, 2)
plt.plot([22,34],[22,34], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_tas_cru_samz, np.nanmean(mon_tas_reg_exp1_samz, axis=1), edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_tas_cru_samz, np.nanmean(mon_tas_reg_exp2_samz, axis=1), edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_tas_cru_samz, mon_tas_had_samz, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_tas_cru_samz, m_tas_reg1_samz*mon_tas_cru_samz+c_tas_reg1_samz, color='blue', linewidth=1.)
#~ plt.plot(mon_tas_cru_samz, m_tas_reg2_samz*mon_tas_cru_samz+c_tas_reg2_samz, color='red', linewidth=1.)
#~ plt.plot(mon_tas_cru_samz, m_tas_had_samz*mon_tas_cru_samz+c_tas_had_samz, color='gray', linewidth=1.)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(22, 34)
plt.xticks(np.arange(22, 36, 2), fontsize=8)
plt.ylim(22, 34)
plt.yticks(np.arange(22, 36, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.text(23, 35, r2_tas_reg1_samz, fontsize=6, color='blue')
plt.text(26, 35, r2_tas_reg2_samz, fontsize=6, color='red')
plt.text(30, 35, r2_tas_had_samz, fontsize=6, color='gray')
#~ plt.text(23, 35, eq_tas_reg1_samz, fontsize=6, color='blue')
#~ plt.text(26, 35, eq_tas_reg2_samz, fontsize=6, color='red')
#~ plt.text(30, 35, eq_tas_had_samz, fontsize=6, color='gray')

ax3 = fig.add_subplot(3, 2, 3)
plt.plot([0,15],[0,15], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_pre_cru_eneb, mon_pre_reg_exp1_eneb, edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_pre_cru_eneb, mon_pre_reg_exp2_eneb, edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_pre_cru_eneb, mon_pre_had_eneb, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_pre_cru_eneb, m_pre_reg1_eneb*mon_pre_cru_eneb+c_pre_reg1_eneb, color='blue', linewidth=1.)
#~ plt.plot(mon_pre_cru_eneb, m_pre_reg2_eneb*mon_pre_cru_eneb+c_pre_reg2_eneb, color='red', linewidth=1.)
#~ plt.plot(mon_pre_cru_eneb, m_pre_had_eneb*mon_pre_cru_eneb+c_pre_had_eneb, color='gray', linewidth=1.)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Simulated Precipitation (mm d⁻¹)', fontsize=8)
plt.xlim(0, 15)
plt.xticks(np.arange(0, 18, 3), fontsize=8)
plt.ylim(0, 15)
plt.yticks(np.arange(0, 18, 3), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.text(1, 16, r2_pre_reg1_eneb, fontsize=6, color='blue')
plt.text(5.5, 16, r2_pre_reg2_eneb, fontsize=6, color='red')
plt.text(10.5, 16, r2_pre_had_eneb, fontsize=6, color='gray')
#~ plt.text(1, 16, eq_pre_reg1_eneb, fontsize=6, color='blue')
#~ plt.text(5.5, 16, eq_pre_reg2_eneb, fontsize=6, color='red')
#~ plt.text(10.5, 16, eq_pre_had_eneb, fontsize=6, color='gray')

ax4 = fig.add_subplot(3, 2, 4)
plt.plot([22,34],[22,34], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_tas_cru_eneb, np.nanmean(mon_tas_reg_exp1_eneb, axis=1), edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_tas_cru_eneb, np.nanmean(mon_tas_reg_exp2_eneb, axis=1), edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_tas_cru_eneb, mon_tas_had_eneb, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_tas_cru_eneb, m_tas_reg1_eneb*mon_tas_cru_eneb+c_tas_reg1_eneb, color='blue', linewidth=1.)
#~ plt.plot(mon_tas_cru_eneb, m_tas_reg2_eneb*mon_tas_cru_eneb+c_tas_reg2_eneb, color='red', linewidth=1.)
#~ plt.plot(mon_tas_cru_eneb, m_tas_had_eneb*mon_tas_cru_eneb+c_tas_had_eneb, color='gray', linewidth=1.)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Simulated Temperature (°C)', fontsize=8)
plt.xlim(22, 34)
plt.xticks(np.arange(22, 36, 2), fontsize=8)
plt.ylim(22, 34)
plt.yticks(np.arange(22, 36, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.text(23, 35, r2_tas_reg1_eneb, fontsize=6, color='blue')
plt.text(26, 35, r2_tas_reg2_eneb, fontsize=6, color='red')
plt.text(30, 35, r2_tas_had_eneb, fontsize=6, color='gray')
#~ plt.text(23, 35, eq_tas_reg1_eneb, fontsize=6, color='blue')
#~ plt.text(26, 35, eq_tas_reg2_eneb, fontsize=6, color='red')
#~ plt.text(30, 35, eq_tas_had_eneb, fontsize=6, color='gray')

ax5 = fig.add_subplot(3, 2, 5)
plt.plot([0,15],[0,15], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_pre_cru_matopiba, mon_pre_reg_exp1_matopiba, edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_pre_cru_matopiba, mon_pre_reg_exp2_matopiba, edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_pre_cru_matopiba, mon_pre_had_matopiba, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_pre_cru_matopiba, m_pre_reg1_matopiba*mon_pre_cru_matopiba+c_pre_reg1_matopiba, color='blue', linewidth=1.)
#~ plt.plot(mon_pre_cru_matopiba, m_pre_reg2_matopiba*mon_pre_cru_matopiba+c_pre_reg2_matopiba, color='red', linewidth=1.)
#~ plt.plot(mon_pre_cru_matopiba, m_pre_had_matopiba*mon_pre_cru_matopiba+c_pre_had_matopiba, color='gray', linewidth=1.)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Observed Precipitation (mm d⁻¹)', fontsize=8)
plt.xlim(0, 15)
plt.xticks(np.arange(0, 18, 3), fontsize=8)
plt.ylim(0, 15)
plt.yticks(np.arange(0, 18, 3), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.text(1, 16, r2_pre_reg1_matopiba, fontsize=6, color='blue')
plt.text(5.5, 16, r2_pre_reg2_matopiba, fontsize=6, color='red')
plt.text(10.5, 16, r2_pre_had_matopiba, fontsize=6, color='gray')
#~ plt.text(1, 16, eq_pre_reg1_matopiba, fontsize=6, color='blue')
#~ plt.text(5.5, 16, eq_pre_reg2_matopiba, fontsize=6, color='red')
#~ plt.text(10.5, 16, eq_pre_had_matopiba, fontsize=6, color='gray')

ax6 = fig.add_subplot(3, 2, 6)
plt.plot([22,34],[22,34], color='black', linewidth=1., linestyle='--')
plt.scatter(mon_tas_cru_matopiba, np.nanmean(mon_tas_reg_exp1_matopiba, axis=1), edgecolors='blue', s=20, color='white', marker='.', label='RegCM4.7_EXP1')
plt.scatter(mon_tas_cru_matopiba, np.nanmean(mon_tas_reg_exp2_matopiba, axis=1), edgecolors='red', s=20, color='white', marker='.', label='RegCM4.7_EXP2')
plt.scatter(mon_tas_cru_matopiba, mon_tas_had_matopiba, edgecolors='gray', s=20, color='white', marker='.', label='HadGEM2-ES')
#~ plt.plot(mon_tas_cru_matopiba, m_tas_reg1_matopiba*mon_tas_cru_matopiba+c_tas_reg1_matopiba, color='blue', linewidth=1.)
#~ plt.plot(mon_tas_cru_matopiba, m_tas_reg2_matopiba*mon_tas_cru_matopiba+c_tas_reg2_matopiba, color='red', linewidth=1.)
#~ plt.plot(mon_tas_cru_matopiba, m_tas_had_matopiba*mon_tas_cru_matopiba+c_tas_had_matopiba, color='gray', linewidth=1.)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Observed Temperature (°C)', fontsize=8)
plt.xlim(22, 34)
plt.xticks(np.arange(22, 36, 2), fontsize=8)
plt.ylim(22, 34)
plt.yticks(np.arange(22, 36, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.text(23, 35, r2_tas_reg1_matopiba, fontsize=6, color='blue')
plt.text(26, 35, r2_tas_reg2_matopiba, fontsize=6, color='red')
plt.text(30, 35, r2_tas_had_matopiba, fontsize=6, color='gray')
#~ plt.text(23, 35, eq_tas_reg1_matopiba, fontsize=6, color='blue')
#~ plt.text(26, 35, eq_tas_reg2_matopiba, fontsize=6, color='red')
#~ plt.text(30, 35, eq_tas_had_matopiba, fontsize=6, color='gray')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_scatter_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()
