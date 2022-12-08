# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot scatter from regcm46 and obs database"

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
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties


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
	mean_obs = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_obs
	
	
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
	mean_sim = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_sim
	
	
def compute_linear_regression(obs, sim):

	m, c, r, p, se1 = stats.linregress(obs, sim)
	eq = 'y = %2.2fx+%2.2f'%(m, c)
	r2 = 'R²=%1.2f'%(r**2)
	p = 'p-value=%1.3f'%(p)
   
	return p, r2

       
# Import regcm exps and obs database 
namz_obs_djf = import_obs(u'namz', u'gpcp_v2.3_obs', 'djf')
namz_exp1_djf = import_sim(u'namz', u'regcm_exp1', 'djf')
namz_exp2_djf = import_sim(u'namz', u'regcm_exp2', 'djf')

samz_obs_djf = import_obs(u'samz', u'gpcp_v2.3_obs', 'djf')
samz_exp1_djf = import_sim(u'samz', u'regcm_exp1', 'djf')
samz_exp2_djf = import_sim(u'samz', u'regcm_exp2', 'djf')

neb_obs_djf = import_obs(u'neb', u'gpcp_v2.3_obs', 'djf')
neb_exp1_djf = import_sim(u'neb', u'regcm_exp1', 'djf')
neb_exp2_djf = import_sim(u'neb', u'regcm_exp2', 'djf')

namz_obs_djf = namz_obs_djf[2]
namz_exp1_djf = namz_exp1_djf[2]
namz_exp2_djf = namz_exp2_djf[2]

samz_obs_djf = samz_obs_djf[2]
samz_exp1_djf = samz_exp1_djf[2]
samz_exp2_djf = samz_exp2_djf[2]

neb_obs_djf = neb_obs_djf[2]
neb_exp1_djf = neb_exp1_djf[2]
neb_exp2_djf = neb_exp2_djf[2]

# Import linear regression
p_namz_exp1_djf, r2_namz_exp1_djf = compute_linear_regression(namz_obs_djf, namz_exp1_djf)
p_samz_exp1_djf, r2_samz_exp1_djf = compute_linear_regression(samz_obs_djf, samz_exp1_djf)
p_neb_exp1_djf, r2_neb_exp1_djf = compute_linear_regression(neb_obs_djf, neb_exp1_djf)

p_namz_exp2_djf, r2_namz_exp2_djf = compute_linear_regression(namz_obs_djf, namz_exp2_djf)
p_samz_exp2_djf, r2_samz_exp2_djf = compute_linear_regression(samz_obs_djf, samz_exp2_djf)
p_neb_exp2_djf, r2_neb_exp2_djf = compute_linear_regression(neb_obs_djf, neb_exp2_djf)

# Plot models and obs database 
fig = plt.figure() 
gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1]) 
time = np.arange(0, 9 + 1)

ax1 = plt.subplot(gs[0])
plt.plot(time, namz_obs_djf, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, namz_exp1_djf, color='blue', marker='.', markerfacecolor='white', label='Reg_Holtslag')
plt.plot(time, namz_exp2_djf, color='red', marker='.', markerfacecolor='white', label='Reg_UW-PBL')
plt.title(u'A) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.legend(fontsize=7, loc=9, ncol=3, shadow=True)

ax2 = plt.subplot(gs[1])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(namz_obs_djf, namz_exp1_djf, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Holtslag')
plt.scatter(namz_obs_djf, namz_exp2_djf, edgecolors='red', s=40, color='white', marker='.', label='Reg_UW-PBL')
plt.title(u'B) NAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xlim(0, 14)
plt.ylim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.text(1, 12, r2_namz_exp1_djf, fontsize=8, color='blue')
plt.text(6, 12, p_namz_exp1_djf, fontsize=8, color='blue')
plt.text(1, 10, r2_namz_exp2_djf, fontsize=8, color='red')
plt.text(6, 10, p_namz_exp2_djf, fontsize=8, color='red')

ax3 = plt.subplot(gs[2])
plt.plot(time, samz_obs_djf, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, samz_exp1_djf, color='blue', marker='.', markerfacecolor='white', label='Reg_Holtslag')
plt.plot(time, samz_exp2_djf, color='red', marker='.', markerfacecolor='white', label='Reg_UW-PBL')
plt.title(u'C) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax3.get_xticklabels(), visible=False)

ax4 = plt.subplot(gs[3])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(samz_obs_djf, samz_exp1_djf, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Holtslag')
plt.scatter(samz_obs_djf, samz_exp2_djf, edgecolors='red', s=40, color='white', marker='.', label='Reg_UW-PBL')
plt.title(u'D) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Simulated precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 14)
plt.ylim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.text(1, 12, r2_samz_exp1_djf, fontsize=8, color='blue')
plt.text(6, 12, p_samz_exp1_djf, fontsize=8, color='blue')
plt.text(1, 10, r2_samz_exp2_djf, fontsize=8, color='red')
plt.text(6, 10, p_samz_exp2_djf, fontsize=8, color='red')

ax5 = plt.subplot(gs[4])
plt.plot(time, neb_obs_djf, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, neb_exp1_djf, color='blue', marker='.', markerfacecolor='white', label='Reg_Holtslag')
plt.plot(time, neb_exp2_djf, color='red', marker='.', markerfacecolor='white', label='Reg_UW-PB')
plt.title(u'E) NEB', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Years', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')

ax6 = plt.subplot(gs[5])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(neb_obs_djf, neb_exp1_djf, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Holtslag')
plt.scatter(neb_obs_djf, neb_exp2_djf, edgecolors='red', s=40, color='white', marker='.', label='Reg_UW-PBL')
plt.title(u'F) NEB', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Observed precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 14)
plt.ylim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.text(1, 12, r2_neb_exp1_djf, fontsize=8, color='blue')
plt.text(6, 12, p_neb_exp1_djf, fontsize=8, color='blue')
plt.text(1, 10, r2_neb_exp2_djf, fontsize=8, color='red')
plt.text(6, 10, p_neb_exp2_djf, fontsize=8, color='red')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_scatter_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
