# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/10/2021"
__description__ = "This script plot scatter from regcm46 and obs database"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib import gridspec

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties


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
	
	obs_ann = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1)
	
	return obs_ann
	

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

	sim_ann = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1)
	
	return sim_ann
	
	
def compute_linear_regression(obs, sim):

	m, c, r, p, se1 = stats.linregress(obs, sim)
	eq = 'y = %2.2fx+%2.2f'%(m, c)
	r2 = 'R² = %1.2f'%(r)
	p = 'p-value = %1.2f'%(p)
   
	return m, c, p, eq, r2

       
# Import regcm exps and obs database 
namz_exp1_ann = import_sim(u'namz', u'regcm_exp1')
samz_exp1_ann = import_sim(u'samz', u'regcm_exp1')
neb_exp1_ann = import_sim(u'neb', u'regcm_exp1')

namz_exp2_ann = import_sim(u'namz', u'regcm_exp2')
samz_exp2_ann = import_sim(u'samz', u'regcm_exp2')
neb_exp2_ann = import_sim(u'neb', u'regcm_exp2')

namz_obs_ann = import_obs(u'namz', u'gpcp_v2.3_obs')
samz_obs_ann = import_obs(u'samz', u'gpcp_v2.3_obs')
neb_obs_ann = import_obs(u'neb', u'gpcp_v2.3_obs')

# Import linear regression
m_namz_exp1, c_namz_exp1, p_namz_exp1, eq_namz_exp1, r2_namz_exp1 = compute_linear_regression(namz_obs_ann, namz_exp1_ann)
m_samz_exp1, c_samz_exp1, p_samz_exp1, eq_samz_exp1, r2_samz_exp1 = compute_linear_regression(samz_obs_ann, samz_exp1_ann)
m_neb_exp1, c_neb_exp1, p_neb_exp1, eq_neb_exp1, r2_neb_exp1 = compute_linear_regression(neb_obs_ann, neb_exp1_ann)

m_namz_exp2, c_namz_exp2,  p_namz_exp2, eq_namz_exp2, r2_namz_exp2 = compute_linear_regression(namz_obs_ann, namz_exp2_ann)
m_samz_exp2, c_samz_exp2,  p_samz_exp2, eq_samz_exp2, r2_samz_exp2 = compute_linear_regression(samz_obs_ann, samz_exp2_ann)
m_neb_exp2, c_neb_exp2,  p_neb_exp2, eq_neb_exp2, r2_neb_exp2 = compute_linear_regression(neb_obs_ann, neb_exp2_ann)

# Plot models and obs database 
fig = plt.figure() 
gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1]) 
time = np.arange(1, 10 + 1)

ax1 = plt.subplot(gs[0])
plt.plot(time, namz_obs_ann, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, namz_exp1_ann, color='blue', marker='.', markerfacecolor='white', label='Reg_Exp1')
plt.plot(time, namz_exp2_ann, color='red', marker='.', markerfacecolor='white', label='Reg_Exp2')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.text(1, 12, r2_namz_exp1, fontsize=8, color='blue')
plt.text(3, 12, p_namz_exp1, fontsize=8, color='blue')
plt.text(1, 10, r2_namz_exp2, fontsize=8, color='red')
plt.text(3, 10, p_namz_exp2, fontsize=8, color='red')

ax2 = plt.subplot(gs[1])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(namz_obs_ann, namz_exp1_ann, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Exp1')
plt.scatter(namz_obs_ann, namz_exp2_ann, edgecolors='red', s=40, color='white', marker='.', label='Reg_Exp2')
plt.plot(namz_obs_ann, m_namz_exp1*namz_obs_ann+c_namz_exp1, color='blue', linewidth=1.)
plt.plot(namz_obs_ann, m_namz_exp2*namz_obs_ann+c_namz_exp2, color='red', linewidth=1.)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.text(0.5, 12, eq_namz_exp1, fontsize=8, color='blue')
plt.text(0.5, 10, eq_namz_exp2, fontsize=8, color='red')

ax3 = plt.subplot(gs[2])
plt.plot(time, samz_obs_ann, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, samz_exp1_ann, color='blue', marker='.', markerfacecolor='white', label='Reg_Exp1')
plt.plot(time, samz_exp2_ann, color='red', marker='.', markerfacecolor='white', label='Reg_Exp2')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.text(1, 2, r2_samz_exp1, fontsize=8, color='blue')
plt.text(3, 2, p_samz_exp1, fontsize=8, color='blue')
plt.text(1, 0, r2_samz_exp2, fontsize=8, color='red')
plt.text(3, 0, p_samz_exp2, fontsize=8, color='red')

ax4 = plt.subplot(gs[3])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(samz_obs_ann, samz_exp1_ann, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Exp1')
plt.scatter(samz_obs_ann, samz_exp2_ann, edgecolors='red', s=40, color='white', marker='.', label='Reg_Exp2')
plt.plot(samz_obs_ann, m_samz_exp1*samz_obs_ann+c_samz_exp1, color='blue', linewidth=1.)
plt.plot(samz_obs_ann, m_samz_exp2*samz_obs_ann+c_samz_exp2, color='red', linewidth=1.)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Simulated Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.text(0.5, 12, eq_samz_exp1, fontsize=8, color='blue')
plt.text(0.5, 10, eq_samz_exp2, fontsize=8, color='red')

ax5 = plt.subplot(gs[4])
plt.plot(time, neb_obs_ann, color='black', marker='.', markerfacecolor='white', label='GPCP')
plt.plot(time, neb_exp1_ann, color='blue', marker='.', markerfacecolor='white', label='Reg_Exp1')
plt.plot(time, neb_exp2_ann, color='red', marker='.', markerfacecolor='white', label='Reg_Exp2')
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Years', fontsize=8, fontweight='bold')
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'), fontsize=8)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.legend(fontsize=8, loc=1, ncol=1, shadow=True)
plt.text(1, 12, r2_neb_exp1, fontsize=8, color='blue')
plt.text(3, 12, p_neb_exp1, fontsize=8, color='blue')
plt.text(1, 10, r2_neb_exp2, fontsize=8, color='red')
plt.text(3, 10, p_neb_exp2, fontsize=8, color='red')

ax6 = plt.subplot(gs[5])
plt.plot([0,14],[0,14], color='black', linewidth=1., linestyle='--')
plt.scatter(neb_obs_ann, neb_exp1_ann, edgecolors='blue', s=40, color='white', marker='.', label='Reg_Exp1')
plt.scatter(neb_obs_ann, neb_exp2_ann, edgecolors='red', s=40, color='white', marker='.', label='Reg_Exp2')
plt.plot(neb_obs_ann, m_neb_exp1*neb_obs_ann+c_neb_exp1, color='blue', linewidth=1.)
plt.plot(neb_obs_ann, m_neb_exp2*neb_obs_ann+c_neb_exp2, color='red', linewidth=1.)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Observed Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 14)
plt.xticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.grid(True, which='major', linestyle='--')
plt.text(0.5, 12, eq_neb_exp1, fontsize=8, color='blue')
plt.text(0.5, 10, eq_neb_exp2, fontsize=8, color='red')

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/figs'
name_out = 'pyplt_scatter_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
