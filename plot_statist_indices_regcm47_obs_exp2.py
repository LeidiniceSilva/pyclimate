# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot annual climatology from regcm47 and hadgem models and obs database"

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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from comp_statist_indices import compute_corr, compute_rmse, compute_mbe


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp2'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	

	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_RegCM4_HadG_{3}_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)
	
	return value
	

def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp2/historical'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:239,:,:], axis=1), axis=1)

	return value
	
	          
# Import models and obs database 
mon_pre_cru_amz = import_obs('pre', 'amz', 'cru_ts4.04', '1986-2005')
mon_pre_gpcp_amz = import_obs('precip', 'amz', 'gpcp_v2.2', '1986-2005')
mon_pre_era5_amz = import_obs('mtpr', 'amz', 'era5', '1986-2005')
mon_pre_reg_amz = import_rcm('pr', 'amz', 'historical', '1986-2005')
mon_pre_had_amz = import_gcm('pr', 'amz', 'historical', '1986-2005')

mon_pre_cru_neb = import_obs('pre', 'neb', 'cru_ts4.04', '1986-2005')
mon_pre_gpcp_neb = import_obs('precip', 'neb', 'gpcp_v2.2', '1986-2005')
mon_pre_era5_neb = import_obs('mtpr', 'neb', 'era5', '1986-2005')
mon_pre_reg_neb = import_rcm('pr', 'neb', 'historical', '1986-2005')
mon_pre_had_neb = import_gcm('pr', 'neb', 'historical', '1986-2005')

mon_tas_cru_amz = import_obs('tmp', 'amz', 'cru_ts4.04', '1986-2005')
mon_tas_era5_amz = import_obs('t2m', 'amz', 'era5', '1986-2005')
mon_tas_reg_amz = import_rcm('tas', 'amz', 'historical', '1986-2005')
mon_tas_had_amz = import_gcm('tas', 'amz', 'historical', '1986-2005')

mon_tas_cru_neb = import_obs('tmp', 'neb', 'cru_ts4.04', '1986-2005')
mon_tas_era5_neb = import_obs('t2m', 'neb', 'era5', '1986-2005')
mon_tas_reg_neb = import_rcm('tas', 'neb', 'historical', '1986-2005')
mon_tas_had_neb = import_gcm('tas', 'neb', 'historical', '1986-2005')

# Compute correlation, root mean square error and mean bias error
corr_pre_reg_cru_amz = compute_corr(mon_pre_cru_amz, mon_pre_reg_amz)
rmse_pre_reg_cru_amz = compute_rmse(mon_pre_reg_amz, mon_pre_cru_amz)
mbe_pre_reg_cru_amz  = compute_mbe(mon_pre_reg_amz,  mon_pre_cru_amz)
corr_pre_had_cru_amz = compute_corr(mon_pre_cru_amz, mon_pre_had_amz)
rmse_pre_had_cru_amz = compute_rmse(mon_pre_had_amz, mon_pre_cru_amz)
mbe_pre_had_cru_amz  = compute_mbe(mon_pre_had_amz,  mon_pre_cru_amz)

corr_pre_reg_gpcp_amz = compute_corr(mon_pre_gpcp_amz, mon_pre_reg_amz)
rmse_pre_reg_gpcp_amz = compute_rmse(mon_pre_reg_amz,  mon_pre_gpcp_amz)
mbe_pre_reg_gpcp_amz  = compute_mbe(mon_pre_reg_amz,   mon_pre_gpcp_amz)
corr_pre_had_gpcp_amz = compute_corr(mon_pre_gpcp_amz, mon_pre_had_amz)
rmse_pre_had_gpcp_amz = compute_rmse(mon_pre_had_amz,  mon_pre_gpcp_amz)
mbe_pre_had_gpcp_amz  = compute_mbe(mon_pre_had_amz,   mon_pre_gpcp_amz)

corr_pre_reg_era5_amz = compute_corr(mon_pre_era5_amz, mon_pre_reg_amz)
rmse_pre_reg_era5_amz = compute_rmse(mon_pre_reg_amz,  mon_pre_era5_amz)
mbe_pre_reg_era5_amz  = compute_mbe(mon_pre_reg_amz,   mon_pre_era5_amz)
corr_pre_had_era5_amz = compute_corr(mon_pre_era5_amz, mon_pre_had_amz)
rmse_pre_had_era5_amz = compute_rmse(mon_pre_had_amz,  mon_pre_era5_amz)
mbe_pre_had_era5_amz  = compute_mbe(mon_pre_had_amz,   mon_pre_era5_amz)

corr_pre_reg_cru_neb = compute_corr(mon_pre_cru_neb, mon_pre_reg_neb)
rmse_pre_reg_cru_neb = compute_rmse(mon_pre_reg_neb, mon_pre_cru_neb)
mbe_pre_reg_cru_neb  = compute_mbe(mon_pre_reg_neb,  mon_pre_cru_neb)
corr_pre_had_cru_neb = compute_corr(mon_pre_cru_neb, mon_pre_had_neb)
rmse_pre_had_cru_neb = compute_rmse(mon_pre_had_neb, mon_pre_cru_neb)
mbe_pre_had_cru_neb  = compute_mbe(mon_pre_had_neb,  mon_pre_cru_neb)

corr_pre_reg_gpcp_neb = compute_corr(mon_pre_gpcp_neb, mon_pre_reg_neb)
rmse_pre_reg_gpcp_neb = compute_rmse(mon_pre_reg_neb,  mon_pre_gpcp_neb)
mbe_pre_reg_gpcp_neb  = compute_mbe(mon_pre_reg_neb,   mon_pre_gpcp_neb)
corr_pre_had_gpcp_neb = compute_corr(mon_pre_gpcp_neb, mon_pre_had_neb)
rmse_pre_had_gpcp_neb = compute_rmse(mon_pre_had_neb,  mon_pre_gpcp_neb)
mbe_pre_had_gpcp_neb  = compute_mbe(mon_pre_had_neb,   mon_pre_gpcp_neb)

corr_pre_reg_era5_neb = compute_corr(mon_pre_era5_neb, mon_pre_reg_neb)
rmse_pre_reg_era5_neb = compute_rmse(mon_pre_reg_neb,  mon_pre_era5_neb)
mbe_pre_reg_era5_neb  = compute_mbe(mon_pre_reg_neb,   mon_pre_era5_neb)
corr_pre_had_era5_neb = compute_corr(mon_pre_era5_neb, mon_pre_had_neb)
rmse_pre_had_era5_neb = compute_rmse(mon_pre_had_neb,  mon_pre_era5_neb)
mbe_pre_had_era5_neb  = compute_mbe(mon_pre_had_neb,   mon_pre_era5_neb)

corr_tas_reg_cru_amz = compute_corr(mon_tas_cru_amz, np.nanmean(mon_tas_reg_amz, axis=1))
rmse_tas_reg_cru_amz = compute_rmse(np.nanmean(mon_tas_reg_amz, axis=1), mon_tas_cru_amz)
mbe_tas_reg_cru_amz  = compute_mbe(np.nanmean(mon_tas_reg_amz, axis=1),  mon_tas_cru_amz)
corr_tas_had_cru_amz = compute_corr(mon_tas_cru_amz, mon_tas_had_amz)
rmse_tas_had_cru_amz = compute_rmse(mon_tas_had_amz, mon_tas_cru_amz)
mbe_tas_had_cru_amz  = compute_mbe(mon_tas_had_amz,  mon_tas_cru_amz)

corr_tas_reg_era5_amz = compute_corr(mon_tas_era5_amz, np.nanmean(mon_tas_reg_amz, axis=1))
rmse_tas_reg_era5_amz = compute_rmse(np.nanmean(mon_tas_reg_amz, axis=1),  mon_tas_era5_amz)
mbe_tas_reg_era5_amz  = compute_mbe(np.nanmean(mon_tas_reg_amz, axis=1),   mon_tas_era5_amz)
corr_tas_had_era5_amz = compute_corr(mon_tas_era5_amz, mon_tas_had_amz)
rmse_tas_had_era5_amz = compute_rmse(mon_tas_had_amz,  mon_tas_era5_amz)
mbe_tas_had_era5_amz  = compute_mbe(mon_tas_had_amz,   mon_tas_era5_amz)

corr_tas_reg_cru_neb = compute_corr(mon_tas_cru_neb, np.nanmean(mon_tas_reg_neb, axis=1))
rmse_tas_reg_cru_neb = compute_rmse(np.nanmean(mon_tas_reg_neb, axis=1), mon_tas_cru_neb)
mbe_tas_reg_cru_neb  = compute_mbe(np.nanmean(mon_tas_reg_neb, axis=1),  mon_tas_cru_neb)
corr_tas_had_cru_neb = compute_corr(mon_tas_cru_neb, mon_tas_had_neb)
rmse_tas_had_cru_neb = compute_rmse(mon_tas_had_neb, mon_tas_cru_neb)
mbe_tas_had_cru_neb  = compute_mbe(mon_tas_had_neb,  mon_tas_cru_neb)

corr_tas_reg_era5_neb = compute_corr(mon_tas_era5_neb, np.nanmean(mon_tas_reg_neb, axis=1))
rmse_tas_reg_era5_neb = compute_rmse(np.nanmean(mon_tas_reg_neb, axis=1),  mon_tas_era5_neb)
mbe_tas_reg_era5_neb  = compute_mbe(np.nanmean(mon_tas_reg_neb, axis=1),   mon_tas_era5_neb)
corr_tas_had_era5_neb = compute_corr(mon_tas_era5_neb, mon_tas_had_neb)
rmse_tas_had_era5_neb = compute_rmse(mon_tas_had_neb,  mon_tas_era5_neb)
mbe_tas_had_era5_neb  = compute_mbe(mon_tas_had_neb,   mon_tas_era5_neb)

# Plot models and obs database 
fig = plt.figure()

print(rmse_pre_reg_cru_amz,rmse_pre_reg_gpcp_amz,rmse_pre_reg_era5_amz,rmse_pre_had_cru_amz,rmse_pre_had_gpcp_amz,rmse_pre_had_era5_amz)
print(rmse_pre_reg_cru_neb,rmse_pre_reg_gpcp_neb,rmse_pre_reg_era5_neb,rmse_pre_had_cru_neb,rmse_pre_had_gpcp_neb,rmse_pre_had_era5_neb)
print(rmse_tas_reg_cru_amz,rmse_tas_reg_era5_amz,rmse_tas_had_cru_amz,rmse_tas_had_era5_amz)
print(rmse_tas_reg_cru_neb,rmse_tas_reg_era5_neb,rmse_tas_had_cru_neb,rmse_tas_had_era5_neb)

ax = fig.add_subplot(2, 2, 1)
plt.plot([0, 2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot([0, -2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot(mbe_pre_reg_cru_amz,  corr_pre_reg_cru_amz,  markersize=rmse_pre_reg_cru_amz+5,  marker='<', markerfacecolor='gray',   color='black', label='RegCM4 (CRU)')
plt.plot(mbe_pre_reg_gpcp_amz, corr_pre_reg_gpcp_amz, markersize=rmse_pre_reg_gpcp_amz+5, marker='s', markerfacecolor='blue',   color='black', label='RegCM4 (GPCP)')
plt.plot(mbe_pre_reg_era5_amz, corr_pre_reg_era5_amz, markersize=rmse_pre_reg_era5_amz+5, marker='^', markerfacecolor='green',  color='black', label='RegCM4 (ERA5)')
plt.plot(mbe_pre_had_cru_amz,  corr_pre_had_cru_amz,  markersize=rmse_pre_had_cru_amz+5,  marker='h', markerfacecolor='red',    color='black', label='HadGEM2-ES (CRU)')
plt.plot(mbe_pre_had_gpcp_amz, corr_pre_had_gpcp_amz, markersize=rmse_pre_had_gpcp_amz+5, marker='>', markerfacecolor='purple', color='black', label='HadGEM2-ES (GPCP)')
plt.plot(mbe_pre_had_era5_amz, corr_pre_had_era5_amz, markersize=rmse_pre_had_era5_amz+5, marker='o', markerfacecolor='orange', color='black', label='HadGEM2-ES (ERA5)')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'PCC', fontsize=8, fontweight='bold')
plt.ylim(-1, 1)
plt.xlim(-2, 2)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.xticks(np.arange(-2, 2.5, 0.5), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(0, linewidth=1., linestyle='-', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(handlelength=0.25, handleheight=0.25, fontsize=6, loc=4, ncol=1)
plt.text(-1.6,  0.85, 'RMSE = 1.26', color='orange', fontsize=8, fontweight='bold')
plt.text(0.5,  0.7, 'RMSE = 2.63', color='red', fontsize=8, fontweight='bold')
       
ax = fig.add_subplot(2, 2, 2)
plt.plot([0, 2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot([0, -2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot(mbe_tas_reg_cru_amz,  corr_tas_reg_cru_amz,  markersize=rmse_tas_reg_cru_amz+5,  marker='<', markerfacecolor='gray',   color='black', label='RegCM4 (CRU)')
plt.plot(mbe_tas_reg_era5_amz, corr_tas_reg_era5_amz, markersize=rmse_tas_reg_era5_amz+5, marker='^', markerfacecolor='green',  color='black', label='RegCM4 (ERA5)')
plt.plot(mbe_tas_had_cru_amz,  corr_tas_had_cru_amz,  markersize=rmse_tas_had_cru_amz+5,  marker='h', markerfacecolor='red',    color='black', label='HadGEM2-ES (CRU)')
plt.plot(mbe_tas_had_era5_amz, corr_tas_had_era5_amz, markersize=rmse_tas_had_era5_amz+5, marker='o', markerfacecolor='orange', color='black', label='HadGEM2-ES (ERA5)')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.ylim(-1, 1)
plt.xlim(-2, 2)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.xticks(np.arange(-2, 2.5, 0.5), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(0, linewidth=1., linestyle='-', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.text(-0.1,  0.55, 'RMSE = 1.05', color='orange', fontsize=8, fontweight='bold')
plt.text(-1.8,  0.35, 'RMSE = 1.77', color='gray', fontsize=8, fontweight='bold')
   
ax = fig.add_subplot(2, 2, 3)
plt.plot([0, 2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot([0, -2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot(mbe_pre_reg_cru_neb,  corr_pre_reg_cru_neb,  markersize=rmse_pre_reg_cru_neb+5,  marker='<', markerfacecolor='gray',   color='black', label='RegCM4 (CRU)')
plt.plot(mbe_pre_reg_gpcp_neb, corr_pre_reg_gpcp_neb, markersize=rmse_pre_reg_gpcp_neb+5, marker='s', markerfacecolor='blue',   color='black', label='RegCM4 (GPCP)')
plt.plot(mbe_pre_reg_era5_neb, corr_pre_reg_era5_neb, markersize=rmse_pre_reg_era5_neb+5, marker='^', markerfacecolor='green',  color='black', label='RegCM4 (ERA5)')
plt.plot(mbe_pre_had_cru_neb,  corr_pre_had_cru_neb,  markersize=rmse_pre_had_cru_neb+5,  marker='h', markerfacecolor='red',    color='black', label='HadGEM2-ES (CRU)')
plt.plot(mbe_pre_had_gpcp_neb, corr_pre_had_gpcp_neb, markersize=rmse_pre_had_gpcp_neb+5, marker='>', markerfacecolor='purple', color='black', label='HadGEM2-ES (GPCP)')
plt.plot(mbe_pre_had_era5_neb, corr_pre_had_era5_neb, markersize=rmse_pre_had_era5_neb+5, marker='o', markerfacecolor='orange', color='black', label='HadGEM2-ES (ERA5)')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'MBE (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.ylabel(u'PCC', fontsize=8, fontweight='bold')
plt.ylim(-1, 1)
plt.xlim(-2, 2)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.xticks(np.arange(-2, 2.5, 0.5), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(0, linewidth=1., linestyle='-', color='black')
plt.text(0.45,  0.85, 'RMSE = 2.06', color='blue', fontsize=8, fontweight='bold')
plt.text(0.6,  0.5, 'RMSE = 2.12', color='gray', fontsize=8, fontweight='bold')
   
ax = fig.add_subplot(2, 2, 4)
plt.plot([0, 2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot([0, -2], [0, 1], linewidth=1., linestyle='--', color='black')
plt.plot(mbe_tas_reg_cru_neb,  corr_tas_reg_cru_neb,  markersize=rmse_tas_reg_cru_neb+5,  marker='<', markerfacecolor='gray',   color='black', label='RegCM4 (CRU)')
plt.plot(mbe_tas_reg_era5_neb, corr_tas_reg_era5_neb, markersize=rmse_tas_reg_era5_neb+5, marker='^', markerfacecolor='green',  color='black', label='RegCM4 (ERA5)')
plt.plot(mbe_tas_had_cru_neb,  corr_tas_had_cru_neb,  markersize=rmse_tas_had_cru_neb+5,  marker='h', markerfacecolor='red',    color='black', label='HadGEM2-ES (CRU)')
plt.plot(mbe_tas_had_era5_neb, corr_tas_had_era5_neb, markersize=rmse_tas_had_era5_neb+5, marker='o', markerfacecolor='orange', color='black', label='HadGEM2-ES (ERA5)')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'MBE (°C)', fontsize=8, fontweight='bold')
plt.ylim(-1, 1)
plt.xlim(-2, 2)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.xticks(np.arange(-2, 2.5, 0.5), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(0, linewidth=1., linestyle='-', color='black')
plt.text(-0.5,  0.8, 'RMSE = 1.02', color='gray', fontsize=8, fontweight='bold')
plt.text(-1.9,  0.6, 'RMSE = 1.31', color='orange', fontsize=8, fontweight='bold')
   
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_statist_indices_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
