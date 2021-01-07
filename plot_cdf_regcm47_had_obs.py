# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute statiscs indices from Reg and Had models"

import os
import netCDF4
import numpy as np
import numpy as np
import scipy.stats as s
import matplotlib.pyplot as plt

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return gcm

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return obs
	
	
# Import regcm exps model end obs database climatology
p_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
p_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
p_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
p_udel_samz = import_obs('pre', 'samz', 'udel_v301', '1986-2005')
p_chirps_samz = import_obs('precip', 'samz', 'chirps-v2.0', '1986-2005')
p_era5_samz = import_obs('mtpr', 'samz', 'era5', '1986-2005')

p_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
p_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
p_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
p_udel_eneb = import_obs('pre', 'eneb', 'udel_v301', '1986-2005')
p_chirps_eneb = import_obs('precip', 'eneb', 'chirps-v2.0', '1986-2005')
p_era5_eneb = import_obs('mtpr', 'eneb', 'era5', '1986-2005')

p_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
p_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
p_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
p_udel_matopiba = import_obs('pre', 'matopiba', 'udel_v301', '1986-2005')
p_chirps_matopiba = import_obs('precip', 'matopiba', 'chirps-v2.0', '1986-2005')
p_era5_matopiba = import_obs('mtpr', 'matopiba', 'era5', '1986-2005')

t_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
t_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
t_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
t_udel_samz = import_obs('temp', 'samz', 'udel_v301', '1986-2005')
t_era5_samz = import_obs('t2m', 'samz', 'era5', '1986-2005')

t_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
t_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
t_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
t_udel_eneb = import_obs('temp', 'eneb', 'udel_v301', '1986-2005')
t_era5_eneb = import_obs('t2m', 'eneb', 'era5', '1986-2005')

t_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
t_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
t_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
t_udel_matopiba = import_obs('temp', 'matopiba', 'udel_v301', '1986-2005')
t_era5_matopiba = import_obs('t2m', 'matopiba', 'era5', '1986-2005')

wei = s.weibull_min(2, 0, 2) # shape, loc, scale - creates weibull object
sample = wei.rvs(1000)
shape, loc, scale = s.weibull_min.fit(sample, floc=0) 

x = np.linspace(np.min(sample), np.max(sample))
dx = x[1]-x[0]
deriv = np.diff(wei.cdf(x))/dx
plt.plot(x, wei.cdf(x), label="cdf")
plt.plot(x, wei.pdf(x), label="pdf")
plt.legend(loc=1)
plt.show()
exit()

