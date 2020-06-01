# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "03/25/2019"
__description__ = "This script plot correlation Portrait Diagram from CMIP5 models"

import os
import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from comp_statist_indices import compute_corr
from comp_statist_indices import compute_rmse
from comp_statist_indices import compute_mae
from comp_statist_indices import compute_r2

def import_cmip5(model):
	
	param = 'tas' # pr or tas
	area  = 'amz' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	month_sim = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return month_sim


def import_obs(database):
	
	param = 'tmp' # pre or tmp
	area  = 'amz' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	month_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return month_obs
	

models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

metrics_gcm = []

for mdl in models:
	
	# Import cmip5 model end obs database monthly	
	month_gcm = import_cmip5(mdl)
		
	obs  = u'cru_ts4.02'
	month_cru = import_obs(obs)
	
	corr = compute_corr(month_gcm, month_cru)
	r2 = compute_r2(month_gcm, month_cru)
	rmse = compute_rmse(month_gcm, month_cru)
	mae = compute_mae(month_gcm, month_cru)
	
	metrics = [corr, r2, rmse, mae]
	metrics_gcm.append(metrics)
	
# Compute TOPSIS to CMIP5 models
# Method 1 (Import topsis)
a = metrics_gcm
w = [0.25, 0.25, 0.25, 0.25]
j = [1, 1, 1, 1]

print(metrics_gcm)

exit()

# Compute TOPSIS to CMIP5 models
# Method 2 (Steps)
x = np.array([[8,7,2,1], [5,3,7,5], [7,5,6,4], 
		[9,9,7,3], [11,10,3,7], [6,9,5,4]])
weigths = np.array([0.4, 0.3, 0.1, 0.2])

# Step 1 - Vertor normalization: cumsum() produces the cumulative sum of the values in the array 
# and the can also be used with a second argument to indicated the axix to use
col_sums = np.array(sum(x**2, 0))

norm_x = np.array([[round(x[i, j] / np.sqrt(col_sums[x.shape[0] 
	- 1, j]), 3) for j in range (4)] for i in range(6)])

# Step  2 - Multiple each evaluation by the associated weigth:
# wnx is the weighted normalizerd x matrix
wnx = np.array([[round(weigth[i] * norm_x[j, i], 3) 
	for i in range(4)] for j in range(6)])

# Step3 - Positive and negative idea solution
pis = np.array([amax(wnx[:, :1]), amax(wnx[:, 1:2]), 
	amax(wnx[:, 2:3]), amax(wnx[:, 3:4])])
nis = np.array([amin(wnx[:, :1]), amin(wnx[:, 1:2]), 
	amin(wnx[:, 2:3]), amin(wnx[:, 3:4])])

# Step 4a - Determine the distance  to the positive ideal solution (dpis)
b1 = np.array([[(wnx[j, i] - pis[i])**2 for i in range(4)] 
	for j in range(6)])
dpis = np.sqrt(sum(b1, 1))

# Step 4b - Determine the distance  to the negative ideal solution (dnis)
b2 = np.array([[(wnx[j, i] - nis[i])**2 for i in range(4)] 
	for j in range(6)])
dnis = np.sqrt(sum(b2, 1))

# Step 5 - Calculate the relative closeness to the ideal solution
final_solution = np.array([round(dnis[i] / (dpis[i] + dnis[i]), 
	3) for i in range(6)])
print('Closeness coeficient = ', final_solution)
exit()

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_topsis_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


