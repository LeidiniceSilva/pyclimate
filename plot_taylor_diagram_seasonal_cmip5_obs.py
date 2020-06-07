# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/15/2019"
__description__ = "This script plot taylor diagram graphics from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from plot_taylor_diagram_monthly_cmip5_obs import TaylorDiagram


def import_cmip5_season(model):
		
	param = 'pr' # pr or tas
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5/hist'
	arq   = '{0}/{1}_amz_neb_Amon_{2}_{3}_{4}.nc'.format(path, param,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	sea_mdl = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	sea1_mdl = sea_mdl[0:120:4]
	sea2_mdl = sea_mdl[1:120:4]
	sea3_mdl = sea_mdl[2:120:4]
	sea4_mdl = sea_mdl[3:120:4]

	return sea1_mdl, sea2_mdl, sea3_mdl, sea4_mdl


def import_obs_season(database):
	
	param = 'pre' # pre or tmp
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_amz_neb_{2}_obs_mon_{3}.nc'.format(path,
	param, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	sea_obs = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	sea1_obs = sea_obs[0:120:4]
	sea2_obs = sea_obs[1:120:4]
	sea3_obs = sea_obs[2:120:4]
	sea4_obs = sea_obs[3:120:4]

	return sea1_obs, sea2_obs, sea3_obs, sea4_obs


# Import cmip5 model end obs database seasonaly
model  = u'BCC-CSM1.1'
sea1_mdl1, sea2_mdl1, sea3_mdl1, sea4_mdl1 = import_cmip5_season(model)

model  = u'BCC-CSM1.1M'
sea1_mdl2, sea2_mdl2, sea3_mdl2, sea4_mdl2 = import_cmip5_season(model)

model  = u'BNU-ESM'
sea1_mdl3, sea2_mdl3, sea3_mdl3, sea4_mdl3 = import_cmip5_season(model)

model  = u'CanESM2'
sea1_mdl4, sea2_mdl4, sea3_mdl4, sea4_mdl4 = import_cmip5_season(model)

model  = u'CNRM-CM5'
sea1_mdl5, sea2_mdl5, sea3_mdl5, sea4_mdl5 = import_cmip5_season(model)

model  = u'CSIRO-ACCESS-1'
sea1_mdl6, sea2_mdl6, sea3_mdl6, sea4_mdl6 = import_cmip5_season(model)

model  = u'CSIRO-ACCESS-3'
sea1_mdl7, sea2_mdl7, sea3_mdl7, sea4_mdl7 = import_cmip5_season(model)

model  = u'CSIRO-MK36'
sea1_mdl8, sea2_mdl8, sea3_mdl8, sea4_mdl8 = import_cmip5_season(model)

model  = u'FIO-ESM'
sea1_mdl9, sea2_mdl9, sea3_mdl9, sea4_mdl9 = import_cmip5_season(model)

model  = u'GISS-E2-H-CC'
sea1_mdl10, sea2_mdl10, sea3_mdl10, sea4_mdl10 = import_cmip5_season(model)

model  = u'GISS-E2-H'
sea1_mdl11, sea2_mdl11, sea3_mdl11, sea4_mdl11 = import_cmip5_season(model)

model  = u'GISS-E2-R'
sea1_mdl12, sea2_mdl12, sea3_mdl12, sea4_mdl12 = import_cmip5_season(model)

model  = u'HadGEM2-AO'
sea1_mdl13, sea2_mdl13, sea3_mdl13, sea4_mdl13 = import_cmip5_season(model)

model  = u'HadGEM2-CC'
sea1_mdl14, sea2_mdl14, sea3_mdl14, sea4_mdl14 = import_cmip5_season(model)

model  = u'HadGEM2-ES'
sea1_mdl15, sea2_mdl15, sea3_mdl15, sea4_mdl15 = import_cmip5_season(model)

model  = u'INMCM4'
sea1_mdl16, sea2_mdl16, sea3_mdl16, sea4_mdl16 = import_cmip5_season(model)

model  = u'IPSL-CM5A-LR'
sea1_mdl17, sea2_mdl17, sea3_mdl17, sea4_mdl17 = import_cmip5_season(model)

model  = u'IPSL-CM5A-MR'
sea1_mdl18, sea2_mdl18, sea3_mdl18, sea4_mdl18 = import_cmip5_season(model)

model  = u'IPSL-CM5B-LR'
sea1_mdl19, sea2_mdl19, sea3_mdl19, sea4_mdl19 = import_cmip5_season(model)

model  = u'LASG-FGOALS-G2'
sea1_mdl20, sea2_mdl20, sea3_mdl20, sea4_mdl20 = import_cmip5_season(model)

model  = u'LASG-FGOALS-S2'
sea1_mdl21, sea2_mdl21, sea3_mdl21, sea4_mdl21 = import_cmip5_season(model)

model  = u'MIROC5'
sea1_mdl22, sea2_mdl22, sea3_mdl22, sea4_mdl22 = import_cmip5_season(model)

model  = u'MIROC-ESM-CHEM'
sea1_mdl23, sea2_mdl23, sea3_mdl23, sea4_mdl23 = import_cmip5_season(model)

model  = u'MIROC-ESM'
sea1_mdl24, sea2_mdl24, sea3_mdl24, sea4_mdl24 = import_cmip5_season(model)

model  = u'MPI-ESM-LR'
sea1_mdl25, sea2_mdl25, sea3_mdl25, sea4_mdl25 = import_cmip5_season(model)

model  = u'MPI-ESM-MR'
sea1_mdl26, sea2_mdl26, sea3_mdl26, sea4_mdl26 = import_cmip5_season(model)

model  = u'MRI-CGCM3'
sea1_mdl27, sea2_mdl27, sea3_mdl27, sea4_mdl27 = import_cmip5_season(model)

model  = u'NCAR-CCSM4'
sea1_mdl28, sea2_mdl28, sea3_mdl28, sea4_mdl28 = import_cmip5_season(model)

model  = u'NCAR-CESM1-BGC'
sea1_mdl29, sea2_mdl29, sea3_mdl29, sea4_mdl29 = import_cmip5_season(model)

model  = u'NCAR-CESM1-CAM5'
sea1_mdl30, sea2_mdl30, sea3_mdl30, sea4_mdl30 = import_cmip5_season(model)

model  = u'NorESM1-ME'
sea1_mdl31, sea2_mdl31, sea3_mdl31, sea4_mdl31 = import_cmip5_season(model)

model  = u'NorESM1-M'
sea1_mdl32, sea2_mdl32, sea3_mdl32, sea4_mdl32 = import_cmip5_season(model)

model  = u'ensmean_cmip5'
sea1_mdl33, sea2_mdl33, sea3_mdl33, sea4_mdl33 = import_cmip5_season(model)

model  = u'cru_ts4.02'
sea1_obs1, sea2_obs1, sea3_obs1, sea4_obs1 = import_obs_season(model)

# Reference database standard desviation
DJF = np.cos(sea1_obs1)
MAM = np.cos(sea2_obs1)
JJA = np.cos(sea3_obs1)
SON = np.cos(sea4_obs1)

refstd1 = DJF.std(ddof=1)
refstd2 = MAM.std(ddof=1)
refstd3 = JJA.std(ddof=1)
refstd4 = SON.std(ddof=1)

stdrefs = dict(DJF=refstd1,
               MAM=refstd2,
               JJA=refstd3,
               SON=refstd4)           

# Sample std, rho: Be sure to check order and that correct numbers are placed!
samples = dict(DJF=[[sea1_mdl1.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl1)[0,1], 'BCC-CSM1.1'],
                       [sea1_mdl2.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl2)[0,1], 'BCC-CSM1.1M'],
                       [sea1_mdl3.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl3)[0,1], 'BNU-ESM'],
                       [sea1_mdl4.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl4)[0,1], 'CanESM2'],
                       [sea1_mdl5.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl5)[0,1], 'CNRM-CM5'],
                       [sea1_mdl6.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl6)[0,1], 'CSIRO-ACCESS-1'],
                       [sea1_mdl7.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl7)[0,1], 'CSIRO-ACCESS-3'],
                       [sea1_mdl8.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl8)[0,1], 'CSIRO-MK36'],
                       [sea1_mdl9.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl9)[0,1], 'FIO-ESM'],
                       [sea1_mdl10.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl10)[0,1], 'GISS-E2-H-CC'],
                       [sea1_mdl11.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl11)[0,1], 'GISS-E2-H'],
                       [sea1_mdl12.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl12)[0,1], 'GISS-E2-R'],
                       [sea1_mdl13.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl13)[0,1], 'HadGEM2-AO'],
                       [sea1_mdl14.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl14)[0,1], 'HadGEM2-CC'],
                       [sea1_mdl15.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl15)[0,1], 'HadGEM2-ES'],
                       [sea1_mdl16.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl16)[0,1], 'INMCM4'],
                       [sea1_mdl17.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl17)[0,1], 'IPSL-CM5A-LR'],
                       [sea1_mdl18.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl18)[0,1], 'IPSL-CM5A-MR'],
                       [sea1_mdl19.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl19)[0,1], 'IPSL-CM5B-LR'],
                       [sea1_mdl20.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl20)[0,1], 'LASG-FGOALS-G2'],
                       [sea1_mdl21.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl21)[0,1], 'LASG-FGOALS-S2'],
                       [sea1_mdl22.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl22)[0,1], 'MIROC5'],
                       [sea1_mdl23.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl23)[0,1], 'MIROC-ESM-CHEM'],
                       [sea1_mdl24.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl24)[0,1], 'MIROC-ESM'],
                       [sea1_mdl25.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl25)[0,1], 'MPI-ESM-LR'],
                       [sea1_mdl26.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl26)[0,1], 'MPI-ESM-MR'],
                       [sea1_mdl27.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl27)[0,1], 'MRI-CGCM3'],
                       [sea1_mdl28.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl28)[0,1], 'NCAR-CCSM4'],
                       [sea1_mdl29.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl29)[0,1], 'NCAR-CESM1-BGC'],
                       [sea1_mdl30.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl30)[0,1], 'NCAR-CESM1-CAM5'],
                       [sea1_mdl31.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl31)[0,1], 'NorESM1-ME'],
                       [sea1_mdl32.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl32)[0,1], 'NorESM1-M'],
                       [sea1_mdl33.std(ddof=1), np.corrcoef(sea1_obs1, sea1_mdl33)[0,1], 'ensmean_cmip5']],              
               MAM=[[sea2_mdl1.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl1)[0,1], 'BCC-CSM1.1'],
                       [sea2_mdl2.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl2)[0,1], 'BCC-CSM1.1M'],
                       [sea2_mdl3.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl3)[0,1], 'BNU-ESM'],
                       [sea2_mdl4.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl4)[0,1], 'CanESM2'],
                       [sea2_mdl5.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl5)[0,1], 'CNRM-CM5'],
                       [sea2_mdl6.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl6)[0,1], 'CSIRO-ACCESS-1'],
                       [sea2_mdl7.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl7)[0,1], 'CSIRO-ACCESS-3'],
                       [sea2_mdl8.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl8)[0,1], 'CSIRO-MK36'],
                       [sea2_mdl9.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl9)[0,1], 'FIO-ESM'],
                       [sea2_mdl10.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl10)[0,1], 'GISS-E2-H-CC'],
                       [sea2_mdl11.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl11)[0,1], 'GISS-E2-H'],
                       [sea2_mdl12.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl12)[0,1], 'GISS-E2-R'],
                       [sea2_mdl13.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl13)[0,1], 'HadGEM2-AO'],
                       [sea2_mdl14.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl14)[0,1], 'HadGEM2-CC'],
                       [sea2_mdl15.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl15)[0,1], 'HadGEM2-ES'],
                       [sea2_mdl16.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl16)[0,1], 'INMCM4'],
                       [sea2_mdl17.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl17)[0,1], 'IPSL-CM5A-LR'],
                       [sea2_mdl18.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl18)[0,1], 'IPSL-CM5A-MR'],
                       [sea2_mdl19.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl19)[0,1], 'IPSL-CM5B-LR'],
                       [sea2_mdl20.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl20)[0,1], 'LASG-FGOALS-G2'],
                       [sea2_mdl21.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl21)[0,1], 'LASG-FGOALS-S2'],
                       [sea2_mdl22.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl22)[0,1], 'MIROC5'],
                       [sea2_mdl23.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl23)[0,1], 'MIROC-ESM-CHEM'],
                       [sea2_mdl24.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl24)[0,1], 'MIROC-ESM'],
                       [sea2_mdl25.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl25)[0,1], 'MPI-ESM-LR'],
                       [sea2_mdl26.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl26)[0,1], 'MPI-ESM-MR'],
                       [sea2_mdl27.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl27)[0,1], 'MRI-CGCM3'],
                       [sea2_mdl28.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl28)[0,1], 'NCAR-CCSM4'],
                       [sea2_mdl29.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl29)[0,1], 'NCAR-CESM1-BGC'],
                       [sea2_mdl30.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl30)[0,1], 'NCAR-CESM1-CAM5'],
                       [sea2_mdl31.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl31)[0,1], 'NorESM1-ME'],
                       [sea2_mdl32.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl32)[0,1], 'NorESM1-M'],
                       [sea2_mdl33.std(ddof=1), np.corrcoef(sea2_obs1, sea2_mdl33)[0,1], 'ensmean_cmip5']],
               JJA=[[sea3_mdl1.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl1)[0,1], 'BCC-CSM1.1'],
                       [sea3_mdl2.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl2)[0,1], 'BCC-CSM1.1M'],
                       [sea3_mdl3.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl3)[0,1], 'BNU-ESM'],
                       [sea3_mdl4.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl4)[0,1], 'CanESM2'],
                       [sea3_mdl5.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl5)[0,1], 'CNRM-CM5'],
                       [sea3_mdl6.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl6)[0,1], 'CSIRO-ACCESS-1'],
                       [sea3_mdl7.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl7)[0,1], 'CSIRO-ACCESS-3'],
                       [sea3_mdl8.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl8)[0,1], 'CSIRO-MK36'],
                       [sea3_mdl9.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl9)[0,1], 'FIO-ESM'],
                       [sea3_mdl10.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl10)[0,1], 'GISS-E2-H-CC'],
                       [sea3_mdl11.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl11)[0,1], 'GISS-E2-H'],
                       [sea3_mdl12.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl12)[0,1], 'GISS-E2-R'],
                       [sea3_mdl13.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl13)[0,1], 'HadGEM2-AO'],
                       [sea3_mdl14.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl14)[0,1], 'HadGEM2-CC'],
                       [sea3_mdl15.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl15)[0,1], 'HadGEM2-ES'],
                       [sea3_mdl16.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl16)[0,1], 'INMCM4'],
                       [sea3_mdl17.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl17)[0,1], 'IPSL-CM5A-LR'],
                       [sea3_mdl18.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl18)[0,1], 'IPSL-CM5A-MR'],
                       [sea3_mdl19.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl19)[0,1], 'IPSL-CM5B-LR'],
                       [sea3_mdl20.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl20)[0,1], 'LASG-FGOALS-G2'],
                       [sea3_mdl21.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl21)[0,1], 'LASG-FGOALS-S2'],
                       [sea3_mdl22.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl22)[0,1], 'MIROC5'],
                       [sea3_mdl23.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl23)[0,1], 'MIROC-ESM-CHEM'],
                       [sea3_mdl24.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl24)[0,1], 'MIROC-ESM'],
                       [sea3_mdl25.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl25)[0,1], 'MPI-ESM-LR'],
                       [sea3_mdl26.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl26)[0,1], 'MPI-ESM-MR'],
                       [sea3_mdl27.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl27)[0,1], 'MRI-CGCM3'],
                       [sea3_mdl28.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl28)[0,1], 'NCAR-CCSM4'],
                       [sea3_mdl29.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl29)[0,1], 'NCAR-CESM1-BGC'],
                       [sea3_mdl30.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl30)[0,1], 'NCAR-CESM1-CAM5'],
                       [sea3_mdl31.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl31)[0,1], 'NorESM1-ME'],
                       [sea3_mdl32.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl32)[0,1], 'NorESM1-M'],
                       [sea3_mdl33.std(ddof=1), np.corrcoef(sea3_obs1, sea3_mdl33)[0,1], 'ensmean_cmip5']],
               SON=[[sea4_mdl1.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl1)[0,1], 'BCC-CSM1.1'],
                       [sea4_mdl2.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl2)[0,1], 'BCC-CSM1.1M'],
                       [sea4_mdl3.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl3)[0,1], 'BNU-ESM'],
                       [sea4_mdl4.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl4)[0,1], 'CanESM2'],
                       [sea4_mdl5.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl5)[0,1], 'CNRM-CM5'],
                       [sea4_mdl6.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl6)[0,1], 'CSIRO-ACCESS-1'],
                       [sea4_mdl7.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl7)[0,1], 'CSIRO-ACCESS-3'],
                       [sea4_mdl8.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl8)[0,1], 'CSIRO-MK36'],
                       [sea4_mdl9.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl9)[0,1], 'FIO-ESM'],
                       [sea4_mdl10.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl10)[0,1], 'GISS-E2-H-CC'],
                       [sea4_mdl11.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl11)[0,1], 'GISS-E2-H'],
                       [sea4_mdl12.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl12)[0,1], 'GISS-E2-R'],
                       [sea4_mdl13.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl13)[0,1], 'HadGEM2-AO'],
                       [sea4_mdl14.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl14)[0,1], 'HadGEM2-CC'],
                       [sea4_mdl15.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl15)[0,1], 'HadGEM2-ES'],
                       [sea4_mdl16.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl16)[0,1], 'INMCM4'],
                       [sea4_mdl17.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl17)[0,1], 'IPSL-CM5A-LR'],
                       [sea4_mdl18.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl18)[0,1], 'IPSL-CM5A-MR'],
                       [sea4_mdl19.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl19)[0,1], 'IPSL-CM5B-LR'],
                       [sea4_mdl20.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl20)[0,1], 'LASG-FGOALS-G2'],
                       [sea4_mdl21.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl21)[0,1], 'LASG-FGOALS-S2'],
                       [sea4_mdl22.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl22)[0,1], 'MIROC5'],
                       [sea4_mdl23.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl23)[0,1], 'MIROC-ESM-CHEM'],
                       [sea4_mdl24.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl24)[0,1], 'MIROC-ESM'],
                       [sea4_mdl25.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl25)[0,1], 'MPI-ESM-LR'],
                       [sea4_mdl26.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl26)[0,1], 'MPI-ESM-MR'],
                       [sea4_mdl27.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl27)[0,1], 'MRI-CGCM3'],
                       [sea4_mdl28.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl28)[0,1], 'NCAR-CCSM4'],
                       [sea4_mdl29.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl29)[0,1], 'NCAR-CESM1-BGC'],
                       [sea4_mdl30.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl30)[0,1], 'NCAR-CESM1-CAM5'],
                       [sea4_mdl31.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl31)[0,1], 'NorESM1-ME'],
                       [sea4_mdl32.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl32)[0,1], 'NorESM1-M'],
                       [sea4_mdl33.std(ddof=1), np.corrcoef(sea4_obs1, sea4_mdl33)[0,1], 'ensmean_cmip5']])

# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
colors = plt.matplotlib.cm.Set1(np.linspace(0,1,len(samples['DJF'])))

# Here set placement of the points marking 95th and 99th significance
# levels. For more than 102 samples (degrees freedom > 100), critical
# correlation levels are 0.195 and 0.254 for 95th and 99th
# significance levels respectively. Set these by eyeball using the
# standard deviation x and y axis.

#~ x95 = [0.01, 0.68] # For Tair, this is for 95th level (r = 0.195)
#~ y95 = [0.0, 3.45]
#~ x99 = [0.01, 0.95] # For Tair, this is for 99th level (r = 0.254)
#~ y99 = [0.0, 3.45]

x95 = [0.05, 13.9] # For Prcp, this is for 95th level (r = 0.195)
y95 = [0.0, 71.0]
x99 = [0.05, 19.0] # For Prcp, this is for 99th level (r = 0.254)
y99 = [0.0, 70.0]

rects = dict(DJF=221,
             MAM=222,
             JJA=223,
             SON=224)

# Plot model end obs data taylor diagram 
var_name  = u'Temperature 2m ($^\circ$C)' # Rainfall (mm/day) or Temperature 2m ($^\circ$C)
area_name = u'NEB (Lat:15S 2N, Lon:46W 34W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

out_var   = u'tmp' # pre or tmp
out_area  = u'neb' # amz or neb

fig = plt.figure(figsize=(20,12))
fig.suptitle('{0} Taylor Diagram - {1} \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Reference period: 1850-2005)'.format(var_name, area_name), fontsize=14, y=0.98)

for season in ['DJF','MAM','JJA','SON']:

    dia = TaylorDiagram(stdrefs[season], fig=fig, rect=rects[season],
                        label='Reference')

    dia.ax.plot(x95,y95,color='k')
    dia.ax.plot(x99,y99,color='k')

    # Add samples to Taylor diagram
    for i,(stddev,corrcoef,name) in enumerate(samples[season]):
        dia.add_sample(stddev, corrcoef,
                       marker='$%d$' % (i+1), ms=10, ls='',
                       #mfc='k', mec='k', # B&W
                       mfc=colors[i], mec=colors[i], # Colors
                       label=name)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=14, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title(season)

# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
# Can also use special options here:
# http://matplotlib.sourceforge.net/users/legend_guide.html

fig.legend(dia.samplePoints,
           [ p.get_label() for p in dia.samplePoints ],
           numpoints=1, prop=dict(size=14), loc='center')

fig.tight_layout()
path_out = '/home/nice'
name_out = 'pyplt_taylor_diagram_cmip5_cru_season_1975-2005.png'

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')




