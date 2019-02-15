# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "02/15/2019"
__description__ = "This script plot Taylor Diagram seasonal from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from plot_taylor_diagram_cmip5_obs import TaylorDiagram


def import_cmip5_clim(model):
	
	param = 'pr' # pr or tas
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
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

	return sea1_mdl, sea2_mdl, sea2_mdl, sea4_mdl


def import_obs_clim(database):
	
	param = 'pre' # pre or tmp
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
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


model  = u'BCC-CSM1.1'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1 = import_cmip5_clim(model)

model  = u'BCC-CSM1.1M'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'BNU-ESM'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'CanESM2'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'CNRM-CM5'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-1'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-3'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'CSIRO-MK36'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'FIO-ESM'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'GISS-E2-H-CC'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'GISS-E2-H'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'GISS-E2-R'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'HadGEM2-AO'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'HadGEM2-CC'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'HadGEM2-ES'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'INMCM4'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'IPSL-CM5A-LR'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'IPSL-CM5A-MR'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'IPSL-CM5B-LR'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'LASG-FGOALS-G2'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'LASG-FGOALS-S2'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MIROC5'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MIROC-ESM-CHEM'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MIROC-ESM'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MPI-ESM-LR'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MPI-ESM-MR'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'MRI-CGCM3'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'NCAR-CCSM4'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'NCAR-CESM1-BGC'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'NCAR-CESM1-CAM5'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'NorESM1-ME'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'NorESM1-M'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

model  = u'ensmean_cmip5'
sea1_mdl1, sea2_mdl1, sea2_mdl1, sea4_mdl1  = import_cmip5_clim(model)

print SON
exit()
	
# Reference std
stdrefs = dict(winter=48.491,
               spring=44.927,
               summer=37.664,
               autumn=41.589)

# Sample std,rho: Be sure to check order and that correct numbers are placed!
samples = dict(winter=[[17.831, 0.360, "CCSM CRCM"],
                       [27.062, 0.360, "CCSM MM5"],
                       [33.125, 0.585, "CCSM WRFG"],
                       [25.939, 0.385, "CGCM3 CRCM"],
                       [29.593, 0.509, "CGCM3 RCM3"],
                       [35.807, 0.609, "CGCM3 WRFG"],
                       [38.449, 0.342, "GFDL ECP2"],
                       [29.593, 0.509, "GFDL RCM3"],
                       [71.215, 0.473, "HADCM3 HRM3"]],
               spring=[[32.174, -0.262, "CCSM CRCM"],
                       [24.042, -0.055, "CCSM MM5"],
                       [29.647, -0.040, "CCSM WRFG"],
                       [22.820, 0.222, "CGCM3 CRCM"],
                       [20.505, 0.445, "CGCM3 RCM3"],
                       [26.917, 0.332, "CGCM3 WRFG"],
                       [25.776, 0.366, "GFDL ECP2"],
                       [18.018, 0.452, "GFDL RCM3"],
                       [79.875, 0.447, "HADCM3 HRM3"]],
               summer=[[35.863, 0.096, "CCSM CRCM"],
                       [43.771, 0.367, "CCSM MM5"],
                       [35.890, 0.267, "CCSM WRFG"],
                       [49.658, 0.134, "CGCM3 CRCM"],
                       [28.972, 0.027, "CGCM3 RCM3"],
                       [60.396, 0.191, "CGCM3 WRFG"],
                       [46.529, 0.258, "GFDL ECP2"],
                       [35.230, -0.014, "GFDL RCM3"],
                       [87.562, 0.503, "HADCM3 HRM3"]],
               autumn=[[27.374, 0.150, "CCSM CRCM"],
                       [20.270, 0.451, "CCSM MM5"],
                       [21.070, 0.505, "CCSM WRFG"],
                       [25.666, 0.517, "CGCM3 CRCM"],
                       [35.073, 0.205, "CGCM3 RCM3"],
                       [25.666, 0.517, "CGCM3 WRFG"],
                       [23.409, 0.353, "GFDL ECP2"],
                       [29.367, 0.235, "GFDL RCM3"],
                       [70.065, 0.444, "HADCM3 HRM3"]])

# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
colors = plt.matplotlib.cm.Set1(np.linspace(0,1,len(samples['winter'])))

# Here set placement of the points marking 95th and 99th significance
# levels. For more than 102 samples (degrees freedom > 100), critical
# correlation levels are 0.195 and 0.254 for 95th and 99th
# significance levels respectively. Set these by eyeball using the
# standard deviation x and y axis.

#x95 = [0.01, 0.68] # For Tair, this is for 95th level (r = 0.195)
#y95 = [0.0, 3.45]
#x99 = [0.01, 0.95] # For Tair, this is for 99th level (r = 0.254)
#y99 = [0.0, 3.45]

x95 = [0.05, 13.9] # For Prcp, this is for 95th level (r = 0.195)
y95 = [0.0, 71.0]
x99 = [0.05, 19.0] # For Prcp, this is for 99th level (r = 0.254)
y99 = [0.0, 70.0]

rects = dict(winter=221,
             spring=222,
             summer=223,
             autumn=224)

fig = plt.figure(figsize=(11,8))
fig.suptitle("Rainfall Taylor Diagram", size='x-large')

for season in ['winter','spring','summer','autumn']:

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
    dia.ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title(season.capitalize())

# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
# Can also use special options here:
# http://matplotlib.sourceforge.net/users/legend_guide.html

fig.legend(dia.samplePoints,
           [ p.get_label() for p in dia.samplePoints ],
           numpoints=1, prop=dict(size='small'), loc='center')

fig.tight_layout()

path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
name_out = 'pyplt_taylor_diagram_pre_amz_neb_cmip5_cru_season_1975-2005.png'

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out))
plt.show()




