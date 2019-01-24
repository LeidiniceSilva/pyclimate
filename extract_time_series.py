# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/21/2019"
__description__ = "This script to extract time series from CMIP5 models"


import netCDF4
import numpy as np
import matplotlib.pyplot as plt
# mpl.use('Agg')


var   = 'pr'
model = 'BCC-CSM1.1'
exp   = 'historical_r1i1p1'
data  = '185001-201212'

mod_path = '/home/nice/Documentos/cmip/cmip5_hist/BCC-CSM1.1'
arq1     = '{0}/{1}_Amon_{2}_{3}_{4}.nc'.format(mod_path, var, model, exp, data)
data1    = netCDF4.Dataset(arq1)

var1     = data1.variables['pr'][:] * 86400
lat      = data1.variables['lat'][:]
lon      = data1.variables['lon'][:]
mod      = var1[:][:,:,:]
mod_ini1 = np.nanmean(mod, axis=1)
mod_end1 = np.nanmean(mod_ini1, axis=1)

#~ print mod_end1
#~ print
#~ exit()

# plot pluma graphs of all models
fig  = plt.figure(figsize=(12, 10))
time = np.arange(1, 1956 + 1)

#~ a = plt.plot(time, mod_end1)
#~ plt.show()
#~ exit()

mod_end1.write(mod_end1)
mod_end1.close()
exit()


