# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/15/2020"
__description__ = "This script compute time series climatology from Xavier basedata"

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# pegando variavel
path_var = '/home/nice/Documents/'
ds = xr.open_mfdataset(path_var + 'Tmin_daily_UT_Brazil_v2.1_19860101_20051231.nc')
var = ds['Tmin']

# Nome dos pontos
Names = ['alto_parnaiba', 'balsas', 'barra_corda', 'porto_nacional', 'pedro_afonso', 'peixe', 
'taguatinga', 'floriano', 'bom_jesus', 'teresina', 'barreiras', 'correntina']

lat_lon = [[-9.06, -45.56],
         [-7.32, -46.02],
         [-5.30, -45.13],
         [-10.43, -48.25],
         [-8.96, -48.18],
         [-12.00, -48.22],
         [-12.24, -46.25],
         [-6.46, -43.01],
         [-9.06, -44.07],
         [-5.08, -42.81],
         [-12.09, -45.00],
         [-13.20, -43.37]]
         
varMonthly2Export = pd.DataFrame(np.empty((12, len(Names))),
                                 columns=Names, index=range(1, 13))
# media mensal
for n, names in enumerate(Names):
    CityDaily = var.sel(latitude=lat_lon[n][0], longitude=lat_lon[n][1], method='nearest')
    CityMonthly = CityDaily.groupby('time.month').mean('time')
	
	# concatenating in pandas
    varMonthly2Export[names] = CityMonthly

# exporting
varMonthly2Export.to_csv('tmin_xavier_clim.csv')

