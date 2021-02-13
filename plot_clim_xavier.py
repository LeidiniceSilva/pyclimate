# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/15/2020"
__description__ = "This script plot climatology graphics from Xavier obs basedata"

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# alto_parnaiba
est1_tmax = [30.512474,30.419691,30.464743,31.096563,31.557661,31.717657,32.104664,33.53589,34.60768,33.70436,31.813725,30.606049]
est1_tmin = [20.03944,20.060205,20.17641,20.147396,19.30523,17.574205,16.72776,17.565334,19.677149,20.611366,20.423203,20.248539]
est1_prec = [219.55652,182.48618,206.84721,107.33055,31.36229,1.9156437,1.3894821,3.8520641,18.733963,68.866394,156.14052,201.20578]

# balsas
est2_tmax = [30.820312,31.083712,31.033165,31.727396,32.34428,32.890835,33.320843,34.742237,35.20422,33.870186,31.98578,30.877369]
est2_tmin = [21.608341,21.653595,21.770237,21.808542,21.000807,19.209322,18.262135,18.984287,21.131094,22.05111,21.993933,21.85635]
est2_prec = [198.82523,155.13885,207.85786,97.01478,41.52467,5.105381,3.0888467,4.108902,23.001392,72.98577,137.5507,177.23888]

# barra_corda
est3_tmax = [31.214441,30.904812,30.754536,31.222109,31.590399,32.084507,32.875404,34.365803,35.126015,34.734375,33.67018,32.46026]
est3_tmin = [22.078201,22.056028,22.130192,22.278881,21.852219,20.67539,19.902672,20.378313,21.775677,22.68188,22.819975,22.61303]
est3_prec = [178.8049,175.72395,276.55304,180.98825,70.6456,20.617634,9.434916,6.207123,15.184965,37.727776,72.05738,118.49251]

# porto_nacional
est4_tmax = [31.512197,31.4185,31.467314,32.51237,33.302494,33.757137,34.336014,35.895138,35.964584,34.08846,32.44547,31.474396]
est4_tmin = [21.885887,21.881914,21.997278,22.124687,21.352219,19.328972,18.333317,19.352028,21.55547,22.162878,22.109192,22.056553]
est4_prec = [267.98273,247.28223,264.62607,159.36533,50.270363,2.6338909,1.3330654,4.4707966,43.07609,127.29675,216.02673,303.61533]

# pedro_afonso
est5_tmax = [31.216003,31.246046,31.075605,32.1168,32.871548,33.657734,34.453,35.761692,35.405834,33.25963,32.08917,31.38682]
est5_tmin = [22.219152,22.273506,22.364315,22.490026,21.841885,19.931406,18.937273,19.828691,21.825808,22.307737,22.36323,22.360811]
est5_prec = [262.36566,243.22185,271.1773,167.7032,55.853302,6.4482484,3.98131,8.699883,53.136467,144.9924,223.14168,263.95172]

# peixe
est6_tmax = [31.386164,31.539547,31.42558,32.273308,32.53891,32.573986,33.198387,35.041836,35.858334,34.474117,32.426952,31.211542]
est6_tmin = [22.13377,22.067284,22.264112,22.399452,21.421219,19.51418,18.704197,19.874294,22.068697,22.666481,22.424454,22.278402]
est6_prec = [260.42596,206.57983,239.53337,115.99778,27.236053,2.3830948,1.051042,2.8281555,30.804892,98.92752,199.05702,293.29327]

# taguatinga
est7_tmax = [29.365828,29.434458,29.391533,30.25888,30.295135,29.87453,29.802067,31.38256,32.768124,32.15237,30.164896,29.240927]
est7_tmin = [18.702608,18.667368,18.879839,18.79694,17.53159,15.554284,14.712828,15.910724,18.291822,19.438482,19.146835,18.974257]
est7_prec = [258.75034,249.04251,231.0055,114.42116,37.47221,2.261372,1.4909382,3.1765301,24.062105,91.354485,227.24345,291.57584]

# floriano
est8_tmax = [31.186592,30.948452,30.828201,31.33198,31.661718,32.138515,32.89569,34.43581,35.662422,35.541935,34.160233,32.810863]
est8_tmin = [21.699648,21.63194,21.702192,21.726171,21.244833,20.19961,19.605822,20.118359,21.706146,22.62132,22.65185,22.281855]
est8_prec = [196.63266,167.78638,237.83398,163.21823,55.185394,8.41461,5.973451,2.1534371,14.528684,34.11431,85.701035,127.15298]

# bom_jesus
est9_tmax = [30.982409,31.224585,31.089163,31.856874,32.38407,32.94302,33.12432,34.212376,35.46185,35.17258,33.521225,32.012726]
est9_tmin = [19.888306,9.87912,20.02379,19.975964,19.188772,17.732813,17.097897,17.781956,19.719284,20.73367,20.560104,20.23755]
est9_prec = [160.73666,129.643,138.19699,87.99954,27.156763,5.2094536,1.1332346,1.2621288,14.759476,47.72509,107.563416,137.15103]

# teresina
est10_tmax = [32.29791,31.870077,31.58309,31.682917,31.980444,32.275288,33.00882,34.724346,36.10357,36.448387,35.535572,34.40159]
est10_tmin = [22.753529,22.74989,22.69128,22.859453,22.576084,21.651302,20.995388,21.25015,22.284687,23.084703,23.37086,23.3312]
est10_prec = [200.00313,216.34497,319.9117,240.34985,99.924576,26.451794,12.884451,6.6968446,13.262314,25.874989,51.878113,106.60359]

# barreiras
est11_tmax = [31.414793,31.534153,31.335836,32.081173,32.206654,31.928072,32.10083,33.601788,35.111485,35.03241,32.642292,31.295893]
est11_tmin = [20.848362,20.839214,21.010157,20.666588,19.242603,17.274088,16.327597,17.42054,19.90022,21.463963,21.2744,21.085232]
est11_prec = [160.98077,135.73529,143.23935,58.42077,18.073734,1.1801082,0.6987953,2.2485366,14.899374,67.98729,173.84201,222.8271]

# correntina
est12_tmax = [31.891027,32.534348,32.121395,32.515392,32.145767,31.372187,31.22578,32.552395,34.253178,34.90872,32.697033,31.629662]
est12_tmin = [21.032082,21.2125,21.19201,20.695312,18.989529,17.158255,16.5454,17.427923,19.559114,21.236996,21.411797,21.227318]
est12_prec = [132.9362,104.04408,125.650856,41.02611,9.44953,2.1881757,2.8906128,2.386552,12.917013,45.84664,140.87846,199.52654]

# Define figure and figure size figsize=(width, height)
fig = plt.figure()
time = np.arange(1, 12+1)

ax1 = fig.add_subplot(4,3,1) 
ax2 = fig.add_subplot(4,3,2)
ax3 = fig.add_subplot(4,3,3) 
ax4 = fig.add_subplot(4,3,4) 
ax5 = fig.add_subplot(4,3,5) 
ax6 = fig.add_subplot(4,3,6)
ax7 = fig.add_subplot(4,3,7) 
ax8 = fig.add_subplot(4,3,8) 
ax9 = fig.add_subplot(4,3,9) 
ax10 = fig.add_subplot(4,3,10)
ax11 = fig.add_subplot(4,3,11) 
ax12 = fig.add_subplot(4,3,12)

#plot data and normalized data
ax1.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax1_1 = ax1.twinx()
tmp_clim1 = ax1_1.plot(time, est1_tmax, time, est1_tmin)

ax2.bar(time, est2_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax2_2 = ax2.twinx()
tmp_clim2 = ax2_2.plot(time, est2_tmax, time, est2_tmin)

ax3.bar(time, est3_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax3_3 = ax3.twinx()
tmp_clim3 = ax3_3.plot(time, est3_tmax, time, est3_tmin)

ax4.bar(time, est4_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax4_4 = ax4.twinx()
tmp_clim4 = ax4_4.plot(time, est4_tmax, time, est4_tmin)

ax5.bar(time, est5_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax5_5 = ax5.twinx()
tmp_clim5 = ax5_5.plot(time, est5_tmax, time, est5_tmin)

ax6.bar(time, est6_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax6_6 = ax6.twinx()
tmp_clim6 = ax6_6.plot(time, est6_tmax, time, est6_tmin)

ax7.bar(time, est7_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax7_7 = ax7.twinx()
tmp_clim7 = ax7_7.plot(time, est7_tmax, time, est7_tmin)

ax8.bar(time, est8_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax8_8 = ax8.twinx()
tmp_clim8 = ax8_8.plot(time, est8_tmax, time, est8_tmin)

ax9.bar(time, est9_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax9_9 = ax9.twinx()
tmp_clim9 = ax9_9.plot(time, est9_tmax, time, est9_tmin)

ax10.bar(time, est10_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax10_10 = ax10.twinx()
tmp_clim10 = ax10_10.plot(time, est10_tmax, time, est1_tmin)

ax11.bar(time, est11_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax11_11 = ax11.twinx()
tmp_clim11 = ax11_11.plot(time, est11_tmax, time, est11_tmin)

prec_clim = ax12.bar(time, est12_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax12_12 = ax12.twinx()
tmp_clim12 = ax12_12.plot(time, est12_tmax, time, est12_tmin)

l1_1, l2_1 = tmp_clim1
l1_2, l2_2 = tmp_clim2
l1_3, l2_3 = tmp_clim3
l1_4, l2_4 = tmp_clim4
l1_5, l2_5 = tmp_clim5
l1_6, l2_6 = tmp_clim6
l1_7, l2_7 = tmp_clim7
l1_8, l2_8 = tmp_clim8
l1_9, l2_9 = tmp_clim9
l1_10, l2_10 = tmp_clim10
l1_11, l2_11 = tmp_clim11
l1_12, l2_12 = tmp_clim12

plt.setp(l1_1, color='red')
plt.setp(l2_1, color='blue')
plt.setp(l1_2, color='red')
plt.setp(l2_2, color='blue')
plt.setp(l1_3, color='red')
plt.setp(l2_3, color='blue')
plt.setp(l1_4, color='red')
plt.setp(l2_4, color='blue')
plt.setp(l1_5, color='red')
plt.setp(l2_5, color='blue')
plt.setp(l1_6, color='red')
plt.setp(l2_6, color='blue')
plt.setp(l1_7, color='red')
plt.setp(l2_7, color='blue')
plt.setp(l1_8, color='red')
plt.setp(l2_8, color='blue')
plt.setp(l1_9, color='red')
plt.setp(l2_9, color='blue')
plt.setp(l1_10, color='red')
plt.setp(l2_10, color='blue')
plt.setp(l1_11, color='red')
plt.setp(l2_11, color='blue')
plt.setp(l1_12, color='red')
plt.setp(l2_12, color='blue')

# plots 1 2 3
ax1.set_ylim(0, 400)
ax2.set_ylim(0, 400)
ax3.set_ylim(0, 400)

ax1.set_yticks(np.arange(0, 400, 80))
ax2.set_yticks(np.arange(0, 400, 80))
ax3.set_yticks(np.arange(0, 400, 80))

ax2.set_yticklabels([])
ax3.set_yticklabels([])

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])

ax1_1.set_ylim(0, 50)
ax2_2.set_ylim(0, 50)
ax3_3.set_ylim(0, 50)

ax1_1.set_yticks(np.arange(0, 50, 10))
ax2_2.set_yticks(np.arange(0, 50, 10))
ax3_3.set_yticks(np.arange(0, 50, 10))

ax1_1.set_yticklabels([])
ax2_2.set_yticklabels([])

# Plots 4 5 6
ax4.set_ylim(0, 400)
ax5.set_ylim(0, 400)
ax6.set_ylim(0, 400)

ax4.set_yticks(np.arange(0, 400, 80))
ax5.set_yticks(np.arange(0, 400, 80))
ax6.set_yticks(np.arange(0, 400, 80))

ax5.set_yticklabels([])
ax6.set_yticklabels([])

ax4.set_xticklabels([])
ax5.set_xticklabels([])
ax6.set_xticklabels([])

ax4_4.set_ylim(0, 50)
ax5_5.set_ylim(0, 50)
ax6_6.set_ylim(0, 50)

ax4_4.set_yticks(np.arange(0, 50, 10))
ax5_5.set_yticks(np.arange(0, 50, 10))
ax6_6.set_yticks(np.arange(0, 50, 10))

ax4_4.set_yticklabels([])
ax5_5.set_yticklabels([])

# Plots 7 8 9
ax7.set_ylim(0, 400)
ax8.set_ylim(0, 400)
ax9.set_ylim(0, 400)

ax7.set_yticks(np.arange(0, 400, 80))
ax8.set_yticks(np.arange(0, 400, 80))
ax9.set_yticks(np.arange(0, 400, 80))

ax8.set_yticklabels([])
ax9.set_yticklabels([])

ax7.set_xticklabels([])
ax8.set_xticklabels([])
ax9.set_xticklabels([])

ax7_7.set_ylim(0, 50)
ax8_8.set_ylim(0, 50)
ax9_9.set_ylim(0, 50)

ax7_7.set_yticks(np.arange(0, 50, 10))
ax8_8.set_yticks(np.arange(0, 50, 10))
ax9_9.set_yticks(np.arange(0, 50, 10))

ax7_7.set_yticklabels([])
ax8_8.set_yticklabels([])

# Plots 10, 11, 12
ax10.set_ylim(0, 400)
ax11.set_ylim(0, 400)
ax12.set_ylim(0, 400)

ax10.set_yticks(np.arange(0, 400, 80))
ax11.set_yticks(np.arange(0, 400, 80))
ax12.set_yticks(np.arange(0, 400, 80))

ax11.set_yticklabels([])
ax12.set_yticklabels([])

ax10_10.set_ylim(0, 50)
ax11_11.set_ylim(0, 50)
ax12_12.set_ylim(0, 50)

ax10_10.set_yticks(np.arange(0, 50, 10))
ax11_11.set_yticks(np.arange(0, 50, 10))
ax12_12.set_yticks(np.arange(0, 50, 10))

ax10_10.set_yticklabels([])
ax11_11.set_yticklabels([])

# Define the location of ticks (yAxis 1)
ax1.yaxis.tick_left()
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.tick_params(axis='both', which='both', labelsize=10, length=0)
ax3.yaxis.tick_right()

ax4.yaxis.tick_left()
plt.setp(ax5.get_yticklabels(), visible=False)
ax5.tick_params(axis='both', which='both', labelsize=10, length=0)
ax6.yaxis.tick_right()

ax7.yaxis.tick_left()
plt.setp(ax8.get_yticklabels(), visible=False)
ax8.tick_params(axis='both', which='both', labelsize=10, length=0)
ax9.yaxis.tick_right()

ax10.yaxis.tick_left()
plt.setp(ax11.get_yticklabels(), visible=False)
ax11.tick_params(axis='y', which='both', labelsize=10, length=0)
ax12.yaxis.tick_right()

# Define the location of ticks (yAxis 2)
plt.setp(ax1_1.get_yticklabels(), visible=False)
ax1_1.tick_params(axis='both', which='both', length=0)

plt.setp(ax2_2.get_yticklabels(), visible=False)
ax2_2.tick_params(axis='both', which='both', length=0)

ax3_3.tick_params(axis='y', which='both', labelsize=10, length=0)

plt.setp(ax4_4.get_yticklabels(), visible=False)
ax4_4.tick_params(axis='both', which='both', length=0)

plt.setp(ax5_5.get_yticklabels(), visible=False)
ax5_5.tick_params(axis='both', which='both', length=0)

ax6_6.tick_params(axis='y', which='both', labelsize=10, length=0)

plt.setp(ax7_7.get_yticklabels(), visible=False)
ax7_7.tick_params(axis='both', which='both', length=0)

plt.setp(ax8_8.get_yticklabels(), visible=False)
ax8_8.tick_params(axis='both', which='both', length=0)

ax9_9.tick_params(axis='y', which='both', labelsize=10, length=0)

plt.setp(ax10_10.get_yticklabels(), visible=False)
ax10_10.tick_params(axis='both', which='both', length=0)

plt.setp(ax11_11.get_yticklabels(), visible=False)
ax11_11.tick_params(axis='both', which='both', length=0)

ax12_12.tick_params(axis='y', which='both', labelsize=10, length=0)

# Define the location of ticks (xAxis)
ax10.set_xticks(np.arange(1, 12+1, 1))
ax10.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
plt.setp(ax10.get_xticklabels(),fontsize=10)

ax11.set_xticks(np.arange(1, 12+1, 1))
ax11.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
plt.setp(ax11.get_xticklabels(),fontsize=10)

ax12.set_xticks(np.arange(1, 12+1, 1))
ax12.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
plt.setp(ax12.get_xticklabels(),fontsize=10)

fig.suptitle('Ciclo Anual - Estações: MATOPIBA (1986-2005)', y=1.01)

myl=[prec_clim]+tmp_clim1
legend = (u'Precipitação (mm)', u'Temperatura máx ($^\circ$C)', u'Temperatura min ($^\circ$C)')    
ax11.legend(myl, legend, loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=True, ncol=3, prop=FontProperties(size=8))

ax1.text(0.5, 300, u'A) Alto Parnaíba', fontsize=8)
ax2.text(0.5, 300, u'B) Balsas', fontsize=8)
ax3.text(0.5, 300, u'C) Barra do Corda', fontsize=8)
ax4.text(0.5, 300, u'D) Porto Nacional', fontsize=8)
ax5.text(0.5, 300, u'E) Pedro Afonso', fontsize=8)
ax6.text(0.5, 300, u'F) Peixe', fontsize=8)
ax7.text(0.5, 300, u'G) Taguatinga', fontsize=8)
ax8.text(0.5, 300, u'H) Floriano', fontsize=8)
ax9.text(0.5, 300, u'I) Bom Jesus', fontsize=8)
ax10.text(0.5, 300, u'J) Teresina', fontsize=8)
ax11.text(0.5, 300, u'L) Barreiras', fontsize=8)
ax12.text(0.5, 300, u'M) Correntina', fontsize=8)

plt.subplots_adjust(left=0.15, bottom=0.06, right=0.93, top=0.97, wspace=0.05, hspace=0)

# Path out to save figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_clim_pr_tmp_xavier_matopiba_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
