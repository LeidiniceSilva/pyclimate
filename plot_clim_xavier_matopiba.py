# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology graphics from Xavier obs basedata"

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

#data
est1_tmax = [31,33,31,34,30,31,33,31,32,32,31,33]
est1_tmin = [21,23,21,24,20,21,23,21,22,22,21,23]
est1_prec = [211,232,213,244,200,110,139,111, 119,102,119,236]

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

ax2.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax2_2 = ax2.twinx()
tmp_clim2 = ax2_2.plot(time, est1_tmax, time, est1_tmin)

ax3.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax3_3 = ax3.twinx()
tmp_clim3 = ax3_3.plot(time, est1_tmax, time, est1_tmin)

ax4.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax4_4 = ax4.twinx()
tmp_clim4 = ax4_4.plot(time, est1_tmax, time, est1_tmin)

ax5.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax5_5 = ax5.twinx()
tmp_clim5 = ax5_5.plot(time, est1_tmax, time, est1_tmin)

ax6.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax6_6 = ax6.twinx()
tmp_clim6 = ax6_6.plot(time, est1_tmax, time, est1_tmin)

ax7.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax7_7 = ax7.twinx()
tmp_clim7 = ax7_7.plot(time, est1_tmax, time, est1_tmin)

ax8.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax8_8 = ax8.twinx()
tmp_clim8 = ax8_8.plot(time, est1_tmax, time, est1_tmin)

ax9.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax9_9 = ax9.twinx()
tmp_clim9 = ax9_9.plot(time, est1_tmax, time, est1_tmin)

ax10.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax10_10 = ax10.twinx()
tmp_clim10 = ax10_10.plot(time, est1_tmax, time, est1_tmin)

ax11.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax11_11 = ax11.twinx()
tmp_clim11 = ax11_11.plot(time, est1_tmax, time, est1_tmin)

prec_clim = ax12.bar(time, est1_prec, alpha=0.8, color='gray', width = 0.40, edgecolor='black')
ax12_12 = ax12.twinx()
tmp_clim12 = ax12_12.plot(time, est1_tmax, time, est1_tmin)

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
ax1.set_ylim(0, 500)
ax2.set_ylim(0, 500)
ax3.set_ylim(0, 500)

ax1.set_yticks(np.arange(0, 500, 100))
ax2.set_yticks(np.arange(0, 500, 100))
ax3.set_yticks(np.arange(0, 500, 100))

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
ax4.set_ylim(0, 500)
ax5.set_ylim(0, 500)
ax6.set_ylim(0, 500)

ax4.set_yticks(np.arange(0, 500, 100))
ax5.set_yticks(np.arange(0, 500, 100))
ax6.set_yticks(np.arange(0, 500, 100))

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
ax7.set_ylim(0, 500)
ax8.set_ylim(0, 500)
ax9.set_ylim(0, 500)

ax7.set_yticks(np.arange(0, 500, 100))
ax8.set_yticks(np.arange(0, 500, 100))
ax9.set_yticks(np.arange(0, 500, 100))

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
ax10.set_ylim(0, 500)
ax11.set_ylim(0, 500)
ax12.set_ylim(0, 500)

ax10.set_yticks(np.arange(0, 500, 100))
ax11.set_yticks(np.arange(0, 500, 100))
ax12.set_yticks(np.arange(0, 500, 100))

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

ax1.text(0.5, 400, u'A) Alto Parnaíba', fontsize=8)
ax2.text(0.5, 400, u'B) Balsas', fontsize=8)
ax3.text(0.5, 400, u'C) Barra do Corda', fontsize=8)
ax4.text(0.5, 400, u'D) Porto Nacional', fontsize=8)
ax5.text(0.5, 400, u'E) Pedro Afonso', fontsize=8)
ax6.text(0.5, 400, u'F) Peixe', fontsize=8)
ax7.text(0.5, 400, u'G) Taguatinga', fontsize=8)
ax8.text(0.5, 400, u'H) Floriano', fontsize=8)
ax9.text(0.5, 400, u'I) Bom Jesus', fontsize=8)
ax10.text(0.5, 400, u'J) Teresina', fontsize=8)
ax11.text(0.5, 400, u'L) Barreiras', fontsize=8)
ax12.text(0.5, 400, u'M) Correntina', fontsize=8)

plt.subplots_adjust(left=0.15, bottom=0.06, right=0.93, top=0.97, wspace=0.05, hspace=0)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_clim_pr_tmp_xavier_matopiba_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
