import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as co

#data
out1 = [11,13,14,15,17,21,21,27]
out2 = [1236,1673,2203,2803,3288,3610,3675,4132]
out3 = [492,626,731,850,865,817,659,550]
out4 = [68,205,325,421,611,883,1343,1646]
out5 = [1612,2274,2894,3643,4206,4287,4934,5475]
out6 = [1346,1833,2383,3576,2919,3306,3517,3364]
out7 = [1252,1634,2149,2761,3071,3205,3088,3278]
out8 = [13595,18271,23755,29429,33968,36602,36706,39761]
out9 = [6341,9619,17873,15622,23714,18384,20295,21127]
out10 = [84,112.6,248.3936402,420.9011487,632,825,971.27,1343.7]
out11 = [374,535.1847201,821.9746154,1147.672,1526.072807,1857.583225,2165.713401,2712.648204]
I_out1 = [1305,1991.881,2671.393,3394.5,4093.24719,4706.72,5243.552,6016.643]
gas1 = [4186,7473,10356,11242,11475,12149,12149,12149]
gas2 = [0,0,0,7,62,601,342,396,]
gas3 = [3,7,7.8,13,16,19,20,23]
flux1 = [0,0,0,1.5001,10.35221696,12.4,25.934,33.727]
flux2 = [6,13.69,22.358,31.7604005,42.0262,47.8298,56.31445151,72.54693]
flux3 = [17,20.048,27.779,31.7394,31.211,29.0254371,23.9779,29.64697]
flux4 = [98,75,91,18,64,88,33,48]

#function to normalize data(right y axis)
def norm(y):
    normf = co.normalize()
    normr = normf(y)
    return normr


#define figure and figure size figsize=(width, height)
fig = plt.figure(figsize=(6.299, 9))

#define subplots 3x6
# 1        2    3
# 4        5    6
# 7        8    9
# 10       11   12



ax1 = fig.add_subplot(6,3,1) #gas1
ax1n = fig.add_subplot(6,3,1, sharex=ax1, frameon=False) #gas1 nomeeritud


ax2 = fig.add_subplot(6,3,2) #out9
ax2n = fig.add_subplot(6,3,2, sharex=ax2, frameon=False) #out9

ax3 = fig.add_subplot(6,3,3) #out8
ax3n = fig.add_subplot(6,3,3, sharex=ax3, frameon=False) #out8


ax4 = fig.add_subplot(6,3,4) #out2
ax4n = fig.add_subplot(6,3,4, sharex=ax4, frameon=False) #out2

ax5 = fig.add_subplot(6,3,5) #out5
ax5n = fig.add_subplot(6,3,5, sharex=ax5, frameon=False) #out5


ax6 = fig.add_subplot(6,3,6) #I_out1
ax6n = fig.add_subplot(6,3,6, sharex=ax6, frameon=False) #I_out1


ax7 = fig.add_subplot(6,3,7) #out11
ax7n = fig.add_subplot(6,3,7, sharex=ax7, frameon=False) #out11

ax8 = fig.add_subplot(6,3,8) #out7
ax8n = fig.add_subplot(6,3,8,sharex=ax8, frameon=False) #out7


ax9 = fig.add_subplot(6,3,9) #out6
ax9n = fig.add_subplot(6,3,9, sharex=ax9, frameon=False) #out6


ax10 = fig.add_subplot(6,3,10) #out3
ax10n = fig.add_subplot(6,3,10, sharex=ax10, frameon=False) #out3


ax11 = fig.add_subplot(6,3,11) #out10
ax11n = fig.add_subplot(6,3,11, sharex=ax11, frameon=False) #out10


ax12 = fig.add_subplot(6,3,12) #out4
ax12n = fig.add_subplot(6,3,12, sharex=ax12, frameon=False) #out4

ax13 = fig.add_subplot(6,3,13) #flux2
ax13n = fig.add_subplot(6,3,13, sharex=ax13, frameon=False) #flux2


ax14 = fig.add_subplot(6,3,14) #flux4
ax14n = fig.add_subplot(6,3,14, sharex=ax14, frameon=False) #flux4

ax15 = fig.add_subplot(6,3,15) #gas2
ax15n = fig.add_subplot(6,3,15, sharex=ax15, frameon=False) #gas2



ax16 = fig.add_subplot(6,3,16) #gas3
ax16n = fig.add_subplot(6,3,16, sharex=ax16, frameon=False) #gas3

ax17 = fig.add_subplot(6,3,17) #flux1
ax17n = fig.add_subplot(6,3,17, sharex=ax17, frameon=False) #flux1


ax18 = fig.add_subplot(6,3,18) #flux3
ax18n = fig.add_subplot(6,3,18, sharex=ax18, frameon=False) #flux3

#plot data and normalized data
ax1.plot(out1, gas1, ls="None", marker="o",mfc="None", color="k")
ax1n.plot(out1, norm(gas1), ls="dotted",  color="k")

ax2.plot(out1, out9, ls="None", marker="s", mfc="None",color="k")
ax2n.plot(out1, norm(out9), ls="dotted",color="k")

ax3.plot(out1, out8, ls="None", marker="d",mfc="None", color="k")
ax3n.plot(out1, norm(out8), ls="dotted",mfc="None", color="k")

ax4.plot(out1, I_out1, ls="None", marker="+", color="k")
ax4n.plot(out1, norm(I_out1), ls="dotted", color="k")

ax5.plot(out1, out5, ls="None", marker="s", color="k")
ax5n.plot(out1, norm(out5), ls="dotted", color="k")


ax6.plot(out1, out2, ls="None", marker="x", color="k")
ax6n.plot(out1, norm(out2), ls="dotted", color="k")

ax7.plot(out1, out11, ls="None", marker=">", color="k")
ax7n.plot(out1, norm(out11), ls="dotted",  color="k")

ax8.plot(out1, out7, ls="None", marker="<", color="k")
ax8n.plot(out1, norm(out7), ls="dotted",  color="k")


ax9.plot(out1, out6, ls="None", marker="^", color="k")
ax9n.plot(out1, norm(out6), ls="dotted", color="k")

ax10.plot(out1, out3, ls="None", marker=">",mfc="None", color="k")
ax10n.plot(out1, norm(out3), ls="dotted",mfc="None", color="k")

ax11.plot(out1, out10, ls="None", marker="<",mfc="None", color="k")
ax11n.plot(out1, norm(out10), ls="dotted", marker="<",mfc="None", color="k")

ax12.plot(out1, out4, ls="None", marker="^",mfc="None", color="k")
ax12n.plot(out1, norm(out4), ls="dotted", mfc="None", color="k")

ax13.plot(out1, flux2, ls="None", marker="o", color="k")
ax13n.plot(out1, norm(flux2), ls="dotted",  color="k")

ax14.plot(out1, flux4, ls="None", marker="d", color="k")
ax14n.plot(out1, norm(flux4), ls="dotted",  color="k")

ax15.plot(out1, gas2, ls="None", marker="*", color="k")
ax15n.plot(out1, norm(gas2), ls="dotted",  color="k")

ax16.plot(out1, gas3, ls="None", marker="_", color="k")
ax16n.plot(out1, norm(gas3), ls="dotted", color="k")

ax17.plot(out1, flux1, ls="None", marker="v", color="k")
ax17n.plot(out1, norm(flux1), ls="dotted", color="k")

ax18.plot(out1, flux3, ls="None", marker="v",mfc="None", color="k")
ax18n.plot(out1, norm(flux3), ls="dotted", mfc="None", color="k")



#configure axis

#123
ax1.set_ylim(3000, 50000)
ax2.set_ylim(3000, 50000)
ax3.set_ylim(3000, 50000)

ax1.set_yticks(np.arange(5000, 50000, 10000))
ax2.set_yticks(np.arange(5000, 50000, 10000))
ax3.set_yticks(np.arange(5000, 50000, 10000))

#hide Y tick labels for some plots(only plots on the left and right have labels and ticklabels
ax2.set_yticklabels([])
ax3.set_yticklabels([])

ax1.set_xlim(5, 35)
ax2.set_xlim(5, 35)
ax3.set_xlim(5, 35)

ax1.set_xticks(np.arange(10, 35, 5))
ax2.set_xticks(np.arange(10, 35, 5))
ax3.set_xticks(np.arange(10, 35, 5))



ax2.set_xticklabels([])
ax3.set_xticklabels([])

#456
ax4.set_ylim(500, 7500)
ax5.set_ylim(500, 7500)
ax6.set_ylim(500, 7500)

ax4.set_yticks(np.arange(1000, 8000, 1500))
ax5.set_yticks(np.arange(1000, 8000, 1500))
ax6.set_yticks(np.arange(1000, 8000, 1500))


ax5.set_yticklabels([])
ax6.set_yticklabels([])

ax4.set_xlim(5, 35)
ax5.set_xlim(5, 35)
ax6.set_xlim(5, 35)

ax4.set_xticks(np.arange(10, 35, 5))
ax5.set_xticks(np.arange(10, 35, 5))
ax6.set_xticks(np.arange(10, 35, 5))

ax4.set_xticklabels([])
ax5.set_xticklabels([])
ax6.set_xticklabels([])

#789
ax7.set_ylim(0, 4000)
ax8.set_ylim(0, 4000)
ax9.set_ylim(0, 4000)

ax7.set_yticks(np.arange(500, 4000, 1000))
ax8.set_yticks(np.arange(500, 4000, 1000))
ax9.set_yticks(np.arange(500, 4000, 1000))

ax8.set_yticklabels([])
ax9.set_yticklabels([])

ax7.set_xlim(5, 35)
ax8.set_xlim(5, 35)
ax9.set_xlim(5, 35)

ax7.set_xticks(np.arange(10, 35, 5))
ax8.set_xticks(np.arange(10, 35, 5))
ax9.set_xticks(np.arange(10, 35, 5))

ax7.set_xticklabels([])
ax8.set_xticklabels([])
ax9.set_xticklabels([])

# 10, 11, 12
ax10.set_ylim(-100, 2000)
ax11.set_ylim(-100, 2000)
ax12.set_ylim(-100, 2000)

ax10.set_yticks(np.arange(0, 1900, 300))
ax11.set_yticks(np.arange(0, 1900, 300))
ax12.set_yticks(np.arange(0, 1900, 300))

ax11.set_yticklabels([])
ax12.set_yticklabels([])

ax10.set_xlim(5, 35)
ax11.set_xlim(5, 35)
ax12.set_xlim(5, 35)

ax10.set_xticks(np.arange(10, 35, 5))
ax11.set_xticks(np.arange(10, 35, 5))
ax12.set_xticks(np.arange(10, 35, 5))

ax10.set_xticklabels([])
ax11.set_xticklabels([])
ax12.set_xticklabels([])

# 13, 14, 15
ax13.set_ylim(-50, 700)
ax14.set_ylim(-50, 700)
ax15.set_ylim(-50, 700)

ax13.set_yticks(np.arange(0, 750, 150))
ax14.set_yticks(np.arange(0, 750, 150))
ax15.set_yticks(np.arange(0, 750, 150))

ax14.set_yticklabels([])
ax15.set_yticklabels([])

ax13.set_xlim(5, 35)
ax14.set_xlim(5, 35)
ax15.set_xlim(5, 35)

ax13.set_xticks(np.arange(10, 35, 5))
ax14.set_xticks(np.arange(10, 35, 5))
ax15.set_xticks(np.arange(10, 35, 5))

ax13.set_xticklabels([])
ax14.set_xticklabels([])
ax15.set_xticklabels([])

# 16, 17, 18
ax16.set_ylim(-5, 45)
ax17.set_ylim(-5, 45)
ax18.set_ylim(-5, 45)

ax16.set_yticks(np.arange(0, 50, 10))
ax17.set_yticks(np.arange(0, 50, 10))
ax18.set_yticks(np.arange(0, 50, 10))

ax17.set_yticklabels([])
ax18.set_yticklabels([])

ax16.set_xlim(5, 35)
ax17.set_xlim(5, 35)
ax18.set_xlim(5, 35)

ax16.set_xticks(np.arange(10, 35, 5))
ax17.set_xticks(np.arange(10, 35, 5))
ax18.set_xticks(np.arange(10, 35, 5))



#disable labels for some plots
ax1n.set_yticklabels([])
ax2n.set_yticklabels([])

ax4n.set_yticklabels([])
ax5n.set_yticklabels([])

ax7n.set_yticklabels([])
ax8n.set_yticklabels([])

ax10n.set_yticklabels([])
ax11n.set_yticklabels([])

ax13n.set_yticklabels([])
ax14n.set_yticklabels([])

ax16n.set_yticklabels([])
ax17n.set_yticklabels([])


#define the location of ticks

ax1n.yaxis.tick_right()
ax2n.yaxis.tick_right()
ax3n.yaxis.tick_right()
ax4n.yaxis.tick_right()
ax5n.yaxis.tick_right()
ax6n.yaxis.tick_right()
ax7n.yaxis.tick_right()
ax8n.yaxis.tick_right()
ax9n.yaxis.tick_right()
ax10n.yaxis.tick_right()
ax11n.yaxis.tick_right()
ax12n.yaxis.tick_right()
ax13n.yaxis.tick_right()
ax14n.yaxis.tick_right()
ax15n.yaxis.tick_right()
ax16n.yaxis.tick_right()
ax17n.yaxis.tick_right()
ax18n.yaxis.tick_right()

ax1.yaxis.tick_left()
ax2.yaxis.tick_left()
ax3.yaxis.tick_left()
ax4.yaxis.tick_left()
ax5.yaxis.tick_left()
ax6.yaxis.tick_left()
ax7.yaxis.tick_left()
ax8.yaxis.tick_left()
ax9.yaxis.tick_left()
ax10.yaxis.tick_left()
ax11.yaxis.tick_left()
ax12.yaxis.tick_left()
ax13.yaxis.tick_left()
ax14.yaxis.tick_left()
ax15.yaxis.tick_left()
ax16.yaxis.tick_left()
ax17.yaxis.tick_left()
ax18.yaxis.tick_left()


ax1n.set_ylim(-0.1,1.1)
ax2n.set_ylim(-0.1,1.1)
ax3n.set_ylim(-0.1,1.1)
ax4n.set_ylim(-0.1,1.1)
ax5n.set_ylim(-0.1,1.1)
ax6n.set_ylim(-0.1,1.1)
ax7n.set_ylim(-0.1,1.1)
ax8n.set_ylim(-0.1,1.1)
ax9n.set_ylim(-0.1,1.1)
ax10n.set_ylim(-0.1,1.1)
ax11n.set_ylim(-0.1,1.1)
ax12n.set_ylim(-0.1,1.1)
ax13n.set_ylim(-0.1,1.1)
ax14n.set_ylim(-0.1,1.1)
ax15n.set_ylim(-0.1,1.1)
ax16n.set_ylim(-0.1,1.1)
ax17n.set_ylim(-0.1,1.1)
ax18n.set_ylim(-0.1,1.1)


ax3n.yaxis.set_label_position("right")
ax6n.yaxis.set_label_position("right")
ax9n.yaxis.set_label_position("right")
ax12n.yaxis.set_label_position("right")
ax15.yaxis.set_label_position("right")
ax18.yaxis.set_label_position("right")


#enable grid lines

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax5.grid(True)
ax6.grid(True)
ax7.grid(True)
ax8.grid(True)
ax9.grid(True)
ax10.grid(True)
ax11.grid(True)
ax12.grid(True)
ax13.grid(True)
ax14.grid(True)
ax15.grid(True)
ax16.grid(True)
ax17.grid(True)
ax18.grid(True)



#set Y labels
ax1.set_ylabel(r"$\mu{g}\/ (h)^{-1}$")
ax4.set_ylabel(r"$\mu{g}\/ (h)^{-1}$")
ax7.set_ylabel(r"$\mu{g}\/ (h)^{-1}$")
ax10.set_ylabel(r"$\mu{g}\/ (h)^{-1}$")
ax13.set_ylabel(r"$\mu{g}\/ (h)^{-1}$")
ax16.set_ylabel(r"$\mu{g}\/(h)^{-1}$")
#set X labels
ax16.set_xlabel(r"$\mu{mol}\/ {L}^{-1}$")
ax17.set_xlabel(r"$\mu{mol}\/ {L}^{-1}$")
ax18.set_xlabel(r"$\mu{mol}\/{L}^{-1}$")

#adjust plot spacing
plt.subplots_adjust(left=0.15, bottom=0.06, right=0.93, top=0.97, wspace=0.05, hspace=0)

#finally draw the plot
plt.show()
