# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/15/2019"
__description__ = "This script plot taylor diagram graphics from Rec_EXP models end OBS basedata"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import scipy.stats as st
import matplotlib.pyplot as plt


class TaylorDiagram(object):
	
    """Taylor diagram: plot model standard deviation and correlation
    to reference (obs database) sample in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=111, label='_'):
        """Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using mpl_toolkits.axisartist.floating_axes. refstd is
        the reference standard deviation to be compared to.
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

		# Reference standard deviation
        self.refstd = refstd            

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.concatenate((np.arange(10)/10.,[0.95,0.99]))
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str,rlocs))))

        # Standard deviation axis extent
        self.smin = 0
        self.smax = 3.5*self.refstd

        ghelper = FA.GridHelperCurveLinear(tr, extremes=(0,np.pi/2, # 1st quadrant
                                           self.smin,self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           )

        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom") # "X axis"
        ax.axis["left"].label.set_text("Standard Deviation")

        ax.axis["right"].set_axis_direction("top")   # "Y axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")

        ax.axis["bottom"].set_visible(False)         # Useless

        # Contours along standard deviations
        ax.grid(True)

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        print("Reference std:", self.refstd)
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=10, label=label)
        t = np.linspace(0, np.pi/2)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t,r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]


    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        
        """Add sample (stddev,corrcoeff) to the Taylor diagram. args
        and kwargs are directly propagated to the Figure.plot
        command."""

        l, = self.ax.plot(np.arccos(corrcoef), stddev,  *args, **kwargs) # (theta,radius)
        self.samplePoints.append(l)

        return l


    def add_contours(self, levels=5, **kwargs):
        
        """Add constant centered RMS difference contours."""

        rs,ts = np.meshgrid(np.linspace(self.smin,self.smax),
                            np.linspace(0,np.pi/2))
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours


def import_sim_season(exp):
	
	param = 'pr' # pr or tas
	area  = 'amz_neb' # amz or neb
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	sea_exp = np.nanmean(np.nanmean(value[2:120:3,:,:], axis=1), axis=1)
	sea1_exp = sea_exp[3:40:4]
	sea2_exp = sea_exp[0:40:4]
	sea3_exp = sea_exp[1:40:4]
	sea4_exp = sea_exp[2:40:4]
	
	return sea1_exp, sea2_exp, sea3_exp, sea4_exp


def import_obs_season(database):
	
	param = 'precip' # precip, pre or tmp
	area  = 'amz_neb' # amz or neb
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	sea_obs = np.nanmean(np.nanmean(value[2:120:3,:,:], axis=1), axis=1)
	sea1_obs = sea_obs[3:40:4]
	sea2_obs = sea_obs[0:40:4]
	sea3_obs = sea_obs[1:40:4]
	sea4_obs = sea_obs[2:40:4]

	return sea_obs, sea1_obs, sea2_obs, sea3_obs, sea4_obs


if __name__=='__main__':

	# Import regcm exps model end obs database seasonaly
	exp  = u'regcm_exp1'
	sea1_exp1, sea2_exp1, sea3_exp1, sea4_exp1 = import_sim_season(exp)
			
	exp  = u'regcm_exp2'
	sea1_exp2, sea2_exp2, sea3_exp2, sea4_exp2 = import_sim_season(exp)

	obs  = u'gpcp_v2.2_obs'
	sea_obs1, sea1_obs1, sea2_obs1, sea3_obs1, sea4_obs1 = import_obs_season(obs)

	# Reference database standard desviation		   
	stdrefs=sea_obs1.std(ddof=1)          
	
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(DJF=[[sea1_exp1.std(ddof=1), st.pearsonr(sea1_exp1, sea1_obs1)[1], 'Exp1_NAMZ'],
						[sea1_exp1.std(ddof=1), st.pearsonr(sea1_exp1, sea1_obs1)[1], 'Exp1_SAMZ'],
						[sea1_exp1.std(ddof=1), st.pearsonr(sea1_exp1, sea1_obs1)[1], 'Exp1_NEB'],
						[sea1_exp2.std(ddof=1), st.pearsonr(sea1_exp2, sea1_obs1)[1], 'Exp2_NAMZ'],
						[sea1_exp2.std(ddof=1), st.pearsonr(sea1_exp2, sea1_obs1)[1], 'Exp2_SAMZ'],
						[sea1_exp2.std(ddof=1), st.pearsonr(sea1_exp2, sea1_obs1)[1], 'Exp2_NEB']],              
				   MAM=[[sea2_exp1.std(ddof=1), st.pearsonr(sea2_exp1, sea2_obs1)[1], 'Exp1_NAMZ'],
						[sea2_exp1.std(ddof=1), st.pearsonr(sea2_exp1, sea2_obs1)[1], 'Exp1_SAMZ'],
						[sea2_exp1.std(ddof=1), st.pearsonr(sea2_exp1, sea2_obs1)[1], 'Exp1_NEB'],
						[sea2_exp2.std(ddof=1), st.pearsonr(sea2_exp2, sea2_obs1)[1], 'Exp2_NAMZ'],
						[sea2_exp2.std(ddof=1), st.pearsonr(sea2_exp2, sea2_obs1)[1], 'Exp2_SAMZ'],
						[sea2_exp2.std(ddof=1), st.pearsonr(sea2_exp2, sea2_obs1)[1], 'Exp2_NEB']],
				   JJA=[[sea3_exp1.std(ddof=1), st.pearsonr(sea3_exp1, sea3_obs1)[1], 'Exp1_NAMZ'],
						[sea3_exp1.std(ddof=1), st.pearsonr(sea3_exp1, sea3_obs1)[1], 'Exp1_SAMZ'],
						[sea3_exp1.std(ddof=1), st.pearsonr(sea3_exp1, sea3_obs1)[1], 'Exp1_NEB'],
						[sea3_exp2.std(ddof=1), st.pearsonr(sea3_exp2, sea3_obs1)[1], 'Exp2_NAMZ'],
						[sea3_exp2.std(ddof=1), st.pearsonr(sea3_exp2, sea3_obs1)[1], 'Exp2_SAMZ'],
						[sea3_exp2.std(ddof=1), st.pearsonr(sea3_exp2, sea3_obs1)[1], 'Exp2_NEB']],
				   SON=[[sea4_exp1.std(ddof=1), st.pearsonr(sea4_exp1, sea4_obs1)[1], 'Exp1_NAMZ'],
						[sea4_exp1.std(ddof=1), st.pearsonr(sea4_exp1, sea4_obs1)[1], 'Exp1_SAMZ'],
						[sea4_exp1.std(ddof=1), st.pearsonr(sea4_exp1, sea4_obs1)[1], 'Exp1_NEB'],
						[sea4_exp2.std(ddof=1), st.pearsonr(sea4_exp2, sea4_obs1)[1], 'Exp2_NAMZ'],
						[sea4_exp2.std(ddof=1), st.pearsonr(sea4_exp2, sea4_obs1)[1], 'Exp2_SAMZ'],
						[sea4_exp2.std(ddof=1), st.pearsonr(sea4_exp2, sea4_obs1)[1], 'Exp2_NEB']])	
						   
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

	x95 = [0.05, 10.9] # For Prcp, this is for 95th level (r = 0.195)
	y95 = [0.0, 80.0]
	x99 = [0.05, 10.0] # For Prcp, this is for 99th level (r = 0.254)
	y99 = [0.0, 80.0]

	rects = dict(DJF=221,
				 MAM=222,
				 JJA=223,
				 SON=224)

	# Plot model end obs data taylor diagram 
	fig = plt.figure()

	fig.suptitle('Diagrama de Taylor de Precipitação (mm) 2001-2010 \n NAMZ () SAMZ () NEB()')

	for season in ['DJF','MAM','JJA','SON']:

		dia = TaylorDiagram(stdrefs, fig=fig, rect=rects[season], label='Reference')

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
		dia.ax.clabel(contours, inline=1, fontsize=6, fmt='%.1f')
		# Tricky: ax is the polar ax (used for plots), _ax is the
		# container (used for layout)
		dia._ax.set_title(season)

	# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
	# Can also use special options here:
	# http://matplotlib.sourceforge.net/users/legend_guide.html

	fig.legend(dia.samplePoints,
			   [ p.get_label() for p in dia.samplePoints ],
			   numpoints=1, prop=dict(size=6), ncol=4, loc='lower center')

	# Path out to save figure
	path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
	name_out = 'pyplt_taylor_diagram_pr_regcm_pbl_obs_2001-2010.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

	plt.show()
	exit()





