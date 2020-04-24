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


def import_sim_season(area, exp, season):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, season, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	season_exp = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
	
	return season_exp


def import_obs_season(area, obs, season):
	
	param = 'precip' # precip, pre or tmp
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, season, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	season_obs = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)

	return season_obs


if __name__=='__main__':

	# Import regcm exps model end obs database seasonaly
	
	nam_exp1_djf = import_sim_season(u'nam', u'regcm_exp1', u'djf')
	sam_exp1_djf = import_sim_season(u'sam', u'regcm_exp1', u'djf')
	neb_exp1_djf = import_sim_season(u'neb', u'regcm_exp1', u'djf')

	nam_exp1_jja = import_sim_season(u'nam', u'regcm_exp1', u'jja')
	sam_exp1_jja = import_sim_season(u'sam', u'regcm_exp1', u'jja')
	neb_exp1_jja = import_sim_season(u'neb', u'regcm_exp1', u'jja')
		
	nam_exp2_djf = import_sim_season(u'nam', u'regcm_exp2', u'djf')
	sam_exp2_djf = import_sim_season(u'sam', u'regcm_exp2', u'djf')
	neb_exp2_djf = import_sim_season(u'neb', u'regcm_exp2', u'djf')

	nam_exp2_jja = import_sim_season(u'nam', u'regcm_exp2', u'jja')
	sam_exp2_jja = import_sim_season(u'sam', u'regcm_exp2', u'jja')
	neb_exp2_jja = import_sim_season(u'neb', u'regcm_exp2', u'jja')
	
	nam_obs_djf = import_obs_season(u'nam', u'gpcp_v2.2_obs', u'djf')
	sam_obs_djf = import_obs_season(u'sam', u'gpcp_v2.2_obs', u'djf')
	neb_obs_djf = import_obs_season(u'neb', u'gpcp_v2.2_obs', u'djf')

	nam_obs_jja = import_obs_season(u'nam', u'gpcp_v2.2_obs', u'jja')
	sam_obs_jja = import_obs_season(u'sam', u'gpcp_v2.2_obs', u'jja')
	neb_obs_jja = import_obs_season(u'neb', u'gpcp_v2.2_obs', u'jja')


	# Reference database standard desviation		   
	stdrefs=sea_obs1.std(ddof=1)          
	
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(DJF=[[nam_exp1_djf.std(ddof=1), st.pearsonr(nam_exp1_djf, nam_obs_djf)[1], 'Exp1_NAMZ_DJF'],
						[sam_exp1_djf.std(ddof=1), st.pearsonr(sam_exp1_djf, sam_obs_djf)[1], 'Exp1_SAMZ_DJF'],
						[neb_exp1_djf.std(ddof=1), st.pearsonr(neb_exp1_djf, neb_obs_djf)[1], 'Exp1_NEB_DJF'],
						[nam_exp2_djf.std(ddof=1), st.pearsonr(nam_exp2_djf, nam_obs_djf)[1], 'Exp2_NAMZ_DJF'],
						[sam_exp2_djf.std(ddof=1), st.pearsonr(sam_exp2_djf, sam_obs_djf)[1], 'Exp2_SAMZ_DJF'],
						[neb_exp2_djf.std(ddof=1), st.pearsonr(neb_exp2_djf, neb_obs_djf)[1], 'Exp2_NEB_DJF']],
					JJA=[[nam_exp1_jja.std(ddof=1), st.pearsonr(nam_exp1_jja, nam_obs_jja)[1], 'Exp1_NAMZ_JJA'],
						[sam_exp1_jja.std(ddof=1), st.pearsonr(sam_exp1_jja, sam_obs_jja)[1], 'Exp1_SAMZ_JJA'],
						[neb_exp1_jja.std(ddof=1), st.pearsonr(neb_exp1_jja, neb_obs_jja)[1], 'Exp1_NEB_JJA'],
						[nam_exp2_jja.std(ddof=1), st.pearsonr(nam_exp2_jja, nam_obs_jja)[1], 'Exp2_NAMZ_JJA'],
						[sam_exp2_jja.std(ddof=1), st.pearsonr(sam_exp2_jja, sam_obs_jja)[1], 'Exp2_SAMZ_JJA'],
						[neb_exp2_jja.std(ddof=1), st.pearsonr(neb_exp2_jja, neb_obs_jja)[1], 'Exp2_NEB_JJA']])	
						   
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

	rects = dict(DJF=121,
				 JJA=122)

	# Plot model end obs data taylor diagram 
	fig = plt.figure()

	fig.suptitle('Diagrama de Taylor de Precipitação (mm/d) 2001-2010')

	for season in ['DJF','JJA']:

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





