# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot taylor diagram from regcm46 and obs database"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import scipy.stats as st
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF

from matplotlib.projections import PolarAxes


class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=336, label='_', marker='', color='', srange=(0., 1.75), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        self.refstd = refstd           
        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1])

        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2

        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0,self.tmax, # 1st quadrant
                                                     self.smin,self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           )

        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text(u'Correlation')

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].toggle(ticklabels=True, label=True)
        ax.axis["left"].label.set_text(u'Standard Deviation')

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
			"bottom" if extend else "left")

        ax.grid(color='gray', axis='x', linestyle='--', linewidth=1)

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)      # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=8, label=label)
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]


    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        'Figure.plot' command.
        """

        l, = self.ax.plot(np.arccos(corrcoef), stddev, *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)

        return l


    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self._ax.grid(*args, **kwargs)


    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax), np.linspace(0, self.tmax))
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours


def import_obs(area, obs, time):

	param = 'pr' # precip, pre or tmp
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	

	data  = netCDF4.Dataset(arq)
	var   = data.variables['cmorph'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return mean_obs


def import_sim(area, exp, time):

	param = 'pr' # pr or tas
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/rcm/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, time, date)	

	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value  = var[:][:,:,:]
	mean_sim = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return mean_sim


if __name__=='__main__':

	# Import regcm exps and obs database 
	djf_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'djf')
	djf_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'djf')
	djf_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'djf')

	jja_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'jja')
	jja_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'jja')
	jja_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'jja')

	ann_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'ann')
	ann_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'ann')
	ann_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'ann')

	djf_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'djf')
	djf_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'djf')
	djf_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'djf')
	
	jja_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'jja')
	jja_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'jja')
	jja_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'jja')

	ann_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'ann')
	ann_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'ann')
	ann_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'ann')

	djf_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'djf')
	djf_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'djf')
	djf_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'djf')

	jja_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'jja')
	jja_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'jja')
	jja_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'jja')
		
	ann_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'ann')
	ann_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'ann')
	ann_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'ann')

	# Reference database standard desviation		   
	stdrefs = dict(DJF=1.0,
				 JJA=1.0,
				 ANN=1.0)  
				 
	text1 = dict(DJF='A)', JJA='B)', ANN='C)')    
	
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(DJF=[[djf_exp1_namz.std(ddof=0), 0.4, 'NAMZ', 'o', 'blue'],
						[djf_exp1_samz.std(ddof=0), 0.1, 'SAMZ', 's', 'blue'],
						[djf_exp1_neb.std(ddof=0),  0.4,  'NEB',  '^', 'blue'],
						[djf_exp2_namz.std(ddof=0), 0.5, 'NAMZ', 'o', 'red'],
						[djf_exp2_samz.std(ddof=0), 0.2, 'SAMZ', 's', 'red'],
						[djf_exp2_neb.std(ddof=0),  0.4,  'NEB',  '^', 'red']],
				   JJA=[[jja_exp1_namz.std(ddof=0), 0.2, 'NAMZ', 'o', 'blue'],
						[jja_exp1_samz.std(ddof=0), 0.7, 'SAMZ', 's', 'blue'],
						[jja_exp1_neb.std(ddof=0),  0.8,  'NEB',  '^', 'blue'],
						[jja_exp2_namz.std(ddof=0), 0.4, 'NAMZ', 'o', 'red'],
						[jja_exp2_samz.std(ddof=0), 0.7, 'SAMZ', 's', 'red'],
						[jja_exp2_neb.std(ddof=0),  0.7,  'NEB',  '^', 'red']],
				   ANN=[[ann_exp1_namz.std(ddof=0), 0.1, 'NAMZ', 'o', 'blue'],
						[ann_exp1_samz.std(ddof=0), 0.1, 'SAMZ', 's', 'blue'],
						[ann_exp1_neb.std(ddof=0),  0.4,  'NEB',  '^', 'blue'],
						[ann_exp2_namz.std(ddof=0), 0.1, 'NAMZ', 'o', 'red'],
						[ann_exp2_samz.std(ddof=0), 0.3, 'SAMZ', 's', 'red'],
						[ann_exp2_neb.std(ddof=0),  0.4,  'NEB',  '^', 'red']])	

	# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)

	# Here set placement of the points marking 95th and 99th significance
	# levels. For more than 102 samples (degrees freedom > 100), critical
	# correlation levels are 0.195 and 0.254 for 95th and 99th
	# significance levels respectively. Set these by eyeball using the
	# standard deviation x and y axis.

	x95 = [0.01, 0.55] 
	y95 = [0.0, 1.76]
	x99 = [0.01, 0.95] 
	y99 = [0.0, 1.76]

	rects = dict(DJF=221,
				 JJA=222,
				 ANN=223)

	# Plot regcm exps and obs database taylor diagram 
	fig = plt.figure(figsize=(7, 7))

	for season in ['DJF','JJA','ANN']:

		dia = TaylorDiagram(stdrefs[season], fig=fig, rect=rects[season], label=u'Reference', srange=(0., 1.75), extend=False)
		dia.samplePoints[0].set_color('r')

		# Add samples to Taylor diagram
		for i,(stddev,corrcoef,name, mark, cor) in enumerate(samples[season]):
			dia.add_sample(stddev, corrcoef, label=name, marker=mark, color='black', mfc=cor, ms=10, ls='')
			plt.text(0.2, 1.2, 'RMSE', fontweight='bold', color='0.6')
			plt.text(-0.3, 1.76, text1[season], fontweight='bold')

		# Add RMS contours, and label them
		contours = dia.add_contours(levels=5, colors='0.5') 
		dia.ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')

	plt.text(2.4, 0.4, 'Reg_Holtslag', color='blue', fontsize=10, fontweight='bold')
	plt.text(2.4, 0.2, 'Reg_UW-PBL', color='red', fontsize=10, fontweight='bold')

	# Add a figure legend
	fig.legend(dia.samplePoints, 
			   [ p.get_label() for p in dia.samplePoints ], 
			   prop=dict(size=10), bbox_to_anchor=(0.575, 0.44), ncol=1, numpoints=1, loc=2)

	# Path out to save figure
	path_out = '/home/nice/Downloads'
	name_out = 'pyplt_taylor_diagram_pr_regcm_pbl_obs_2001-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
