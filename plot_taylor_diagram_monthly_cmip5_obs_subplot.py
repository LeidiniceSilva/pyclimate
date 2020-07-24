# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot taylor diagram from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
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

    def __init__(self, refstd, fig=None, rect=336, label='_', srange=(0., 6.5), extend=False):
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
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1])
        
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
        ax.axis["left"].label.set_text(u'                Standard deviation')
        
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
			"bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)      # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=10, label=label)
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


def import_cmip5_clim(param, area, model):

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5/hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_historical_r1i1p1_197512-200511.nc'.format(path, param, area,	model)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mdl_data


def import_obs_clim(param, area, database):

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_197512-200511.nc'.format(path, param, area, database)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return obs_data	


if __name__=='__main__':
	
	# Import cmip5 model end obs database climatology
	pre_amz_gcm1 = import_cmip5_clim(u'pr', u'amz', u'BCC-CSM1.1')
	tmp_amz_gcm1 = import_cmip5_clim(u'tas', u'amz', u'BCC-CSM1.1')
	pre_neb_gcm1 = import_cmip5_clim(u'pr', u'neb', u'BCC-CSM1.1')
	tmp_neb_gcm1 = import_cmip5_clim(u'tas', u'neb', u'BCC-CSM1.1')
	pre_mato_gcm1 = import_cmip5_clim(u'pr', u'matopiba', u'BCC-CSM1.1')
	tmp_mato_gcm1 = import_cmip5_clim(u'tas', u'matopiba', u'BCC-CSM1.1')

	pre_amz_gcm2 = import_cmip5_clim(u'pr', u'amz', u'BCC-CSM1.1M')
	tmp_amz_gcm2 = import_cmip5_clim(u'tas', u'amz', u'BCC-CSM1.1M')
	pre_neb_gcm2 = import_cmip5_clim(u'pr', u'neb', u'BCC-CSM1.1M')
	tmp_neb_gcm2 = import_cmip5_clim(u'tas', u'neb', u'BCC-CSM1.1M')
	pre_mato_gcm2 = import_cmip5_clim(u'pr', u'matopiba', u'BCC-CSM1.1M')
	tmp_mato_gcm2 = import_cmip5_clim(u'tas', u'matopiba', u'BCC-CSM1.1M')

	pre_amz_gcm3 = import_cmip5_clim(u'pr', u'amz', u'BNU-ESM')
	tmp_amz_gcm3 = import_cmip5_clim(u'tas', u'amz', u'BNU-ESM')
	pre_neb_gcm3 = import_cmip5_clim(u'pr', u'neb', u'BNU-ESM')
	tmp_neb_gcm3 = import_cmip5_clim(u'tas', u'neb', u'BNU-ESM')
	pre_mato_gcm3 = import_cmip5_clim(u'pr', u'matopiba', u'BNU-ESM')
	tmp_mato_gcm3 = import_cmip5_clim(u'tas', u'matopiba', u'BNU-ESM')

	pre_amz_gcm4 = import_cmip5_clim(u'pr', u'amz', u'CanESM2')
	tmp_amz_gcm4 = import_cmip5_clim(u'tas', u'amz', u'CanESM2')
	pre_neb_gcm4 = import_cmip5_clim(u'pr', u'neb', u'CanESM2')
	tmp_neb_gcm4 = import_cmip5_clim(u'tas', u'neb', u'CanESM2')
	pre_mato_gcm4 = import_cmip5_clim(u'pr', u'matopiba', u'CanESM2')
	tmp_mato_gcm4 = import_cmip5_clim(u'tas', u'matopiba', u'CanESM2')

	pre_amz_gcm5 = import_cmip5_clim(u'pr', u'amz', u'CNRM-CM5')
	tmp_amz_gcm5 = import_cmip5_clim(u'tas', u'amz', u'CNRM-CM5')
	pre_neb_gcm5 = import_cmip5_clim(u'pr', u'neb', u'CNRM-CM5')
	tmp_neb_gcm5 = import_cmip5_clim(u'tas', u'neb', u'CNRM-CM5')
	pre_mato_gcm5 = import_cmip5_clim(u'pr', u'matopiba', u'CNRM-CM5')
	tmp_mato_gcm5 = import_cmip5_clim(u'tas', u'matopiba', u'CNRM-CM5')

	pre_amz_gcm6 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-ACCESS-1')
	tmp_amz_gcm6 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-ACCESS-1')
	pre_neb_gcm6 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-ACCESS-1')
	tmp_neb_gcm6 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-ACCESS-1')
	pre_mato_gcm6 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-ACCESS-1')
	tmp_mato_gcm6 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-ACCESS-1')

	pre_amz_gcm7 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-ACCESS-3')
	tmp_amz_gcm7 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-ACCESS-3')
	pre_neb_gcm7 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-ACCESS-3')
	tmp_neb_gcm7 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-ACCESS-3')
	pre_mato_gcm7 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-ACCESS-3')
	tmp_mato_gcm7 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-ACCESS-3')

	pre_amz_gcm8 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-MK36')
	tmp_amz_gcm8 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-MK36')
	pre_neb_gcm8 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-MK36')
	tmp_neb_gcm8 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-MK36')
	pre_mato_gcm8 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-MK36')
	tmp_mato_gcm8 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-MK36')

	pre_amz_gcm9 = import_cmip5_clim(u'pr', u'amz', u'FIO-ESM')
	tmp_amz_gcm9 = import_cmip5_clim(u'tas', u'amz', u'FIO-ESM')
	pre_neb_gcm9 = import_cmip5_clim(u'pr', u'neb', u'FIO-ESM')
	tmp_neb_gcm9 = import_cmip5_clim(u'tas', u'neb', u'FIO-ESM')
	pre_mato_gcm9 = import_cmip5_clim(u'pr', u'matopiba', u'FIO-ESM')
	tmp_mato_gcm9 = import_cmip5_clim(u'tas', u'matopiba', u'FIO-ESM')

	pre_amz_gcm10 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H')
	tmp_amz_gcm10 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H')
	pre_neb_gcm10 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H')
	tmp_neb_gcm10 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H')
	pre_mato_gcm10 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H')
	tmp_mato_gcm10 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H')

	pre_amz_gcm11 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H-CC')
	tmp_amz_gcm11 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H-CC')
	pre_neb_gcm11 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H-CC')
	tmp_neb_gcm11 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H-CC')
	pre_mato_gcm11 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H-CC')
	tmp_mato_gcm11 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H-CC')

	pre_amz_gcm12 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-AO')
	tmp_amz_gcm12 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-AO')
	pre_neb_gcm12 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-AO')
	tmp_neb_gcm12 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-AO')
	pre_mato_gcm12 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-AO')
	tmp_mato_gcm12 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-AO')

	pre_amz_gcm13 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-CC')
	tmp_amz_gcm13 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-CC')
	pre_neb_gcm13 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-CC')
	tmp_neb_gcm13 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-CC')
	pre_mato_gcm13 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-CC')
	tmp_mato_gcm13 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-CC')

	pre_amz_gcm14 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-ES')
	tmp_amz_gcm14 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-ES')
	pre_neb_gcm14 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-ES')
	tmp_neb_gcm14 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-ES')
	pre_mato_gcm14 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-ES')
	tmp_mato_gcm14 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-ES')

	pre_amz_gcm15 = import_cmip5_clim(u'pr', u'amz', u'INMCM4')
	tmp_amz_gcm15 = import_cmip5_clim(u'tas', u'amz', u'INMCM4')
	pre_neb_gcm15 = import_cmip5_clim(u'pr', u'neb', u'INMCM4')
	tmp_neb_gcm15 = import_cmip5_clim(u'tas', u'neb', u'INMCM4')
	pre_mato_gcm15 = import_cmip5_clim(u'pr', u'matopiba', u'INMCM4')
	tmp_mato_gcm15 = import_cmip5_clim(u'tas', u'matopiba', u'INMCM4')

	pre_amz_gcm16 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-LR')
	tmp_amz_gcm16 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-LR')
	pre_neb_gcm16 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-LR')
	tmp_neb_gcm16 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-LR')
	pre_mato_gcm16 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-LR')
	tmp_mato_gcm16 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-LR')

	pre_amz_gcm17 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-MR')
	tmp_amz_gcm17 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-MR')
	pre_neb_gcm17 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-MR')
	tmp_neb_gcm17 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-MR')
	pre_mato_gcm17 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-MR')
	tmp_mato_gcm17 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-MR')

	pre_amz_gcm18 = import_cmip5_clim(u'pr', u'amz', u'LASG-FGOALS-G2')
	tmp_amz_gcm18 = import_cmip5_clim(u'tas', u'amz', u'LASG-FGOALS-G2')
	pre_neb_gcm18 = import_cmip5_clim(u'pr', u'neb', u'LASG-FGOALS-G2')
	tmp_neb_gcm18 = import_cmip5_clim(u'tas', u'neb', u'LASG-FGOALS-G2')
	pre_mato_gcm18 = import_cmip5_clim(u'pr', u'matopiba', u'LASG-FGOALS-G2')
	tmp_mato_gcm18 = import_cmip5_clim(u'tas', u'matopiba', u'LASG-FGOALS-G2')

	pre_amz_gcm19 = import_cmip5_clim(u'pr', u'amz', u'LASG-FGOALS-S2')
	tmp_amz_gcm19 = import_cmip5_clim(u'tas', u'amz', u'LASG-FGOALS-S2')
	pre_neb_gcm19 = import_cmip5_clim(u'pr', u'neb', u'LASG-FGOALS-S2')
	tmp_neb_gcm19 = import_cmip5_clim(u'tas', u'neb', u'LASG-FGOALS-S2')
	pre_mato_gcm19 = import_cmip5_clim(u'pr', u'matopiba', u'LASG-FGOALS-S2')
	tmp_mato_gcm19 = import_cmip5_clim(u'tas', u'matopiba', u'LASG-FGOALS-S2')

	pre_amz_gcm20 = import_cmip5_clim(u'pr', u'amz', u'MIROC5')
	tmp_amz_gcm20 = import_cmip5_clim(u'tas', u'amz', u'MIROC5')
	pre_neb_gcm20 = import_cmip5_clim(u'pr', u'neb', u'MIROC5')
	tmp_neb_gcm20 = import_cmip5_clim(u'tas', u'neb', u'MIROC5')
	pre_mato_gcm20 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC5')
	tmp_mato_gcm20 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC5')

	pre_amz_gcm21 = import_cmip5_clim(u'pr', u'amz', u'MIROC-ESM')
	tmp_amz_gcm21 = import_cmip5_clim(u'tas', u'amz', u'MIROC-ESM')
	pre_neb_gcm21 = import_cmip5_clim(u'pr', u'neb', u'MIROC-ESM')
	tmp_neb_gcm21 = import_cmip5_clim(u'tas', u'neb', u'MIROC-ESM')
	pre_mato_gcm21 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC-ESM')
	tmp_mato_gcm21 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC-ESM')

	pre_amz_gcm22 = import_cmip5_clim(u'pr', u'amz', u'MIROC-ESM-CHEM')
	tmp_amz_gcm22 = import_cmip5_clim(u'tas', u'amz', u'MIROC-ESM-CHEM')
	pre_neb_gcm22 = import_cmip5_clim(u'pr', u'neb', u'MIROC-ESM-CHEM')
	tmp_neb_gcm22 = import_cmip5_clim(u'tas', u'neb', u'MIROC-ESM-CHEM')
	pre_mato_gcm22 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC-ESM-CHEM')
	tmp_mato_gcm22 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC-ESM-CHEM')

	pre_amz_gcm23 = import_cmip5_clim(u'pr', u'amz', u'MPI-ESM-LR')
	tmp_amz_gcm23 = import_cmip5_clim(u'tas', u'amz', u'MPI-ESM-LR')
	pre_neb_gcm23 = import_cmip5_clim(u'pr', u'neb', u'MPI-ESM-LR')
	tmp_neb_gcm23 = import_cmip5_clim(u'tas', u'neb', u'MPI-ESM-LR')
	pre_mato_gcm23 = import_cmip5_clim(u'pr', u'matopiba', u'MPI-ESM-LR')
	tmp_mato_gcm23 = import_cmip5_clim(u'tas', u'matopiba', u'MPI-ESM-LR')

	pre_amz_gcm24 = import_cmip5_clim(u'pr', u'amz', u'MPI-ESM-MR')
	tmp_amz_gcm24 = import_cmip5_clim(u'tas', u'amz', u'MPI-ESM-MR')
	pre_neb_gcm24 = import_cmip5_clim(u'pr', u'neb', u'MPI-ESM-MR')
	tmp_neb_gcm24 = import_cmip5_clim(u'tas', u'neb', u'MPI-ESM-MR')
	pre_mato_gcm24 = import_cmip5_clim(u'pr', u'matopiba', u'MPI-ESM-MR')
	tmp_mato_gcm24 = import_cmip5_clim(u'tas', u'matopiba', u'MPI-ESM-MR')

	pre_amz_gcm25 = import_cmip5_clim(u'pr', u'amz', u'MRI-CGCM3')
	tmp_amz_gcm25 = import_cmip5_clim(u'tas', u'amz', u'MRI-CGCM3')
	pre_neb_gcm25 = import_cmip5_clim(u'pr', u'neb', u'MRI-CGCM3')
	tmp_neb_gcm25 = import_cmip5_clim(u'tas', u'neb', u'MRI-CGCM3')
	pre_mato_gcm25 = import_cmip5_clim(u'pr', u'matopiba', u'MRI-CGCM3')
	tmp_mato_gcm25 = import_cmip5_clim(u'tas', u'matopiba', u'MRI-CGCM3')

	pre_amz_gcm26 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CCSM4')
	tmp_amz_gcm26 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CCSM4')
	pre_neb_gcm26 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CCSM4')
	tmp_neb_gcm26 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CCSM4')
	pre_mato_gcm26 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CCSM4')
	tmp_mato_gcm26 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CCSM4')

	pre_amz_gcm27 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CESM1-BGC')
	tmp_amz_gcm27 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CESM1-BGC')
	pre_neb_gcm27 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CESM1-BGC')
	tmp_neb_gcm27 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CESM1-BGC')
	pre_mato_gcm27 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CESM1-BGC')
	tmp_mato_gcm27 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CESM1-BGC')

	pre_amz_gcm28 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CESM1-CAM5')
	tmp_amz_gcm28 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CESM1-CAM5')
	pre_neb_gcm28 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CESM1-CAM5')
	tmp_neb_gcm28 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CESM1-CAM5')
	pre_mato_gcm28 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CESM1-CAM5')
	tmp_mato_gcm28 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CESM1-CAM5')

	pre_amz_gcm29 = import_cmip5_clim(u'pr', u'amz', u'NorESM1-M')
	tmp_amz_gcm29 = import_cmip5_clim(u'tas', u'amz', u'NorESM1-M')
	pre_neb_gcm29 = import_cmip5_clim(u'pr', u'neb', u'NorESM1-M')
	tmp_neb_gcm29 = import_cmip5_clim(u'tas', u'neb', u'NorESM1-M')
	pre_mato_gcm29 = import_cmip5_clim(u'pr', u'matopiba', u'NorESM1-M')
	tmp_mato_gcm29 = import_cmip5_clim(u'tas', u'matopiba', u'NorESM1-M')

	pre_amz_gcm30 = import_cmip5_clim(u'pr', u'amz', u'NorESM1-ME')
	tmp_amz_gcm30 = import_cmip5_clim(u'tas', u'amz', u'NorESM1-ME')
	pre_neb_gcm30 = import_cmip5_clim(u'pr', u'neb', u'NorESM1-ME')
	tmp_neb_gcm30 = import_cmip5_clim(u'tas', u'neb', u'NorESM1-ME')
	pre_mato_gcm30 = import_cmip5_clim(u'pr', u'matopiba', u'NorESM1-ME')
	tmp_mato_gcm30 = import_cmip5_clim(u'tas', u'matopiba', u'NorESM1-ME')

	pre_amz_gcm31 = import_cmip5_clim(u'pr', u'amz', u'ensmean_cmip5')
	tmp_amz_gcm31 = import_cmip5_clim(u'tas', u'amz', u'ensmean_cmip5')
	pre_neb_gcm31 = import_cmip5_clim(u'pr', u'neb', u'ensmean_cmip5')
	tmp_neb_gcm31 = import_cmip5_clim(u'tas', u'neb', u'ensmean_cmip5')
	pre_mato_gcm31 = import_cmip5_clim(u'pr', u'matopiba', u'ensmean_cmip5')
	tmp_mato_gcm31 = import_cmip5_clim(u'tas', u'matopiba', u'ensmean_cmip5')

	pre_amz_obs  = import_obs_clim(u'pre', u'amz', u'cru_ts4.02')
	tmp_amz_obs  = import_obs_clim(u'tmp', u'amz', u'cru_ts4.02')
	pre_neb_obs  = import_obs_clim(u'pre', u'neb', u'cru_ts4.02')
	tmp_neb_obs  = import_obs_clim(u'tmp', u'neb', u'cru_ts4.02')
	pre_mato_obs  = import_obs_clim(u'pre', u'matopiba', u'cru_ts4.02')
	tmp_mato_obs  = import_obs_clim(u'tmp', u'matopiba', u'cru_ts4.02')
	
	# Reference database standard desviation
	stdrefs = dict(PRE1=1,
				 TMP1=1,
				 PRE2=1,
				 TMP2=1,
				 PRE3=1,
				 TMP3=1)       

	text1 = dict(PRE1='A)',
				 TMP1='D)',
				 PRE2='B)',
				 TMP2='E)',
				 PRE3='C)',
				 TMP3='F')  
				 
	# Compute stddev and correlation coefficient of models
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(PRE1=[[pre_amz_gcm1.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm1)[0,1], 'BCC-CSM1.1'],
                       [pre_amz_gcm2.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [pre_amz_gcm3.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm3)[0,1], 'BNU-ESM'],
                       [pre_amz_gcm4.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm4)[0,1], 'CanESM2'],
                       [pre_amz_gcm5.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm5)[0,1], 'CNRM-CM5'],
                       [pre_amz_gcm6.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [pre_amz_gcm7.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [pre_amz_gcm8.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm8)[0,1], 'CSIRO-MK36'],
                       [pre_amz_gcm9.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm9)[0,1], 'FIO-ESM'],
                       [pre_amz_gcm10.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm10)[0,1], 'GISS-E2-H'],
                       [pre_amz_gcm11.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [pre_amz_gcm12.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm12)[0,1], 'HadGEM2-AO'],
                       [pre_amz_gcm13.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm13)[0,1], 'HadGEM2-CC'],
                       [pre_amz_gcm14.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm14)[0,1], 'HadGEM2-ES'],
                       [pre_amz_gcm15.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm15)[0,1], 'INMCM4'],
                       [pre_amz_gcm16.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [pre_amz_gcm17.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [pre_amz_gcm18.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [pre_amz_gcm19.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [pre_amz_gcm20.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm20)[0,1], 'MIROC5'],
                       [pre_amz_gcm21.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm21)[0,1], 'MIROC-ESM'],
                       [pre_amz_gcm22.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [pre_amz_gcm23.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm23)[0,1], 'MPI-ESM-LR'],
                       [pre_amz_gcm24.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm24)[0,1], 'MPI-ESM-MR'],
                       [pre_amz_gcm25.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm25)[0,1], 'MRI-CGCM3'],
                       [pre_amz_gcm26.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm26)[0,1], 'NCAR-CCSM4'],
                       [pre_amz_gcm27.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [pre_amz_gcm28.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [pre_amz_gcm29.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm29)[0,1], 'NorESM1-M'],
                       [pre_amz_gcm30.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm30)[0,1], 'NorESM1-ME'],
                       [pre_amz_gcm31.std(ddof=1), np.corrcoef(pre_amz_obs, pre_amz_gcm31)[0,1], 'ensmean_cmip5']],     
                 TMP1=[[tmp_amz_gcm1.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm1)[0,1], 'BCC-CSM1.1'],
                       [tmp_amz_gcm2.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [tmp_amz_gcm3.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm3)[0,1], 'BNU-ESM'],
                       [tmp_amz_gcm4.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm4)[0,1], 'CanESM2'],
                       [tmp_amz_gcm5.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm5)[0,1], 'CNRM-CM5'],
                       [tmp_amz_gcm6.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [tmp_amz_gcm7.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [tmp_amz_gcm8.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm8)[0,1], 'CSIRO-MK36'],
                       [tmp_amz_gcm9.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm9)[0,1], 'FIO-ESM'],
                       [tmp_amz_gcm10.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm10)[0,1], 'GISS-E2-H'],
                       [tmp_amz_gcm11.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [tmp_amz_gcm12.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm12)[0,1], 'HadGEM2-AO'],
                       [tmp_amz_gcm13.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm13)[0,1], 'HadGEM2-CC'],
                       [tmp_amz_gcm14.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm14)[0,1], 'HadGEM2-ES'],
                       [tmp_amz_gcm15.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm15)[0,1], 'INMCM4'],
                       [tmp_amz_gcm16.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [tmp_amz_gcm17.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [tmp_amz_gcm18.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [tmp_amz_gcm19.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [tmp_amz_gcm20.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm20)[0,1], 'MIROC5'],
                       [tmp_amz_gcm21.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm21)[0,1], 'MIROC-ESM'],
                       [tmp_amz_gcm22.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [tmp_amz_gcm23.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm23)[0,1], 'MPI-ESM-LR'],
                       [tmp_amz_gcm24.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm24)[0,1], 'MPI-ESM-MR'],
                       [tmp_amz_gcm25.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm25)[0,1], 'MRI-CGCM3'],
                       [tmp_amz_gcm26.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm26)[0,1], 'NCAR-CCSM4'],
                       [tmp_amz_gcm27.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [tmp_amz_gcm28.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [tmp_amz_gcm29.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm29)[0,1], 'NorESM1-M'],
                       [tmp_amz_gcm30.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm30)[0,1], 'NorESM1-ME'],
                       [tmp_amz_gcm31.std(ddof=1), np.corrcoef(pre_amz_obs, tmp_amz_gcm31)[0,1], 'ensmean_cmip5']],
                 PRE2=[[pre_neb_gcm1.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm1)[0,1], 'BCC-CSM1.1'],
                       [pre_neb_gcm2.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [pre_neb_gcm3.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm3)[0,1], 'BNU-ESM'],
                       [pre_neb_gcm4.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm4)[0,1], 'CanESM2'],
                       [pre_neb_gcm5.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm5)[0,1], 'CNRM-CM5'],
                       [pre_neb_gcm6.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [pre_neb_gcm7.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [pre_neb_gcm8.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm8)[0,1], 'CSIRO-MK36'],
                       [pre_neb_gcm9.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm9)[0,1], 'FIO-ESM'],
                       [pre_neb_gcm10.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm10)[0,1], 'GISS-E2-H'],
                       [pre_neb_gcm11.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [pre_neb_gcm12.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm12)[0,1], 'HadGEM2-AO'],
                       [pre_neb_gcm13.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm13)[0,1], 'HadGEM2-CC'],
                       [pre_neb_gcm14.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm14)[0,1], 'HadGEM2-ES'],
                       [pre_neb_gcm15.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm15)[0,1], 'INMCM4'],
                       [pre_neb_gcm16.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [pre_neb_gcm17.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [pre_neb_gcm18.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [pre_neb_gcm19.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [pre_neb_gcm20.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm20)[0,1], 'MIROC5'],
                       [pre_neb_gcm21.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm21)[0,1], 'MIROC-ESM'],
                       [pre_neb_gcm22.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [pre_neb_gcm23.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm23)[0,1], 'MPI-ESM-LR'],
                       [pre_neb_gcm24.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm24)[0,1], 'MPI-ESM-MR'],
                       [pre_neb_gcm25.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm25)[0,1], 'MRI-CGCM3'],
                       [pre_neb_gcm26.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm26)[0,1], 'NCAR-CCSM4'],
                       [pre_neb_gcm27.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [pre_neb_gcm28.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [pre_neb_gcm29.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm29)[0,1], 'NorESM1-M'],
                       [pre_neb_gcm30.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm30)[0,1], 'NorESM1-ME'],
                       [pre_neb_gcm31.std(ddof=1), np.corrcoef(pre_neb_obs, pre_neb_gcm31)[0,1], 'ensmean_cmip5']],
                 TMP2=[[tmp_neb_gcm1.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm1)[0,1], 'BCC-CSM1.1'],
                       [tmp_neb_gcm2.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [tmp_neb_gcm3.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm3)[0,1], 'BNU-ESM'],
                       [tmp_neb_gcm4.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm4)[0,1], 'CanESM2'],
                       [tmp_neb_gcm5.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm5)[0,1], 'CNRM-CM5'],
                       [tmp_neb_gcm6.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [tmp_neb_gcm7.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [tmp_neb_gcm8.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm8)[0,1], 'CSIRO-MK36'],
                       [tmp_neb_gcm9.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm9)[0,1], 'FIO-ESM'],
                       [tmp_neb_gcm10.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm10)[0,1], 'GISS-E2-H'],
                       [tmp_neb_gcm11.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [tmp_neb_gcm12.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm12)[0,1], 'HadGEM2-AO'],
                       [tmp_neb_gcm13.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm13)[0,1], 'HadGEM2-CC'],
                       [tmp_neb_gcm14.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm14)[0,1], 'HadGEM2-ES'],
                       [tmp_neb_gcm15.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm15)[0,1], 'INMCM4'],
                       [tmp_neb_gcm16.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [tmp_neb_gcm17.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [tmp_neb_gcm18.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [tmp_neb_gcm19.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [tmp_neb_gcm20.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm20)[0,1], 'MIROC5'],
                       [tmp_neb_gcm21.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm21)[0,1], 'MIROC-ESM'],
                       [tmp_neb_gcm22.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [tmp_neb_gcm23.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm23)[0,1], 'MPI-ESM-LR'],
                       [tmp_neb_gcm24.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm24)[0,1], 'MPI-ESM-MR'],
                       [tmp_neb_gcm25.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm25)[0,1], 'MRI-CGCM3'],
                       [tmp_neb_gcm26.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm26)[0,1], 'NCAR-CCSM4'],
                       [tmp_neb_gcm27.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [tmp_neb_gcm28.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [tmp_neb_gcm29.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm29)[0,1], 'NorESM1-M'],
                       [tmp_neb_gcm30.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm30)[0,1], 'NorESM1-ME'],
                       [tmp_neb_gcm31.std(ddof=1), np.corrcoef(tmp_neb_obs, tmp_neb_gcm31)[0,1], 'ensmean_cmip5']],
                 PRE3=[[pre_mato_gcm1.std(ddof=1), np.corrcoef(pre_mato_obs, pre_neb_gcm1)[0,1], 'BCC-CSM1.1'],
                       [pre_mato_gcm2.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [pre_mato_gcm3.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm3)[0,1], 'BNU-ESM'],
                       [pre_mato_gcm4.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm4)[0,1], 'CanESM2'],
                       [pre_mato_gcm5.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm5)[0,1], 'CNRM-CM5'],
                       [pre_mato_gcm6.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [pre_mato_gcm7.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [pre_mato_gcm8.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm8)[0,1], 'CSIRO-MK36'],
                       [pre_mato_gcm9.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm9)[0,1], 'FIO-ESM'],
                       [pre_mato_gcm10.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm10)[0,1], 'GISS-E2-H'],
                       [pre_mato_gcm11.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [pre_mato_gcm12.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm12)[0,1], 'HadGEM2-AO'],
                       [pre_mato_gcm13.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm13)[0,1], 'HadGEM2-CC'],
                       [pre_mato_gcm14.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm14)[0,1], 'HadGEM2-ES'],
                       [pre_mato_gcm15.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm15)[0,1], 'INMCM4'],
                       [pre_mato_gcm16.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [pre_mato_gcm17.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [pre_mato_gcm18.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [pre_mato_gcm19.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [pre_mato_gcm20.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm20)[0,1], 'MIROC5'],
                       [pre_mato_gcm21.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm21)[0,1], 'MIROC-ESM'],
                       [pre_mato_gcm22.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [pre_mato_gcm23.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm23)[0,1], 'MPI-ESM-LR'],
                       [pre_mato_gcm24.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm24)[0,1], 'MPI-ESM-MR'],
                       [pre_mato_gcm25.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm25)[0,1], 'MRI-CGCM3'],
                       [pre_mato_gcm26.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm26)[0,1], 'NCAR-CCSM4'],
                       [pre_mato_gcm27.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [pre_mato_gcm28.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [pre_mato_gcm29.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm29)[0,1], 'NorESM1-M'],
                       [pre_mato_gcm30.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm30)[0,1], 'NorESM1-ME'],
                       [pre_mato_gcm31.std(ddof=1), np.corrcoef(pre_mato_obs, pre_mato_gcm31)[0,1], 'ensmean_cmip5']],
                 TMP3=[[tmp_mato_gcm1.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm1)[0,1], 'BCC-CSM1.1'],
                       [tmp_mato_gcm2.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm2)[0,1], 'BCC-CSM1.1M'],
                       [tmp_mato_gcm3.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm3)[0,1], 'BNU-ESM'],
                       [tmp_mato_gcm4.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm4)[0,1], 'CanESM2'],
                       [tmp_mato_gcm5.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm5)[0,1], 'CNRM-CM5'],
                       [tmp_mato_gcm6.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm6)[0,1], 'CSIRO-ACCESS-1'],
                       [tmp_mato_gcm7.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm7)[0,1], 'CSIRO-ACCESS-3'],
                       [tmp_mato_gcm8.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm8)[0,1], 'CSIRO-MK36'],
                       [tmp_mato_gcm9.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm9)[0,1], 'FIO-ESM'],
                       [tmp_mato_gcm10.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm10)[0,1], 'GISS-E2-H'],
                       [tmp_mato_gcm11.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm11)[0,1], 'GISS-E2-H-CC'],
                       [tmp_mato_gcm12.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm12)[0,1], 'HadGEM2-AO'],
                       [tmp_mato_gcm13.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm13)[0,1], 'HadGEM2-CC'],
                       [tmp_mato_gcm14.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm14)[0,1], 'HadGEM2-ES'],
                       [tmp_mato_gcm15.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm15)[0,1], 'INMCM4'],
                       [tmp_mato_gcm16.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm16)[0,1], 'IPSL-CM5A-LR'],
                       [tmp_mato_gcm17.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm17)[0,1], 'IPSL-CM5A-MR'],
                       [tmp_mato_gcm18.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm18)[0,1], 'LASG-FGOALS-G2'],
                       [tmp_mato_gcm19.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm19)[0,1], 'LASG-FGOALS-S2'],
                       [tmp_mato_gcm20.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm20)[0,1], 'MIROC5'],
                       [tmp_mato_gcm21.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm21)[0,1], 'MIROC-ESM'],
                       [tmp_mato_gcm22.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm22)[0,1], 'MIROC-ESM-CHEM'],
                       [tmp_mato_gcm23.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm23)[0,1], 'MPI-ESM-LR'],
                       [tmp_mato_gcm24.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm24)[0,1], 'MPI-ESM-MR'],
                       [tmp_mato_gcm25.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm25)[0,1], 'MRI-CGCM3'],
                       [tmp_mato_gcm26.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm26)[0,1], 'NCAR-CCSM4'],
                       [tmp_mato_gcm27.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm27)[0,1], 'NCAR-CESM1-BGC'],
                       [tmp_mato_gcm28.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm28)[0,1], 'NCAR-CESM1-CAM5'],
                       [tmp_mato_gcm29.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm29)[0,1], 'NorESM1-M'],
                       [tmp_mato_gcm30.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm30)[0,1], 'NorESM1-ME'],
                       [tmp_mato_gcm31.std(ddof=1), np.corrcoef(tmp_mato_obs, tmp_mato_gcm31)[0,1], 'ensmean_cmip5']])
                       
	# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)

	# Here set placement of the points marking 95th and 99th significance
	# levels. For more than 102 samples (degrees freedom > 100), critical
	# correlation levels are 0.195 and 0.254 for 95th and 99th
	# significance levels respectively. Set these by eyeball using the
	# standard deviation x and y axis.

	#~ x95 = [0.01, 0.55] # For Tair, this is for 95th level (r = 0.195)
	#~ y95 = [0.0, 8.45]
	#~ x99 = [0.01, 0.95] # For Tair, this is for 99th level (r = 0.254)
	#~ y99 = [0.0, 8.45]

	x95 = [0.05, 10.9] # For Prcp, this is for 95th level (r = 0.195)
	y95 = [0.0, 80.0]
	x99 = [0.05, 10.0] # For Prcp, this is for 99th level (r = 0.254)
	y99 = [0.0, 80.0]
	
	rects = dict(PRE1=321,
				 TMP1=322,
				 PRE2=323,
				 TMP2=324,
				 PRE3=325,
				 TMP3=326)

	# Plot model end obs data taylor diagram 			 
	fig = plt.figure(figsize=(10, 10))
	
	for var in ['PRE1', 'TMP1', 'PRE2', 'TMP2', 'PRE3', 'TMP3']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label='Reference', srange=(0., 6.5), extend=False)
		dia.samplePoints[0].set_color('r')
		dia.ax.plot(x95,y95,color='blue')
		dia.ax.plot(x99,y99,color='blue')
		
		# Add samples to Taylor diagram
		for i,(stddev,corrcoef,name) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   marker='$%d$' % (i+1), ms=10, ls='',
						   mfc='k', mec='k', # Colors
						   label=name)
			plt.text(1.5, 5, text1[var])

		# Add RMS contours, and label them
		contours = dia.add_contours(colors='0.5')
		plt.clabel(contours, inline=2, fontsize=10, fmt='%.1f')

		# Tricky: ax is the polar ax (used for plots), _ax is the container (used for layout)
		dia.add_grid()                                  
		dia._ax.axis[:].major_ticks.set_tick_out(True) 

	# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
	# Can also use special options here: http://matplotlib.sourceforge.net/users/legend_guide.html
	
	# Add a figure legend
	fig.legend(dia.samplePoints, 
			   [ p.get_label() for p in dia.samplePoints ], 
			   prop=dict(size=10), numpoints=1, loc=(0.75, 0.15))

	plt.subplots_adjust(left=0.10, bottom=0.10, right=0.70, top=0.90, wspace=0.25, hspace=0.20)
    
	# Path out to save figure
	path_out = '/home/nice'
	name_out = 'pyplt_taylor_diagram_cmip5_cru_1975-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()


