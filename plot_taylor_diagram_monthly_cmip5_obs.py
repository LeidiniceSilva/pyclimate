# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot taylor diagram from CMIP5 models end OBS basedata"


import os
import netCDF4
import numpy as NP
import matplotlib.pyplot as PLT
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

    def __init__(self, refstd, fig=None, rect=111, label='_', srange=(0, 5), extend=False):
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
        rlocs = NP.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1])
        
        if extend:
            # Diagram extended to negative correlations
            self.tmax = NP.pi
            rlocs = NP.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = NP.pi/2
            
        tlocs = NP.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd


        ghelper = FA.GridHelperCurveLinear(
            tr, extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)

        if fig is None:
            fig = PLT.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text(u'Correlação')

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text(u'Desvio Padrão')
        
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)      # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=10, label=label)
        t = NP.linspace(0, self.tmax)
        r = NP.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        l, = self.ax.plot(NP.arccos(corrcoef), stddev, *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)
        return l
        

    def add_grid(self, *args, **kwargs):
        """Add a grid."""
        self._ax.grid(*args, **kwargs)
        

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """
        
        rs, ts = NP.meshgrid(NP.linspace(self.smin, self.smax), NP.linspace(0, self.tmax))
        
        # Compute centered RMS difference
        rms = NP.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*NP.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        return contours


def import_cmip5_clim(model):
	
	param = 'pr' # pr or tas
	area  = 'amz' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = NP.nanmean(NP.nanmean(value, axis=1), axis=1)

	return mdl_data


def import_obs_clim(database):
	
	param = 'pre' # pre or tmp
	area  = 'amz' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, 
	database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = NP.nanmean(NP.nanmean(value, axis=1), axis=1)

	return obs_data


if __name__=='__main__':

	# Reference database
	database  = u'cru_ts4.02'
	x = import_obs_clim(database)
	
	# Models
	model  = u'BCC-CSM1.1'
	m1 = import_cmip5_clim(model)
	
	model  = u'BCC-CSM1.1M'
	m2 = import_cmip5_clim(model)
	
	model  = u'BNU-ESM'
	m3 = import_cmip5_clim(model)
	
	model  = u'CanESM2'
	m4 = import_cmip5_clim(model)
	
	model  = u'CNRM-CM5'
	m5 = import_cmip5_clim(model)
	
	model  = u'CSIRO-ACCESS-1'
	m6 = import_cmip5_clim(model)
	
	model  = u'CSIRO-ACCESS-3'
	m7 = import_cmip5_clim(model)
	
	model  = u'CSIRO-MK36'
	m8 = import_cmip5_clim(model)
	
	model  = u'FIO-ESM'
	m9 = import_cmip5_clim(model)
	
	model  = u'GISS-E2-H-CC'
	m10 = import_cmip5_clim(model)
	
	model  = u'GISS-E2-H'
	m11 = import_cmip5_clim(model)
	
	model  = u'GISS-E2-R'
	m12 = import_cmip5_clim(model)

	model  = u'HadGEM2-AO'
	m13 = import_cmip5_clim(model)

	model  = u'HadGEM2-CC'
	m14 = import_cmip5_clim(model)
	
	model  = u'HadGEM2-ES'
	m15 = import_cmip5_clim(model)

	model  = u'INMCM4'
	m16 = import_cmip5_clim(model)

	model  = u'IPSL-CM5A-LR'
	m17 = import_cmip5_clim(model)

	model  = u'IPSL-CM5A-MR'
	m18 = import_cmip5_clim(model)

	model  = u'IPSL-CM5B-LR'
	m19 = import_cmip5_clim(model)

	model  = u'LASG-FGOALS-G2'
	m20 = import_cmip5_clim(model)

	model  = u'LASG-FGOALS-S2'
	m21 = import_cmip5_clim(model)
	
	model  = u'MIROC5'
	m22 = import_cmip5_clim(model)

	model  = u'MIROC-ESM-CHEM'
	m23 = import_cmip5_clim(model)

	model  = u'MIROC-ESM'
	m24 = import_cmip5_clim(model)

	model  = u'MPI-ESM-LR'
	m25 = import_cmip5_clim(model)

	model  = u'MPI-ESM-MR'
	m26 = import_cmip5_clim(model)

	model  = u'MRI-CGCM3'
	m27 = import_cmip5_clim(model)

	model  = u'NCAR-CCSM4'
	m28 = import_cmip5_clim(model)
	
	model  = u'NCAR-CESM1-BGC'
	m29 = import_cmip5_clim(model)

	model  = u'NCAR-CESM1-CAM5'
	m30 = import_cmip5_clim(model)

	model  = u'NorESM1-ME'
	m31 = import_cmip5_clim(model)

	model  = u'NorESM1-M'
	m32 = import_cmip5_clim(model)

	model  = u'ensmean_cmip5'
	m33 = import_cmip5_clim(model)
	
		
def plot_diagram():

	# Compute stddev and correlation coefficient of models
	samples = NP.array([[m.std(ddof=1), NP.corrcoef(x, m)[0,1]]
					     for m in (m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33)])
	
	stdref = 1.
	fig = PLT.figure(figsize=(11.5,5))
	dia = TaylorDiagram(stdref, fig=fig, label='Reference', extend=True)
	dia.samplePoints[0].set_color('r')
	
	for i, (stddev, corrcoef) in enumerate(samples):
		mld_dic = {0:'BCC-CSM1.1',1:'BCC-CSM1.1M',2:'BNU-ESM',3:'CanESM2',4:'CNRM-CM5',5:'CSIRO-ACCESS-1',6:'CSIRO-ACCESS-3',
		7:'CSIRO-MK36',8:'FIO-ESM',9:'GISS-E2-H-CC',10:'GISS-E2-H',11:'GISS-E2-R',12:'HadGEM2-AO',13:'HadGEM2-CC',14:'HadGEM2-ES',
		15:'INMCM4',16:'IPSL-CM5A-LR',17:'IPSL-CM5A-MR',18:'IPSL-CM5B-LR',19:'LASG-FGOALS-G2',20:'LASG-FGOALS-S2',21:'MIROC',
		22:'MIROC-ESM-CHEM',23:'MIROC-ESM',24:'MPI-ESM-LR',25:'MPI-ESM-MR',26:'MRI-CGCM3',27:'NCAR-CCSM4',28:'NCAR-CESM1-BGC',
		29:'NCAR-CESM1-CAM5',30:'NorESM1-ME',31:'NorESM1-M',32:'ENSMEAN_CMIP5'}
		colors = PLT.matplotlib.cm.Set1(NP.linspace(0,1,len(samples)))
		
		dia.add_sample(stddev, corrcoef, marker='$%d$' % (i+1), ms=13, ls='', mfc=colors[i], mec=colors[i], label=mld_dic[i])
	
	# Add RMS contours, and label them
	contours = dia.add_contours(levels=5, colors='0.4')
	PLT.clabel(contours, inline=1, fontsize=13, fmt='%.0f')
	
	dia.add_grid()
	dia._ax.axis[:].major_ticks.set_tick_out(True)
	
	out_var    = u'pre' # pre or tmp
	out_area   = u'amz' # amz or neb
	area_name  = u'AMZ (Lat:16S 4N, Lon:74W 48W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)
	
	if out_var == 'pre':
		var_name   = u'Precipitação'
	else:
		var_name   = u'Temperatura' 
	
	# Add a figure legend and title
	fig.legend(dia.samplePoints, [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size=8), loc='right')		
	PLT.title(u'Diagrama de Taylor de {0} - {1}  \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Período de Referência: 1850-2005)'.format(var_name, area_name), fontsize=15, y=1.1)
		
	return dia, out_var, out_area


if __name__ == '__main__':

    dia, out_var, out_area = plot_diagram()
    
    path_out = '/home/nice'
    name_out = 'pyplt_taylor_diagram_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)
    if not os.path.exists(path_out): create_path(path_out)
    
    PLT.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    PLT.show()
    exit()


