# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/15/2019"
__description__ = "This function compute Taylor Diagram"

import os
import netCDF4
import numpy as np
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
        print "Reference std:", self.refstd
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


def import_cmip5_clim(model):
	
	param = 'pr' # pr or tas
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_amz_neb_Amon_{2}_{3}_{4}.nc'.format(path, param,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mdl_data


def import_obs_clim(database):
	
	param = 'pre' # pre or tmp
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_amz_neb_{2}_obs_mon_{3}.nc'.format(path,
	param, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return obs_data


if __name__=='__main__':

	# Reference database
	database  = u'cru_ts4.02'
	x = import_obs_clim(database)
	
	# Reference database standard desviation
	data = np.cos(x)
	refstd = data.std(ddof=1)
	
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
	
	# Compute stddev and correlation coefficient of models
	samples = np.array([[m.std(ddof=1), np.corrcoef(x, m)[0,1]]
					     for m in (m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33)])
	
	print samples
				      
	fig = plt.figure(figsize=(18,12)) 
	ax1 = fig.add_subplot(1,2,1, xlabel='X', ylabel='Y')

	fig.suptitle(u'Rainfall Taylor Diagram - AMZ_NEB (Lat:85S 15N, Lon:20E 10W) \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Reference period: 1850-2005)', size='x-large')

	# Taylor diagram
	dia = TaylorDiagram(refstd, fig=fig, rect=122, label="Reference")
	colors = plt.matplotlib.cm.jet(np.linspace(0,1,len(samples)))
	ax1.plot(x,data,'ko', label='Data')
	
	for i,m in enumerate([np.cos(m1),np.cos(m2),np.cos(m3),np.cos(m4),np.cos(m5),np.cos(m6),np.cos(m7),np.cos(m8),np.cos(m9),np.cos(m10),np.cos(m11),np.cos(m12),np.cos(m13),np.cos(m14),np.cos(m15),np.cos(m16),np.cos(m17),np.cos(m18),np.cos(m19),np.cos(m20),np.cos(m21),np.cos(m22),np.cos(m23),np.cos(m24),np.cos(m25),np.cos(m26),np.cos(m27),np.cos(m28),np.cos(m29),np.cos(m30),np.cos(m31),np.cos(m32),np.cos(m33)]):
		ax1.plot(x,m, c=colors[i])
	
	# Add samples to Taylor diagram
	for i,(stddev,corrcoef) in enumerate(samples):
		
		mld_dic = {0:'BCC-CSM1.1',1:'BCC-CSM1.1M',2:'BNU-ESM',3:'CanESM2',4:'CNRM-CM5',5:'CSIRO-ACCESS-1',6:'CSIRO-ACCESS-3',
		7:'CSIRO-MK36',8:'FIO-ESM',9:'GISS-E2-H-CC',10:'GISS-E2-H',11:'GISS-E2-R',12:'HadGEM2-AO',13:'HadGEM2-CC',14:'HadGEM2-ES',
		15:'INMCM4',16:'IPSL-CM5A-LR',17:'IPSL-CM5A-MR',18:'IPSL-CM5B-LR',19:'LASG-FGOALS-G2',20:'LASG-FGOALS-S2',21:'MIROC',
		22:'MIROC-ESM-CHEM',23:'MIROC-ESM',24:'MPI-ESM-LR',25:'MPI-ESM-MR',26:'MRI-CGCM3',27:'NCAR-CCSM4',28:'NCAR-CESM1-BGC',
		29:'NCAR-CESM1-CAM5',30:'NorESM1-ME',31:'NorESM1-M',32:'ENSMEAN_CMIP5'}
		
		dia.add_sample(stddev, corrcoef, marker='s', ls='', c=colors[i], label=mld_dic[i])
	
	# Add RMS contours, and label them
	contours = dia.add_contours(colors='0.5')
	plt.clabel(contours, inline=1, fontsize=10)
	
	# Add a figure legend
	fig.legend(dia.samplePoints, [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), loc='upper right')
	
	path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
	name_out = 'pyplt_taylor_diagram_pre_amz_neb_cmip5_cru_1975-2005.png'

	if not os.path.exists(path_out):
		create_path(path_out)
		
	plt.savefig(os.path.join(path_out, name_out))
	plt.show()
	exit()


