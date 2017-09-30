#!/usr/bin/env python
"""-----------------------------------------------------------------
  Python file for plotting Finesse ouput VaryMMTT_95MM_wSqz_MMTTon_maxtem4.out
  created automatically Tue Sep 13 22:23:36 2016

  Run from command line as: python VaryMMTT_95MM_wSqz_MMTTon_maxtem4.py
  Load from python script as: import VaryMMTT_95MM_wSqz_MMTTon_maxtem4
  And then use:
  VaryMMTT_95MM_wSqz_MMTTon_maxtem4.run() for plotting only
  x,y=VaryMMTT_95MM_wSqz_MMTTon_maxtem4.run() for plotting and loading the data
  x,y=VaryMMTT_95MM_wSqz_MMTTon_maxtem4.run(1) for only loading the data
-----------------------------------------------------------------"""

__author__ = "Finesse, http://www.gwoptics.org/finesse"

import numpy as np
import matplotlib
BACKEND = 'Qt4Agg'
matplotlib.use(BACKEND)
from matplotlib import rc
import matplotlib.pyplot as plt
formatter = matplotlib.ticker.EngFormatter(unit='')
formatter.ENG_PREFIXES[-6] = 'u'
def run(noplot=None):
	data = np.loadtxt('VaryMMTT_95MM_wSqz_MMTTon_maxtem4.out',comments='%')
	rows,cols=data.shape
	x=data[:,0]
	y=data[:,1:cols]
	mytitle='VaryMMTT_95MM_wSqz_MMTTon_maxtem4                Tue Sep 13 22:23:36 2016'
	if (noplot==None):
		# setting default font sizes
		rc('font',**pp.font)
		rc('xtick',labelsize=pp.TICK_SIZE)
		rc('ytick',labelsize=pp.TICK_SIZE)
		rc('text', usetex=pp.USETEX)
		rc('axes', labelsize = pp.LABEL_SIZE)
		fig=plt.figure()
		fig.set_size_inches(pp.fig_size)
		fig.set_dpi(pp.FIG_DPI)
		from itertools import cycle
		clist = matplotlib.rcParams['axes.color_cycle']
		colorcycler= cycle(clist)
		ax1 = fig.add_subplot(111)
		ax1.set_xlim(5,5000)
		ax1.set_xscale('log',nonposz='clip')
		ax1.set_xlabel('f [Hz] (darm)')
		ax1.set_ylabel('Re ')
		ax1.yaxis.set_major_formatter(formatter)
		ax2 = ax1.twinx()
		trace1=ax1.plot(x, y[:,0], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'NSR_with_RP nOMC_AROC_trans : Re ')
		trace2=ax2.plot(x, y[:,1], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'NSR_with_RP nOMC_AROC_trans : Im ')
		trace3=ax1.plot(x, y[:,2], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM00 nSRMHRa : Re ')
		trace4=ax2.plot(x, y[:,3], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM00 nSRMHRa : Im ')
		trace5=ax1.plot(x, y[:,4], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM01 nSRMHRa : Re ')
		trace6=ax2.plot(x, y[:,5], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM01 nSRMHRa : Im ')
		trace7=ax1.plot(x, y[:,6], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM02 nSRMHRa : Re ')
		trace8=ax2.plot(x, y[:,7], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'nSRMHRaTEM02 nSRMHRa : Im ')
		trace9=ax1.plot(x, y[:,8], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRCoutx nIBAin : Re ')
		trace10=ax2.plot(x, y[:,9], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRCoutx nIBAin : Im ')
		trace11=ax1.plot(x, y[:,10], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRCouty nIBAin : Re ')
		trace12=ax2.plot(x, y[:,11], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRCouty nIBAin : Im ')
		trace13=ax1.plot(x, y[:,12], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRMYqx nSRMHRa : Re ')
		trace14=ax2.plot(x, y[:,13], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRMYqx nSRMHRa : Im ')
		trace15=ax1.plot(x, y[:,14], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRMYqy nSRMHRa : Re ')
		trace16=ax2.plot(x, y[:,15], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'SRMYqy nSRMHRa : Im ')
		trace17=ax1.plot(x, y[:,16], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMXqx nITMX2 : Re ')
		trace18=ax2.plot(x, y[:,17], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMXqx nITMX2 : Im ')
		trace19=ax1.plot(x, y[:,18], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMXqy nITMX2 : Re ')
		trace20=ax2.plot(x, y[:,19], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMXqy nITMX2 : Im ')
		trace21=ax1.plot(x, y[:,20], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMYqx nITMY2 : Re ')
		trace22=ax2.plot(x, y[:,21], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMYqx nITMY2 : Im ')
		trace23=ax1.plot(x, y[:,22], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMYqy nITMY2 : Re ')
		trace24=ax2.plot(x, y[:,23], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'ITMYqy nITMY2 : Im ')
		trace25=ax1.plot(x, y[:,24], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OMCqx nOMC_HROC_refl : Re ')
		trace26=ax2.plot(x, y[:,25], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OMCqx nOMC_HROC_refl : Im ')
		trace27=ax1.plot(x, y[:,26], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OMCqy nOMC_HROC_refl : Re ')
		trace28=ax2.plot(x, y[:,27], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OMCqy nOMC_HROC_refl : Im ')
		trace29=ax1.plot(x, y[:,28], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OFIqx nIBAin : Re ')
		trace30=ax2.plot(x, y[:,29], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OFIqx nIBAin : Im ')
		trace31=ax1.plot(x, y[:,30], '-', linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OFIqy nIBAin : Re ')
		trace32=ax2.plot(x, y[:,31], '-', dashes=[8, 4, 2, 4, 2, 4], linewidth=pp.LINEWIDTH, color = next(colorcycler), label = 'OFIqy nIBAin : Im ')
		ax2.set_ylabel('Im ')
		ax2.yaxis.set_major_formatter(formatter)
		traces = trace1 + trace2 + trace3 + trace4 + trace5 + trace6 + trace7 + trace8 + trace9 + trace10 + trace11 + trace12 + trace13 + trace14 + trace15 + trace16 + trace17 + trace18 + trace19 + trace20 + trace21 + trace22 + trace23 + trace24 + trace25 + trace26 + trace27 + trace28 + trace29 + trace30 + trace31 + trace32
		traces_a = trace1 + trace3 + trace5 + trace7 + trace9 + trace11 + trace13 + trace15 + trace17 + trace19 + trace21 + trace23 + trace25 + trace27 + trace29 + trace31
		traces_p = trace2 + trace4 + trace6 + trace8 + trace10 + trace12 + trace14 + trace16 + trace18 + trace20 + trace22 + trace24 + trace26 + trace28 + trace30 + trace32
		legends = [t.get_label() for t in traces]
		ax2.legend(traces, legends, loc=0, shadow=pp.SHADOW,prop={'size':pp.LEGEND_SIZE})
		ax1.grid(pp.GRID)
		if pp.PRINT_TITLE:
			plt.title(mytitle)
		if pp.SCREEN_TITLE:
			fig.canvas.manager.set_window_title(mytitle)
		else:
			fig.canvas.manager.set_window_title('')
	return (x,y)
class pp():
	# set some gobal settings first
	BACKEND = 'Qt4Agg' # matplotlib backend
	FIG_DPI=90 # DPI of on sceen plot
	# Some help in calculating good figure size for Latex
	# documents. Starting with plot size in pt,
	# get this from LaTeX using \showthe\columnwidth
	fig_width_pt = 484.0
	inches_per_pt = 1.0/72.27  # Convert TeX pt to inches
	golden_mean = (np.sqrt(5)-1.0)/2.0   # Aesthetic ratio
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*golden_mean      # height in inches
	fig_size = [fig_width,fig_height]
	# some plot options:
	LINEWIDTH = 1 # linewidths of traces in plot
	AA = True # antialiasing of traces
	USETEX = False # use Latex encoding in text
	SHADOW = False # shadow of legend box
	GRID = True # grid on or off
	# font sizes for normal text, tick labels and legend
	FONT_SIZE = 10 # size of normal text
	TICK_SIZE = 10 # size of tick labels
	LABEL_SIZE = 10 # size of axes labels
	LEGEND_SIZE = 10 # size of legend
	# font family and type
	font = {'family':'sans-serif','size':FONT_SIZE}
	DPI=300 # DPI for saving via savefig
	# print options given to savefig command:
	print_options = {'dpi':DPI, 'transparent':True, 'bbox_inches':'tight', 'pad_inches':0.1}
	# for Palatino and other serif fonts use:
	#font = {'family':'serif','serif':['Palatino']}
	SCREEN_TITLE = True # show title on screen?
	PRINT_TITLE = False # show title in saved file?

if __name__=="__main__":
	run()
