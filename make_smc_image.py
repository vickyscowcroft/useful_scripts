#!/usr/bin/env python

import matplotlib.pyplot as mp
#matplotlib.use('Agg')
import sys
import os
sys.path.append('/usr/Montage_v3.3')

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

import aplpy
import numpy
import montage_wrapper as montage



usage = '''Usage: make_galaxy_image.py input_image output_image
			    Plots the galaxy image with a coordinate grid and outputs to eps file'''

#if len(sys.argv) < 3:
#	print usage
#	sys.exit(1)
#output = sys.argv[2]

input = 'smc_dss.fits'
inset = 'smc_dss_huge.fits'
output = 'smc_with_inset.pdf'


fig = mp.figure(figsize=(10,10))

mosaic = aplpy.FITSFigure(input,north='True', figure = fig)
mosaic.show_grayscale(vmin=12000,vmax=16500,invert='true')
mosaic.tick_labels.set_font(size='small')
mosaic.tick_labels.set_xformat("hh:mm:ss")
mosaic.set_theme('publication')

#ic1613_coverage = aplpy.FITSFigure('/Users/vs/IC1613/mosaic_cov.fits')

#ic1613.show_contour('/Users/vs/IC1613/mosaic_cov.fits', colors='black')

#mosaic.show_grayscale(vmin=10,vmax=185,invert='true')

data = numpy.loadtxt('/Users/vs/Dropbox/SMC/smc_cepheids_degrees')
ra, dec = data[:,0], data[:,1]
print ra, dec

mosaic.show_markers(ra,dec, edgecolor='mediumblue',facecolor='darkorange', marker='o',s=30,alpha=1.0)
#mosaic.show_markers(ra,dec, edgecolor='mediumblue',facecolor='none', marker='o',s=30,alpha=1.0)

f1 = aplpy.FITSFigure(inset,north='True', figure = fig, subplot=[0.63,0.13,0.25,0.125])
f1.show_grayscale(vmin=12000,vmax=16500,invert='true')
f1.set_theme('publication')

f1.show_markers(ra,dec, edgecolor='mediumblue',facecolor='darkorange', marker='o',s=8,alpha=1.0)
#f1.show_markers(ra,dec, edgecolor='mediumblue',facecolor='none', marker='o',s=8,alpha=1.0)


f1.hide_yaxis_label()
f1.hide_ytick_labels()
f1.hide_xaxis_label()
f1.hide_xtick_labels()


#mosaic.show_rectangles(16.20727, 2.1179277, 0.201738, 0.200354)



mosaic.save(output)

