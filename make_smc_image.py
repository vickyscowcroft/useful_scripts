#!/usr/bin/env python

import matplotlib.pyplot as mp
import sys
import os
sys.path.append('/usr/Montage_v3.3')

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

import aplpy
import numpy
import montage_wrapper as montage

### aplpy needs montage and montage_wrapper installed to produce finder charts
### montage is not a python package - it's a separate program --- see aplpy documentation
### montage_wrapper is the python package that interacts with it.

#usage = '''Usage: make_galaxy_image.py input_image output_image
#			    Plots the galaxy image with a coordinate grid and outputs to eps file'''

#if len(sys.argv) < 3:
#	print usage
#	sys.exit(1)
#output = sys.argv[2]

### Have these hardcoded right now, but you can change them to be whatever

input = 'finder_chart_example/smc_dss.fits'
inset = 'finder_chart_example/smc_dss_huge.fits'
output = 'finder_chart_example/smc_with_inset.pdf'


fig = mp.figure(figsize=(10,10))

mosaic = aplpy.FITSFigure(input,north='True', figure = fig)
mosaic.show_grayscale(vmin=12000,vmax=16500,invert='true')
mosaic.tick_labels.set_font(size='small')
mosaic.tick_labels.set_xformat("hh:mm:ss")
mosaic.set_theme('publication')

### Easy to change this using astropy wcs packages to read in whatever format (e.g. hours etc.)
## Should update this
data = numpy.loadtxt('finder_chart_example/smc_cepheids_degrees')
ra, dec = data[:,0], data[:,1]
print ra, dec

mosaic.show_markers(ra,dec, edgecolor='mediumblue',facecolor='darkorange', marker='o',s=30,alpha=1.0)

f1 = aplpy.FITSFigure(inset,north='True', figure = fig, subplot=[0.63,0.13,0.25,0.125])
f1.show_grayscale(vmin=12000,vmax=16500,invert='true')
f1.set_theme('publication')

f1.show_markers(ra,dec, edgecolor='mediumblue',facecolor='darkorange', marker='o',s=8,alpha=1.0)


f1.hide_yaxis_label()
f1.hide_ytick_labels()
f1.hide_xaxis_label()
f1.hide_xtick_labels()

mosaic.save(output)

