#!/usr/bin/env/python

import numpy as np
import matplotlib.pyplot as mp
import sys
import gloess_fits as gf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import os
os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']



du = []
db = []
dv = []
dr = []
di = []
dj = []
dh = []
dk = []
dir1 = []
dir2 = []
dir3 = []
dir4 = []
deu = []
deb = []
dev = []
der = []
dei = []
dej = []
deh = []
dek = []
deir1 = []
deir2 = []
deir3 = []
deir4 = []
dmjd = []


## Converting the gloess fourtran/pgplot code to python/matplotlib
## June 15 2012

## Version 1.0
## last edit - June 19 2012

## Next thing to add:
##Print fits to an output text file


## Open the input data file and read the info

input = sys.argv[1]
counter = 0

## Want to know whether the IRAC data is phased or not. 
## If it is phased, must reduce the uncertainty by another factor of sqrt(N)
## if phased == 1 then true. if phased == 0, false

print input

for line in open(input):
	data = line.split()
	if counter == 0:	
		cepname = data[0]
	if counter == 1:
		period = float(data[0])
		if period > 0:
			phased = 1
		else:
			phased = 0
	if counter == 2:
		nlines = float(data[0])
	if counter == 3:
		xu = float(data[0])
		xb = float(data[1])
		xv = float(data[2])
		xr = float(data[3])
		xi = float(data[4])
		xj = float(data[5])
		xh = float(data[6])
		xk = float(data[7])
		xir1 = float(data[8])
		xir2 = float(data[9])
		xir3 = float(data[10])
		xir4 = float(data[11])
	if counter > 3:
		dmjd.append(float(data[0]))
		du.append(float(data[1]))
		deu.append(float(data[2]))
		db.append(float(data[3]))
		deb.append(float(data[4]))
		dv.append(float(data[5]))
		dev.append(float(data[6]))
		dr.append(float(data[7]))
		der.append(float(data[8]))
		di.append(float(data[9]))
		dei.append(float(data[10]))
		dj.append(float(data[11]))
		dej.append(float(data[12]))
		dh.append(float(data[13]))
		deh.append(float(data[14]))
		dk.append(float(data[15]))
		dek.append(float(data[16]))
		dir1.append(float(data[17]))
		deir1.append(float(data[18]))
		dir2.append(float(data[19]))
		deir2.append(float(data[20]))
		dir3.append(float(data[21]))
		deir3.append(float(data[22]))
		dir4.append(float(data[23]))
		deir4.append(float(data[24]))		
	counter  = counter + 1	
		
## Read in all the data from the file and filled the arrays. Need to convert these to numpy arrays.

number = counter - 4 # Number data lines in the file
#print number

u = np.array(du)
b = np.array(db)
v = np.array(dv)
r = np.array(dr)
i = np.array(di)
j = np.array(dj)
h = np.array(dh)
k = np.array(dk)
ir1 = np.array(dir1)
ir2 = np.array(dir2)
ir3 = np.array(dir3)
ir4 = np.array(dir4)
eu = np.array(deu)
eb = np.array(deb)
ev = np.array(dev)
er = np.array(der)
ei = np.array(dei)
ej = np.array(dej)
eh = np.array(deh)
ek = np.array(dek)
eir1 = np.array(deir1)
eir2 = np.array(deir2)
eir3 = np.array(deir3)
eir4 = np.array(deir4)
mjd = np.array(dmjd)

nu = sum(u<50)
nb = sum(b<50)
nv = sum(v<50)
nr = sum(r<50)
ni = sum(i<50)
nj = sum(j<50)
nh = sum(h<50)
nk = sum(k<50)
nir1 = sum(ir1<50)
nir2= sum(ir2<50)
nir3= sum(ir3<50)
nir4= sum(ir4<50)

# Phases don't need to be done individually by band - only depends on P
phase = (mjd / period) - np.floor(mjd / period)
phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

# Usage:  fit_one_band(data,err,phases,n,smooth):
maxvals = []
minvals = []
if nu > 0:
	maxvals.append(np.amax(u[u<50])+3.0)
	minvals.append(np.amin(u[u<50])+3.0)
if nb > 0:
	maxvals.append(np.amax(b[b<50])+1.5)
	minvals.append(np.amin(b[b<50])+1.5)
if nv > 0:
	maxvals.append(np.amax(v[v<50])+1.2)
	minvals.append(np.amin(v[v<50])+1.2)
if nr > 0:
	maxvals.append(np.amax(r[r<50])+0.7)
	minvals.append(np.amin(r[r<50])+0.7)
if ni > 0:
	maxvals.append(np.amax(i[i<50])+0.2)
	minvals.append(np.amin(i[i<50])+0.2)
if nj > 0:
	maxvals.append(np.amax(j[j<50]))
	minvals.append(np.amin(j[j<50]))
if nh > 0:
	maxvals.append(np.amax(h[h<50])-0.4)
	minvals.append(np.amin(h[h<50])-0.4)
if nk > 0:
	maxvals.append(np.amax(k[k<50])-0.8)
	minvals.append(np.amin(k[k<50])-0.8)
if nir1 > 0:
	maxvals.append(np.amax(ir1[ir1<50])-1.4)
	minvals.append(np.amin(ir1[ir1<50])-1.4)
if nir2 > 0:
	maxvals.append(np.amax(ir2[ir2<50])-1.8)
	minvals.append(np.amin(ir2[ir2<50])-1.8)
if nir3 > 0:
	maxvals.append(np.amax(ir3[ir3<50])-2.2)
	minvals.append(np.amin(ir3[ir3<50])-2.2)
if nir4 > 0:
	maxvals.append(np.amax(ir4[ir4<50])-2.6)
	minvals.append(np.amin(ir4[ir4<50])-2.6)


maxvals = np.array(maxvals)
minvals = np.array(minvals)

max = np.max(maxvals)
min = np.min(minvals)
print cepname, ' ---- Period =', period, 'days'
print '------------------------------------------------------'

# Set up names for output files

#fitname = cepname + '.glo_fits'
avname = cepname + '.glo_avs'

avsout = open(avname,'w')
#fitout = open(fitname,'w')

maxlim = max + 0.5
minlim = min - 0.5



mp.clf()

#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#mp.figure(figsize=(16.0,10.0))


gs = gridspec.GridSpec(3, 3)
#ax1 = plt.subplot(gs[:, 0:2])
ax2 = plt.subplot(gs[0,0:3])
ax3 = plt.subplot(gs[1, 0:3])
ax4 = plt.subplot(gs[2, 0:3])
#ax1.axis([1,3.5,(maxlim),(minlim)])
titlestring = cepname + ', P = ' + str(period) + ' days'
#print titlestring
mp.suptitle(titlestring, fontsize=20)

#ax1.set_ylabel('Magnitude')
#ax1.set_xlabel('Phase $\phi$')


## Fitting and plotting for each band
print nu, nb, nv, nr, ni, nj, nh, nk, nir1, nir2, nir3, nir4

if nir1 > 0:
	ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(ir1,eir1,phase,nir1,xir1)
	#ax1.plot(ir1x,ir11-1.4,'k-')
 	#ax1.plot(xphaseir1,yir1-1.4,color='MediumVioletRed',marker='o',ls='None', label='$[3.6]-1.4$')
## for RRLyrae WISE plots:
#	ax1.plot(ir1x,ir11+1.,'k-')
# 	ax1.plot(xphaseir1,yir1+1.,color='Turquoise',marker='o',ls='None', label='W1+1.0')
	aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)
	if phased == 1:
		factor = sqrt(nir1)
	if phased == 0:
		factor = 1 
	if nir1 > 1:
		print >> avsout, '<[3.6]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f} N I1 = {3}'.format(aveir1, sdevir1/factor, ampir1,nir1)
		print  '<[3.6]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir1, sdevir1/factor, ampir1)
	if nir1 == 1:
		print >> avsout, '[3.6] = {0:.3f} --- single point'.format(aveir1)
		print  '[3.6] = {0:.3f} --- single point'.format(aveir1)

if nir2 > 0:
	ir21, ir2x, yir2, yeir2, xphaseir2 = gf.fit_one_band(ir2,eir2,phase,nir2,xir2)
	#ax1.plot(ir2x,ir21-1.8,'k-')
 	#ax1.plot(xphaseir2,yir2-1.8,color='DeepPink',marker='o',ls='None', label='$[4.5]-1.8$')
## For RRLyrae WISE plots:
#	ax1.plot(ir2x,ir21,'k-')
# 	ax1.plot(xphaseir2,yir2,color='Gold',marker='o',ls='None', label='W2')
	aveir2, adevir2, sdevir2, varir2, skewir2, kurtosisir2, ampir2= gf.moment(ir21[200:300],100)
	if phased == 1:
		factor = sqrt(nir2)
	if phased == 0:
		factor = 1

	if nir2 > 1:
		print >> avsout, '<[4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f} N I2 = {3}' .format(aveir2, sdevir2/factor, ampir2,nir2)
		print '<[4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir2, sdevir2/factor, ampir2)
	if nir2 == 1:
		print >> avsout, '[4.5] = {0:.3f} --- single point'.format(aveir2)
		print '[4.5] = {0:.3f} --- single point'.format(aveir2)


## gloess differential

dy1dp_max_light = 1000
dy1dp_min_light = 1000

for ptop in range(301,400):
	dp = 0.01
	dy1 = ir11[ptop] - ir11[ptop - 1]
	dy2 = ir21[ptop] - ir21[ptop - 2]
	diff = np.abs(dy1/dp)
	if (diff < dy1dp_max_light) and (ir11[ptop] > ir11[ptop -1]):
		dy1dp_max_light = diff
		max_light_val = ptop

for ptop in range(351,400):
	dp = 0.01
	dy1 = ir11[ptop] - ir11[ptop - 1]
	dy2 = ir21[ptop] - ir21[ptop - 2]
	diff = np.abs(dy1/dp)

	if (diff < dy1dp_min_light) and (ir11[ptop] > ir11[ptop -1]):
		dy1dp_min_light = diff
		min_light_val = ptop


print "max light ", max_light_val, dy1dp_max_light, ir11[max_light_val], ir1x[max_light_val]
print "min light ", min_light_val, dy1dp_min_light, ir11[min_light_val], ir1x[min_light_val]

min_phase = ir1x[min_light_val] - floor(ir1x[min_light_val])


ir1x = ir1x - min_phase
ir2x = ir2x - min_phase
xphaseir1 = xphaseir1 - min_phase
xphaseir2 = xphaseir2 - min_phase


### Define the colour curve
colour_curve = ir11 - ir21
## Define the colour points
ch1_points = yir1[yir1<99]
ch2_points = yir2[yir2<99]
colour_points = ch1_points - ch2_points
colour_phases = xphaseir1[yir1<99]

colour_points = np.concatenate((colour_points,colour_points,colour_points,colour_points,colour_points))
colour_phases = np.concatenate((colour_phases,(colour_phases+1.),(colour_phases+2.),(colour_phases+3.),(colour_phases+4.)))


avecol, adevcol, sdevcol, varcol, skewcol, kurtosiscol, ampcol = gf.moment(colour_curve[200:300],100)

print >> avsout, '<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol)
print  '<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol)

print np.average(ir11[200:300]) + 0.3
print np.average(ir11[200:300]) - 0.3

mp.clf()

ax2 = subplot(311)

ax2.axis([0,2.5,(np.average(ir11[200:300]) + 0.4),(np.average(ir11[200:300]) - 0.4)])
#ax2.yaxis.tick_right()
ax2.yaxis.set_major_locator(plt.FixedLocator([10.8, 10.6, 10.4, 10.2]))
ax2.plot(ir1x,ir11,'k-')
ax2.plot(xphaseir1,yir1,color='DodgerBlue',marker='o',ls='None', label='$[3.6]$')
ax2.vlines(ir1x[min_light_val], 11.0, 10.0, 'grey', 'dashed')
ax2.vlines(ir1x[max_light_val], 11.0, 10.0,  'grey', 'dashed')
ax2.vlines(ir1x[min_light_val-100], 11.0, 10.0, 'grey', 'dashed')
ax2.vlines(ir1x[max_light_val-100], 11.0, 10.0,  'grey', 'dashed')
mp.ylabel("[3.6]")
ax2.annotate('Smooth transition \n at 3.6~$\mu$m maximum', xy=(2.0, 10.7), xycoords='data',fontsize=14, ha='center', bbox=dict(boxstyle="round", fc="ivory", ec='lightgrey'))
ax2.arrow(1.69, 10.603, -0.19, -0.213, head_width=0.05, head_length=0.05, fc='dimgrey', ec='dimgrey')



ax5 = ax2.twiny()
ax5.axis([0,2.5,(np.average(ir11[200:300]) + 0.4),(np.average(ir11[200:300]) - 0.4)])

ax5.set_xticks([0.5, 1.0, 1.5, 2.0])
ax5.tick_params(axis='x', which='both', top='off')
ax5.set_xticklabels(['Max [3.6]', 'Min [3.6]', 'Max [3.6]', 'Min [3.6]'])

mp.title("HV 900 P = 47.515 d", y=1.2)

#ax2.yaxis.label("[3.6]")
mp.ylabel("[3.6]")

#ax2.annotate('$[3.6]$', xy=(0.04, 0.8375), xycoords='axes fraction', fontsize=16)
ax3 = subplot(312)

ax3.axis([0,2.5,(np.average(ir21[200:300]) + 0.4),(np.average(ir21[200:300]) - 0.4)])
ax3.yaxis.set_major_locator(plt.FixedLocator([10.8, 10.6, 10.4, 10.2]))
ax3.plot(ir2x,ir21,'k-')
ax3.plot(xphaseir2,yir2,color='DeepPink',marker='o',ls='None', label='$[4.5]$')
#ax3.annotate('$[4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)
ax3.vlines(ir1x[min_light_val], 11.0, 10.0, 'grey', 'dashed')
ax3.vlines(ir1x[max_light_val], 11.0, 10.0,  'grey', 'dashed')
ax3.vlines(ir1x[min_light_val-100], 11.0, 10.0, 'grey', 'dashed')
ax3.vlines(ir1x[max_light_val-100], 11.0, 10.0,  'grey', 'dashed')

ax3.annotate('CO recombination \n begins', xy=(0.85, 10.3), xycoords='data',fontsize=14, ha='center', bbox=dict(boxstyle="round", fc="ivory", ec='lightgrey'))
ax3.arrow(1.08, 10.20, 0.20, 0.15, head_width=0.05, head_length=0.05, fc='dimgrey', ec='dimgrey')


mp.ylabel("[4.5]")

ax4 = subplot(313)

#divider = make_axes_locatable(ax1)
#axcol = divider.append_axes("bottom",1.2,pad=0.1,sharex=ax1)
myaxis2 = [0,2.5,-0.2,0.2]
ax4.axis(myaxis2)
ax4.yaxis.set_major_locator(plt.FixedLocator([-0.1,0,0.1]))
ax4.plot(ir1x,colour_curve,'k-')
ax4.plot(colour_phases,colour_points,color='Black',marker='o',ls='None', label='$[3.6]-[4.5]$')

ax4.set_xlabel('Phase $\phi$')
#ax4.annotate('$[3.6] - [4.5]$', xy=(1.1, 0.135), xycoords='data')
#ax4.annotate('$[3.6] - [4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)

ax4.hlines(0,0,2.5,'k','dashdot')
ax4.vlines(ir1x[min_light_val], -0.3, 0.3, 'grey', 'dashed', zorder=0)
ax4.vlines(ir1x[max_light_val], -0.3, 0.3,  'grey', 'dashed', zorder=0)
ax4.vlines(ir1x[min_light_val-100], -0.3, 0.3, 'grey', 'dashed', zorder=0)
ax4.vlines(ir1x[max_light_val-100], -0.3, 0.3,  'grey', 'dashed', zorder=0)

ax4.fill_between(ir1x, colour_curve, 0, where=colour_curve<0, facecolor='#f8e782', alpha=1.0)

ax4.annotate('CO dissociation \n colour plateau',  xy=(0.6, 0.075), xycoords='data', fontsize=14, ha='center', bbox=dict(boxstyle="round", fc="ivory", ec='lightgrey'))
ax4.arrow(0.86, 0.06, 0.15, -0.02, head_width=0.03, head_length=0.03, fc='dimgrey', ec='dimgrey')



ax4.annotate('CO suppression \n of 4.5~$\mu$m flux', xy=(1.80, 0.075), xycoords='data',fontsize=14, ha='center', bbox=dict(boxstyle="round", fc="ivory", ec='lightgrey'))
ax4.arrow(1.81, 0.06, 0.0, -0.08, head_width=0.03, head_length=0.03, fc='dimgrey', ec='dimgrey')



mp.ylabel("$[3.6]-[4.5]$")

mp.setp(ax2.get_xticklabels(),visible=False)
mp.setp(ax3.get_xticklabels(),visible=False)


#mp.savefig('annotated_cepheid.pdf', transparent='True')

avsout.close()
mp.show()
#fitout.close()
	
	







															
		

		
	


