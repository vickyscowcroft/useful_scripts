#!/usr/bin/env python

import numpy as n
import sys
import matplotlib.pyplot as mp
import os
import reddening_laws
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'

#matplotlib.rc('ps',usedistiller='xpdf')
matplotlib.rc('text',usetex=True)
rcParams['font.serif'] = ['Garamond']
font = {'fontname':'Garamond'}

mp.clf()
mp.close('all')

nbands = 9
wlength = []
mu = []
err = []
invwlen = []
name = []
Alam = []
Rv = 3.1
extinc = []
inverr = []
# Need Ak / A(V) for Indebetouw conversion later
#Want to get the IR ones in terms of Av instead of Ak
# Indebetouw definition of K filter is 2.164 um

Ak = reddening_laws.ccm_nearir(2.164,3.1)

#input = "/Users/vs/Dropbox/IC1613/wavelength_mu_data"
#input = "/Users/vs/Dropbox/IC1613/template_fitted_data"
#input = "/Users/vs/Dropbox/IC1613/auracaria_data"
input = "/Users/vs/Dropbox/SMC/av_mags/SMC_reddening_data"

##
#This is the file with the r band sample
#input = "/Users/vs/Dropbox/SMC/av_mags/SMC_reddening_data_rband"


#input = "/Users/vs/NGC6822/extinction_curve_gloess_m09"

c = 0
for line in open(input):
	data = line.split()
	name.append(data[0])
	wlength.append(float(data[1]))
	mu.append(float(data[2]))
	err.append(float(data[3]))
	inv = 1.0 / wlength[c]
	inverr.append(1.0 / (err[c]**2))
	invwlen.append(inv)
	#print wlength[c], mu[c], err[c]
	c = c + 1

# Define CCM and Indebetouw extinction laws here
# CCM optical
for band in range(0,nbands):
	if wlength[band] < 1.0:
		#print name[band], " CCM optical"
		AwAv = reddening_laws.ccm_optical(wlength[band],Rv)
		#print name[band], AwAv
		Alam.append(AwAv)
		
## CCM near IR			
	if wlength[band] > 1.0 and wlength[band] < 3.0:
		#print name[band], " CCM infrared"
		AwAv = reddening_laws.ccm_nearir(wlength[band],Rv)
		#print name[band], AwAv
		Alam.append(AwAv)


# Indebetouw IR
	if wlength[band] > 3.0:
		#print name[band], "Indebetouw infrared"
		AwAk = reddening_laws.indebetouw_ir(wlength[band])
		#print name[band], AwAk
		AwAv = AwAk * Ak
		Alam.append(AwAv)
	#print Alam[band]
## Now guess a value of Av to see if it fits. Need to multiply all Alam
## values by that number to get the extinction.


# Guess a value of Av?
#Av = 0.00
#stddevmin = 1000
#sdev = []
#avals = n.arange(0.5,1.5,0.001)
#for c in range(0,1000):

#	mu0 = []
#	for b in range(0,nbands):
#		A = Alam[b] * Av
#		mu0.append(mu[b] - A)
#	avmu = n.average(mu0,weights=(inverr))
	#print avmu
#	avdmu = mu0 - avmu
#	stddev = n.std(avdmu)
	#print stddev
#	sdev.append(stddev)
#	if stddev < stddevmin:
#		Avmin = Av
#		stddevmin = stddev
#		mumin = avmu

#	Av = Av + 0.001

mumin, Avmin, muminerr, Averr = reddening_laws.fit_reddening(mu[0], mu[1], 0.0, mu[3], mu[4], mu[5], mu[6], mu[7], 0.0)

if len(sys.argv) == 2:
	print "Plotting reddedning law for Av = ",sys.argv[1]
	Avmin = float(sys.argv[1])
	mu0 = []
	for b in range(0,nbands):
		A = Alam[b] * Avmin
		mu0.append(mu[b] - A)
	mumin = n.average(mu0,weights=(inverr))

Ebv = Avmin / Rv


print "Rv = ", Rv, " E(B-V) = ", Ebv, " Av = ", Avmin, " mu = ", mumin, "Averr = ", Averr, "muerr = ", muminerr
Atrue = []
mugood = []

for d in range(0,nbands):
	Atrue.append((Alam[d] * Avmin) + mumin)
	mugood.append(mu[d] - (Alam[d]*Avmin))
	
# Calculate 1 sigma errors
# stddev = stddevmin

evi = (Alam[1]*Avmin) - (Alam[2]*Avmin)

#print invwlen_theor

#print reddening_curve_good

Ap1sig = Avmin + 2*Averr
Am1sig = Avmin - 2*Averr
print Am1sig, Ap1sig
Alow = []
Ahigh = []
mulow = []
muhigh = []
inverrcut = []
for e in [0, 1,  3, 4, 5, 6, 7]:
	Al = Alam[e] * Am1sig
	Ah = Alam[e] * Ap1sig
	mulow.append(mu[e] - Al)
	muhigh.append(mu[e] - Ah)
	inverrcut.append(inverr[e])

mulowav = n.average(mulow,weights=(inverrcut))
muhighav = n.average(muhigh,weights=(inverrcut))

print mulowav, muhighav
for f in range(0,nbands):	
	Alow.append((Alam[f] * Am1sig) + mulowav)
	Ahigh.append((Alam[f] * Ap1sig) + muhighav)
	
#print Alow
	
#musig = sqrt(((mulowav - muhighav) / 2.0) )/ sqrt(8.0)
mugood = np.array(mugood)
musig  = mugood.std()
print musig

mp.clf()

ax1 = subplot(211)
myaxis = [0,2.4,18.8,19.4]
mp.axis(myaxis)

ax1.errorbar(invwlen,mu,yerr=err,ls='None',c='k')
ax1.plot(invwlen,mu,'ko',ms=5)
ax1.plot(invwlen,Atrue)
ax1.errorbar(invwlen,mugood,yerr=err,ls='None',c='r')
ax1.plot(invwlen,mugood,'ro',ms=5)
ax1.hlines(mumin,0,2.5)
ax1.plot(invwlen,Alow,'k--')
ax1.plot(invwlen,Ahigh,'k--')
mp.xlabel('1 / $\lambda (\mu m^{-1})$')
mp.ylabel('$\mu$')
#mp.xlabel('1 / $\lambda (\mu m^{-1})$')

#ax2 = subplot(212)
#ax2.plot(avals,sdev)
#mp.xlabel('$A_{V}$')
#mp.ylabel('$\sigma$')

#mp.savefig('inverse_w_mu.eps')
mp.clf()

## Plot for publication

fig = plt.figure()
ax = fig.add_subplot(111)

myaxis = [0,2.4,18.83, 19.35]
ax.axis(myaxis)
ax.yaxis.set_major_locator(plt.FixedLocator([18.9,19.0,19.1, 19.2, 19.3]))
mutext = r"\langle \mu_0 \rangle"
print mutext
annotation = 'SMC Cepheids \n \n$E(B-V) = {0:.3f} \pm {1:.3f} $\n \n $ {4} = {2:.2f} \pm {3:.2f}$'.format(Ebv, (Averr / 3.1) , mumin, muminerr,mutext)
#print >> annotation, '$E(B-V) = {0:.2f} $\n$\mu = {1:.2f} \pm {2:.2f} mag$ \n'.format(Ebv, mumin, musig)

#annotation = '$A_{V} = ' + str(Avmin) + '$\n$\mu = ' + str(mumin) + ' \pm ' + str(musig) + ' mag$ \n'
print annotation

ax.errorbar(invwlen,mu,yerr=err,ls='None',c='k')
ax.plot(invwlen,mu,'ko',ms=5)
## Greying out points that weren't used in the fit
ax.plot(invwlen[2], mu[2], 'white',marker='o',ms=5,markeredgecolor='black')
ax.errorbar(invwlen[2],mu[2],yerr=err[2],ls='None',c='black')
ax.plot(invwlen[8], mu[8], 'white',marker='o',ms=5,markeredgecolor='black')
ax.errorbar(invwlen[8],mu[8],yerr=err[8],ls='None',c='black',)
ax.plot(invwlen,Atrue,'k-')
ax.hlines(mumin,0,2.5,'k','dashdot')
ax.plot(invwlen,Alow,'k--')
ax.plot(invwlen,Ahigh,'k--')
mp.ylabel('$\mu$')
#ax.yaxis.set_major_locator(plt.FixedLocator([24.2,24.4,24.6,24.8]))
ax.annotate(annotation,xy=(0.5, 19.2), xycoords='data', ha='center', bbox=dict(boxstyle="round", fc="w"))
"""ax.annotate('$B$',xy=((1.0/0.445)-0.025,18.855),xycoords='data')
ax.annotate('$V$',xy=((1.0/0.551)-0.025,18.855),xycoords='data')
ax.annotate('$R$',xy=((1.0/0.658)-0.025,18.855),xycoords='data')
ax.annotate('$I$',xy=((1.0/0.806)-0.025,18.855),xycoords='data')
ax.annotate('$J$',xy=((1.0/1.22)-0.025,18.855),xycoords='data')
ax.annotate('$H$',xy=((1.0/1.63)-0.025,18.855),xycoords='data')
ax.annotate('$K_S$',xy=((1.0/2.19)-0.025,18.85),xycoords='data')
ax.annotate('$[3.6]$',xy=((1.0/3.545)-0.01,18.85),xycoords='data')
ax.annotate('$[4.5]$',xy=((1.0/4.442)-0.10,18.85),xycoords='data')"""


ax3 = ax.twiny()
## make the ticks at the band passes

ax3.axis(myaxis)

ax3.set_xticks([invwlen[0], invwlen[1], invwlen[2], invwlen[3], invwlen[4], invwlen[5], invwlen[6], invwlen[7], invwlen[8]])
plt.setp(ax3.get_xticklabels(), visible=False)
## make the labels as the band passes
ax4 = ax.twiny()
ax4.axis(myaxis)

ax4.set_xticks([invwlen[0], invwlen[1], invwlen[2], invwlen[3], invwlen[4], invwlen[5], invwlen[6], invwlen[7]+0.026, invwlen[8]-0.026])
ax4.tick_params(axis='x', which='both', top='off')
ax4.set_xticklabels(['$B$', '$V$', '$R$', '$I$', '$J$', '$H$', '$K$', '[3.6]', '[4.5]'])

mp.setp(ax.get_xticklabels(),visible=False)

def form2(x, pos):
    """ This function returns a string with 2 decimal places, given the input x"""
    return '%.2f' % x

from matplotlib.ticker import FuncFormatter
formatter = FuncFormatter(form2)

divider = make_axes_locatable(ax)
axcorr = divider.append_axes("bottom",1.2,pad=0.1,sharex=ax)
myaxis2 = [0,2.4,mumin-5*musig,mumin+5*musig]
axcorr.axis(myaxis2)
axcorr.yaxis.set_major_locator(plt.FixedLocator([mumin-3*musig, mumin, mumin+3*musig]))

## Have to put the next line after the tick locator definition to make sure it has numbers to format
gca().yaxis.set_major_formatter(FuncFormatter(formatter))

#axcorr.yaxis.set_major_locator(plt.FixedLocator([24.2,24.3,24.4]))
axcorr.plot(invwlen,mugood,'ko',ms=5)
axcorr.errorbar(invwlen,mugood,yerr=err,ls='None',c='k')
axcorr.plot(invwlen[2], mugood[2], 'white',marker='o',ms=5,markeredgecolor='black')
axcorr.errorbar(invwlen[2],mugood[2],yerr=err[2],ls='None',c='black')
axcorr.plot(invwlen[8], mugood[8], 'white',marker='o',ms=5,markeredgecolor='black')
axcorr.errorbar(invwlen[8],mugood[8],yerr=err[8],ls='None',c='black')

axcorr.hlines(mumin,0,2.5,'k','-')
axcorr.hlines(mumin+2*musig,0,2.5,'k','--')
axcorr.hlines(mumin-2*musig,0,2.5,'k','--')

mp.xlabel('$1 / \lambda (\mu m^{-1})$')
mp.ylabel('$\mu_{0}$')

mp.show()
mp.savefig('smc_reddening.pdf',transparent='True')

for band in range(0,nbands):
	print name[band], mugood[band]

print evi
#reddening_curve_good = []
#invwlen_theor = []
#reddening_curve_p1sig = []
#reddening_curve_m1sig = []
#for z in range(445,1000):
#	z = float(z) / 1000.
#	AwAv = reddening_laws.ccm_optical(z,Rv)
#	reddening_curve_good.append((AwAv*Avmin) + mumin)
#	invwlen_theor.append(1./z)
#	reddening_curve_p1sig.append((AwAv*Avmin+stddevmin) + mumin)
#	reddening_curve_m1sig.append((AwAv*Avmin-stddevmin) + mumin)
#for z in range(100,300):
#	z = float(z) / 100.
#	AwAv = reddening_laws.ccm_nearir(z,Rv)
#	reddening_curve_good.append((AwAv*Avmin) + mumin)
#	invwlen_theor.append(1./z)
#	reddening_curve_p1sig.append((AwAv*Avmin+stddevmin) + mumin)
#	reddening_curve_m1sig.append((AwAv*Avmin-stddevmin) + mumin)
#for z in range(300, 450):
#	z = float(z) / 100.
#	AwAv = reddening_laws.indebetouw_ir(z)
#	reddening_curve_good.append((AwAv*Avmin*Ak) + mumin)
#	invwlen_theor.append(1./z)
#	reddening_curve_p1sig.append((AwAv*Avmin*Ak+stddevmin) + mumin)
#	reddening_curve_m1sig.append((AwAv*Avmin*Ak-stddevmin) + mumin)
#ax.plot(invwlen_theor,reddening_curve_good, 'k-')
#ax.plot(invwlen_theor,reddening_curve_p1sig, 'k--')
#ax.plot(invwlen_theor,reddening_curve_m1sig, 'k--')

