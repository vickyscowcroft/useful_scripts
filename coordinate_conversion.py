#!/usr/bin/env python

import numpy as np


def todeg(ra, dec):	
	if ra.find(":")!= -1 :
		raspl = ra.split(':')
	else:
		raspl = ra.split()
	rah = float(raspl[0])
	ram = float(raspl[1])
	ras = float(raspl[2])
	if dec.find(":") != -1 :
		decspl = dec.split(':')
	else:
		decspl = dec.split()
	decd = float(decspl[0])
	decm = float(decspl[1])
	decs = float(decspl[2])
	
	radeg = (rah + (ram / 60.) + (ras / 3600.))*15.

	## Can't just test for positive or negative for decd
	## This misses the negative zeroes. 

	if str(decd)[0] == '-':
		sign = -1.
	else: 
			sign = 1.
	
	decdeg = decd + (sign * decm / 60.) + (sign * decs / 3600.)


	return(radeg, decdeg)
	
def tosex(ra, dec):
	
	ra = float(ra)
	dec = float(dec)
	rahours = ra /15.
	rah = np.floor(rahours)
	ram = np.floor((rahours - rah)*60.)
	ras = (rahours - rah - (ram/60.))*3600.
	if dec >= 0.0:
		decd = np.floor(dec)
		decm = np.floor((dec - decd)*60.)
		decs = (dec - decd - (decm / 60.0))*3600.
	if dec < 0.0:
		dec = -dec
		decd = np.floor(dec)
		decm = np.floor((dec - decd)*60.)
		decs = (dec - decd - (decm / 60.0))*3600.
		decd = -decd
	rasex = "{0:.0f}:{1:.0f}:{2:.0f}".format(rah, ram, ras)
	decsex = "{0:.0f}:{1:.0f}:{2:.0f}".format(decd, decm, decs)
	#rasex = str(rah) + ":" + str(ram) + ":" + str(ras)
	#decsex = str(decd) + ":" + str(decm) + ":" + str(decs)
	
	return(rasex, decsex)
	