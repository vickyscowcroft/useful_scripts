import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
from scipy.optimize import curve_fit
import reddening_laws as red

def linear_pl(lp, slope, intercept):
    """ simple linear pl relation
    lp = array of (log10 P - pivot)
    fits function of the form mag = a * (lp) + b
    where a = slope, b = intercept
    
    """
    mag = slope * (lp) + intercept

    return(mag)

def find_pl_coefficients(period, mags, pl_type=linear_pl, errs=np.NaN, weighted=False, pivot=0.0, logp=False):
    """ find coefficiencts for pl relation. 
    can do error weighted fitting if you pass phot errors to errs and weighted=True
    default is unweighted
    """

    if logp == True:
        lp = period - pivot
    else:
        lp = np.log10(period) - pivot

    popt, pcov = curve_fit(linear_pl, lp, mags)
    return(popt, pcov)

def wesenheit_3_band_calc(mag1, mag2, mag3, bands=['G', 'BP', 'RP']):
    """make this general. only works for gaia right now"""
    G_lam = 0.622
    BP_lam = 0.511
    RP_lam = 0.777
    R_CCM = 3.1
    
    A_G = red.ccm_optical(G_lam, R_CCM)
    A_BP = red.ccm_optical(BP_lam, R_CCM)
    A_RP = red.ccm_optical(RP_lam, R_CCM)
    alpha = A_G / (A_BP - A_RP)

    wes_mag = mag1 - alpha*(mag2 - mag3)
    return(wes_mag)


