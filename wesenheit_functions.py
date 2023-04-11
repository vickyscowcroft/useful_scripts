import numpy as np
import pandas as pd
import reddening_laws as red

""" Create multi-band wesenheit magnitudes"""

def get_filter_wavelength(filt):
    """get the wavelength in microns of the specified filter

    Args:
        filt (str): filter name 
    """
    data_file_path = '/Users/vs522/Dropbox/Python/useful_scripts/'
    lam_df = pd.read_csv(data_file_path + 'filter_wavelengths.csv')
    
    # convert the filter name to upper case to look it up in the table 
    filt = np.char.upper(filt)
    
    wlen = np.array([lam_df.loc[lam_df['band']==f, 'wavelength_um'].item() for f in filt])
# need to fix the error handling here
    # if len(wlen) == 1:
    #     wlen = wlen.item()        
    # else:
    #     raise ValueError("filter doesn't appear in reference table")
    return(wlen)
        
def get_extinction_coeffs(lam, Rv=3.1, ext_law='CCM89'):
    """Get the extinction coefficients, A_lambda for the given wavelengths
    ** Dependencies:
    reddening_laws (my version)
    
    TO DO: Add option to select which extinction law to use. 
    Right now defaults to CCM89 for optical/near-IR
    Indebetouw for IR
    Need to add option for F99 (Gaia uses this)
    Args:
        lambda (float or array): Wavelength(s) that you want the extinction coefficients for
        Rv (float): defaults to 3.1
        ext_law (str): extinction law. Current options are 'CCM89'
    """
    
    ## decide whether to use CCM optical or CCM IR
    
    x_ccm = 1./lam
    ## optical case
    if (x_ccm >= 1.1 and x_ccm <= 3.3):
        AlAv = red.ccm_optical(lam, Rv)
    ## ccm nearir
    elif (x_ccm > 1./1.5 and x_ccm < 1.1):
        AlAv = red.ccm_nearir(lam, Rv)
    ### indebetouw IR
    elif(lam > 1.5):
        AlAk = red.indebetouw_ir(lam)
        Ak = red.ccm_nearir(2.164,Rv)
        AlAv = AlAk * Ak
    else:
        print("wavelength out of range")
        return(-1)
    return(AlAv)
        
#def get_wesenheit_alpha(band_1, band_2, band_3, Rv=3.1):
def get_wesenheit_alpha(bands, Rv=3.1):
    """ Get alpha parameter for wesenheit function
    takes 3 bands as arguments but two of them can be the same

    Args:
        band_1 (str): filter 1
        band_2 (str): filter 2
        band_3 (str): filter 3
        ** testing array of bands ** 
        Rv (float): defaults to 3.1
    """
    # l_1 =  get_filter_wavelength(band_1)
    # l_2 = get_filter_wavelength(band_2)
    # l_3 = get_filter_wavelength(band_3)
    bands = np.array(bands)
    lams = get_filter_wavelength(bands)
    As = np.array([get_extinction_coeffs(l, Rv) for l in lams])

    #As = get_extinction_coeffs(lams, Rv)
    
    # A_1 = get_extinction_coeffs(l_1,Rv = Rv)
    # A_2 = get_extinction_coeffs(l_2,Rv = Rv)
    # A_3 = get_extinction_coeffs(l_3,Rv = Rv)
    
    #alpha =  A_1 / (A_2 - A_3)
    alpha = As[0] / (As[1] - As[2])
    return(alpha)

def wesenheit_3band_mag(mags, bands, Rv=3.1):
    """Get the wesenheit 3-band magnitude for an object, given 3 mag measurements and 3 band names

    Args:
        mags (array[float]): array containing the 3 magnitudes
        bands (array[str]): array containint the 3 band names
        Rv (float, optional): . Defaults to 3.1.
    """
    bands = np.array(bands)
    mags = np.array(mags)
    alpha = get_wesenheit_alpha(bands)
        
    W = mags[0] - alpha*(mags[1] - mags[2])
    return(W)

def wesenheit_3band_mag_err(errs, bands, Rv=3.1):
    """ get the uncertainty on the wesenheit mag
    Need to check my maths here

    Args:
        errs (array):  uncertainties on the mags
        bands (array): bands array
        Rv (float, optional):  Defaults to 3.1.
    """
    alpha = get_wesenheit_alpha(bands, Rv)

    err_W = np.sqrt(errs[0]**2 + (alpha**2)*(errs[1]**2 + errs[2]**2))
    return(err_W)
    
    
    