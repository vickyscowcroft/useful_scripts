""" useful functions to do annoying things"""

import pandas as pd
import pyvo as vo
import re
import numpy as np
from astroquery.gaia import Gaia

def remove_whitespace_from_cols(df):
    """ removes leading and trailing whitespace from dataframe column names """

    df.columns = df.columns.str.strip()
    return(df)

def print_columns(df):
    """ prints all the column names in a dataframe
        useful for when there's too many and df.columns
        doesn't print them all. """

    for i in range(len(df.columns)):
        print(df.columns[i])
    return()

def get_gaia_source_id(name, gaia_release='DR3'):
    """ gets the gaia source id for an object from simbad
        defaults to DR3 source id (identical to EDR3)
    """
    tap_service = vo.dal.TAPService("http://simbad.cds.unistra.fr/simbad/sim-tap")
    query_string = f"SELECT id2.id FROM ident AS id1 JOIN ident AS id2 USING(oidref) WHERE id1.id = \
        '{name}' and id2.id like 'Gaia {gaia_release}%';"
    tap_results = tap_service.search(query_string)
    if len(tap_results) > 0:
        value = tap_results['id', 0]
        source_id = value.split(' ')[2]
    else:
        source_id = np.nan
    return(source_id)

def get_ogle_id(name):
    """ get the ogle id for a star """
    tap_service = vo.dal.TAPService("http://simbad.cds.unistra.fr/simbad/sim-tap")
    query_string = f"SELECT id2.id FROM ident AS id1 JOIN ident AS id2 USING(oidref) WHERE id1.id = \
        '{name}' and id2.id like 'OGLE%';"
    tap_results = tap_service.search(query_string)
    value = tap_results['id']
    pattern = re.compile('-')
    ogle_id = np.nan
    for i in range(len(value)):
        if re.search(pattern, value[i]):
            ogle_id = value[i]
    return(ogle_id)
    
def mu_to_kpc(mu):
	""" convert distance modulus to kpc """
	d_pc = 10**((mu + 5)/5.)
	d_kpc = d_pc / 1000.
	return(d_kpc)
	
def kpc_to_mu(d_kpc):
	""" convert distance in kpc to distance modulus """
	d_pc = d_kpc * 1000.
	mu = 5*np.log10(d_pc) - 5.0
	return(mu)

def read_gaia_epoch_photometry_from_query(source_id):
    """ read in the gaia epoch photometry files from datalink
    and convert to the right type of file for gloess
    """
    retrieval_type = 'EPOCH_PHOTOMETRY'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
    data_structure = 'INDIVIDUAL'   # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
    data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'


    datalink = Gaia.load_data(ids=source_id, data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = False, output_file = None)

    dl_key = 'EPOCH_PHOTOMETRY-Gaia DR3 ' + str(source_id) + '.xml'
    vot_df = datalink[dl_key][0].to_table().to_pandas()
    #vot_df = vot.parse_single_table(filename).to_table().to_pandas()
    if vot_df.source_id.nunique() > 1:
        print('more than one source_id in this file.')
        return(1)
    vot_df.dropna(subset='time', inplace=True)
    piv_df = vot_df[['band', 'time', 'mag', 'flux_over_error', 'source_id']].pivot(index="time", columns="band", values=["mag", 'flux_over_error', 'source_id'])
    """ check it's just a single object """
    
    filters = vot_df.band.unique()
    """ times are JD-2455197.5"""
    names = set_up_dataframe_cols(filters)
    names = np.append(names, 'source_id')
    df = pd.DataFrame(columns=names, index=vot_df.time.dropna().values)
    df['Gaia_JD'] = df.index.copy()
    df['MJD'] = get_gaia_jds(df, jd_col='Gaia_JD')
    for filt in filters:
        mag_col = 'mag_' + filt
        err_col = 'err_' + filt
        df[mag_col] = piv_df[('mag', filt)]
        df[err_col] = piv_df.apply(lambda x: get_gaia_errs(x[('flux_over_error', filt)], filt), axis=1)
    df['source_id'] = vot_df['source_id'][0]
    df.reset_index(inplace=True, drop=True)

    return(df)


    