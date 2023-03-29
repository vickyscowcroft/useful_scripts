""" useful functions to do annoying things"""

import pandas as pd
import pyvo as vo
import re
import numpy as np

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
	

    