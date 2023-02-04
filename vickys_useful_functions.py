""" useful functions to do annoying things"""

import pandas as pd
import pyvo as vo

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
    value = tap_results['id', 0]
    source_id = value.split(' ')[2]
    return(source_id)

    