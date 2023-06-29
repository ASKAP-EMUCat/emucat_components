###cross match positions with external
##Data central
##CDS?
###make command line friendly with votable i/o
import numpy as np
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus

from pyvo.dal import TAPService



######################################################
######################################################
###parameters

qschema = f"SELECT * FROM TAP_SCHEMA.columns"

#infile = 'test_positions.xml'
infile = 'test_positions_gama.xml'

tap_datacentral = 'https://datacentral.org.au/vo/tap'

######################################################
######################################################
###functionality

def print_tables(tap_url):
    'print available tables in TAP client'
    ###set up tap client
    client = TAPService(tap_url)
    
    ###print available tables
    tables = list(client.tables.keys())
    for table in tables:
        print(table)
    return
    
    
def pyvo_query(query, tap_url, upload_data={},
               language='SQL'):
    'ADQL query of public CADC YouCat database'
    
    ###set up tap client
    client = TAPService(tap_url)
    
    ###query data
    data = client.search(query, language=language).to_table()
    
    return(data)


def xmatch_with_data_central(upload_data, acol, dcol,
                             acol_dc='ra', dcol_dc='dec',
                             qtable='dc_conesearch.gama_dr2',
                             searchrad_arcsec=1,
                             timeout=600):
    'perform positional cross match with tables from astro data central'
    client = TAPService("https://datacentral.org.au/vo/tap")
    searchrad_deg = np.round(searchrad_arcsec/3600, 5)
    
    ###make sure upload positional columns dont conflict with those in data central
    if acol==acol_dc:
        acol_new = f'{acol}_upload'
        upload_data.rename_column(name=acol, new_name=acol_new)
        acol = acol_new
    if dcol==dcol_dc:
        dcol_new = f'{dcol}_upload'
        upload_data.rename_column(name=dcol, new_name=dcol_new)
        dcol = dcol_new
    
    ###write ADQL query
    adql = f"SELECT q3c_dist({acol_dc},{dcol_dc}, tup.{acol},tup.{dcol}) *3600 AS angDist, * FROM {qtable}, tap_upload.upload_table as tup WHERE 't' = q3c_radial_query({acol_dc},{dcol_dc},tup.{acol},tup.{dcol}, {searchrad_deg}) ORDER BY angDist ASC"
    
    ##perform query
    uploads = {"upload_table" : upload_data}
#    results = client.search(adql, uploads=uploads).to_table()
    results = client.run_async(adql, uploads=uploads, timeout=timeout).to_table()
    
    return results


######################################################
######################################################
###main

###issues querying table names starting with number, e.g. noirlab_test2.1_allwise, but working for table names not starting with a number -- use "" around table name: noirlab_test2."1_allwise"tq

####use q3cdist with dc_conesearch tables (except total)
###add in CDS functionality too


data = Table.read(infile)


####testing

test1 = xmatch_with_data_central(upload_data=data[:3],
                                 acol='RAJ2000', dcol='DEJ2000',
                                 acol_dc='ra', dcol_dc='dec',
                                 qtable='dc_conesearch."2mass_fdr"',
                                 searchrad_arcsec=1)
