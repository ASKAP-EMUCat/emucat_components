###cross match positions with external
##Data central
##CDS?
###make command line friendly with votable i/o
import numpy as np
from astropy.table import Table
#from astroquery.utils.tap.core import TapPlus
from astroquery.xmatch import XMatch
from astropy import units as u

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

class searchDC:
    'perform positional cross match with tables from astro data central'
    def __init__ (self, data, acol='RAJ2000', dcol='DEJ2000',
                  acol_dc='ra', dcol_dc='dec', namecol='AllWISE',
                  searchrad_arcsec=1*u.arcsec,
                  timeout=600):
            'initial class setup'
            self.client = TAPService("https://datacentral.org.au/vo/tap")
            self.timeout=timeout
            self.namecol=namecol
            self.searchrad_deg = np.round(searchrad_arcsec.value/3600, 5)
            ###make sure upload positional columns dont conflict with those in data central
            self.acol_dc = acol_dc
            self.dcol_dc = dcol_dc
            self.acol = acol
            self.dcol = dcol
            if acol==acol_dc:
                acol_new = f'{acol}_upload'
                data.rename_column(name=acol,
                                   new_name=acol_new)
                self.acol = acol_new
            if dcol==dcol_dc:
                dcol_new = f'{dcol}_upload'
                data.rename_column(name=dcol,
                                   new_name=dcol_new)
                self.dcol = dcol_new

            ###define uploads
            self.uploads = {"upload_table" : data}
    
    def gama(self):
        'query GAMA'
        adql = f"SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t2.CATAID, t2.NQ, t2.Z FROM dc_conesearch.gama_dr2 AS t1 INNER JOIN gama_dr2.SpecObj AS t2 ON t1.name=t2.CATAID, tap_upload.upload_table as tup WHERE 't' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC"

        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()

        return results_table

    def twodf(self):
        'query 2dFGRS'
#        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t2.serial, t2.quality, t2.z FROM dc_conesearch."2dfgrs_fdr" as t1 INNER JOIN "2dfgrs".spec_best AS t2 ON t1.name=t2.serial, tap_upload.upload_table as tup WHERE '+"'t'"+ f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC' ###need to fix join query
        
        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, * FROM dc_conesearch."2dfgrs_fdr" as t1, tap_upload.upload_table as tup WHERE '+"'t'"+f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC'
        
        print(adql)
        print('')
        
        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()

        return results_table
        

def cds_xmatch(data, racol='RAJ2000', decol='DEJ2000',
               maxsep=1*u.arcsec,
               catcols='*',
               namecol='AllWISE', colsuff=None,
               cat2='vizier:V/154/sdss16',
               timeout=600):
    'use cds xmatch to query sdss'
    ###try to replace with async query (might not need to)
    xm = XMatch()
    xm.TIMEOUT = timeout
    
    xmatch = xm.query(cat1=data, cat2=cat2, max_distance=maxsep,
                      colRA1=racol, colDec1=decol)
    
    ###reduce to only spec data (and unique)
    if catcols != '*':
        outcols = [namecol] + catcols
    else:
        outcols = xmatch.colnames
    
    xmatch.sort('angDist')
    
    return xmatch



######################################################
######################################################
###main
###write specific queries for each dataset and make argument of xmatch_with_data_central()


data = Table.read(infile)
data = data[['AllWISE', 'RAJ2000', 'DEJ2000', 'CATAID']]


####testing

test_data = data[:3]
xdc = searchDC(test_data, searchrad_arcsec=3600*u.arcsec) ###wide search to test output where there are no matches



#client = TAPService("https://datacentral.org.au/vo/tap")
#acol_dc, dcol_dc = 'ra', 'dec'
#acol, dcol = 'RAJ2000', 'DEJ2000'
#qtable = 'dc_conesearch."2mass_fdr"'
#searchrad_deg = 2
#uploads = {"upload_table" : test_data}
#
#adql = f"SELECT q3c_dist({acol_dc},{dcol_dc}, tup.{acol},tup.{dcol}) *3600 AS angDist, * FROM {qtable}, tap_upload.upload_table as tup WHERE 't' = q3c_radial_query({acol_dc},{dcol_dc},tup.{acol},tup.{dcol}, {searchrad_deg}) ORDER BY angDist ASC"
#
#t = "'t'"
#adql = f'SELECT q3c_dist(t1.{acol_dc},t1.{dcol_dc}, tup.{acol},tup.{dcol}) *3600 AS angDist, tup.AllWISE, * FROM dc_conesearch."2dfgrs_fdr" as t1, tap_upload.upload_table as tup WHERE {t} = q3c_radial_query(t1.{acol_dc},t1.{dcol_dc},tup.{acol},tup.{dcol}, {searchrad_deg}) ORDER BY angDist ASC'
#
#atest = 'SELECT q3c_dist(t1.ra,t1.dec, tup.RAJ2000,tup.DEJ2000) *3600 AS angDist, tup.AllWISE, * FROM dc_conesearch."2dfgrs_fdr" as t1, tap_upload.upload_table as tup WHERE '+"'t'"+' =q3c_radial_query(t1.ra,t1.dec,tup.RAJ2000,tup.DEJ2000, 0.00028 arcsec) ORDER BY angDist ASC'


