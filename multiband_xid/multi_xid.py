###cross match positions with external
##Data central
##CDS?
###make command line friendly with votable i/o
import numpy as np, argparse
from distutils.util import strtobool
from astropy.table import Table, unique
from astroquery.xmatch import XMatch
from astropy import units as u
from pyvo.dal import TAPService


######################################################
######################################################
###functionality

def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="query external databases for multi-band matches to EMU")
    parser.add_argument("targets",
                        help="catalog of positions to query")
    parser.add_argument("--racol", action='store',
                        type=str, default='RAJ2000',
                        help="RA column in targets")
    parser.add_argument("--decol", action='store',
                        type=str, default='DEJ2000',
                        help="Decl. column in targets")
    parser.add_argument("--namecol", action='store',
                        type=str, default='AllWISE',
                        help="ID column in targets")
    parser.add_argument("--timeout", action='store',
                        type=int, default=600,
                        help="query timeout limit")
    parser.add_argument("--search_rad", action='store',
                        type=str, default='1arcsec',
                        help="output directory")
    parser.add_argument("--only_find_closest", action='store',
                        type=str, default='True',
                        help="output directory")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="output directory")
    
    args = parser.parse_args()
    
    ###make strings quantities where appropriate
    args.search_rad = u.Quantity(args.search_rad)
    args.only_find_closest = strtobool(args.only_find_closest)
    args.outdir = args.outdir.removesuffix('/')
    
    return args


class searchDC:
    'perform positional cross match with tables from astro data central'
    def __init__ (self, data, acol='RAJ2000', dcol='DEJ2000',
                  acol_dc='ra', dcol_dc='dec', namecol='AllWISE',
                  searchrad_arcsec=1*u.arcsec,
                  timeout=600, closest_only=True):
            'initial class setup'
            self.client = TAPService("https://datacentral.org.au/vo/tap")
            self.timeout = timeout
            self.namecol = namecol
            self.searchrad_deg = np.round(searchrad_arcsec.value/3600, 5)
            self.closest_only =  closest_only
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
    
    def twomass(self):
        'queries 2mass for psf and Kron magnitudes (and uncertainties) in J, H and K bands'
        
        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t1.name, t2.j_m, t2.j_cmsig, t2.h_m, t2.h_cmsig, t2.k_m, t2.k_cmsig, t3.j_m_e, t3.j_msig_e, t3.h_m_e, t3.h_msig_e, t3.k_m_e, t3.k_msig_e FROM dc_conesearch."2mass_fdr" as t1 LEFT JOIN "2mass_fdr".psc as t2 ON t1.name=t2.designation LEFT JOIN "2mass_fdr".xsc as t3 ON t1.name=t3.designation, tap_upload.upload_table as tup WHERE '+"'t'"+f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC'
                
        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()
        if self.closest_only==True and len(results_table)>1:
            results_table = unique(results_table, self.namecol,
                                   keep='first')
        
        return results_table
    
    def gama(self):
        'queries GAMA for redshift and spec quality'
        
        adql = f"SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t2.CATAID, t2.NQ, t2.Z FROM dc_conesearch.gama_dr2 AS t1 LEFT JOIN gama_dr2.SpecObj AS t2 ON t1.name=t2.CATAID, tap_upload.upload_table as tup WHERE 't' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC"

        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()
        if self.closest_only==True and len(results_table)>1:
            results_table = unique(results_table, self.namecol,
                                   keep='first')

        return results_table

    def twodf(self):
        'query 2dFGRS for redshift'
        
        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t1.name, t2.Z FROM dc_conesearch."2dfgrs_fdr" as t1 LEFT JOIN "2dfgrs_fdr".spec_all as t2 ON t1.name=t2.name, tap_upload.upload_table as tup WHERE '+"'t'"+f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC'
        
        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()
        if self.closest_only==True and len(results_table)>1:
            results_table = unique(results_table, self.namecol,
                                   keep='first')
                                   
        return results_table
        
    def sixdf(self):
        'query 6dFGS for redshift'
        
        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t1.name, t2.Z FROM dc_conesearch."6dfgs_fdr" as t1 LEFT JOIN "6dfgs_fdr".SPECTRA as t2 ON t1.name=t2.TARGETNAME, tap_upload.upload_table as tup WHERE '+"'t'"+f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC'
                
        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()
        if self.closest_only==True and len(results_table)>1:
            results_table = unique(results_table, self.namecol,
                                   keep='first')
                                   
        return results_table
        
    def wigglez(self):
        'query WiggleZ for redshift and quality'
        
        adql = f'SELECT q3c_dist(t1.{self.acol_dc},t1.{self.dcol_dc}, tup.{self.acol},tup.{self.dcol}) *3600 AS angDist, tup.{self.namecol}, t1.name, t2.redshift, t2.Q  FROM dc_conesearch.wigglez_final as t1 LEFT JOIN wigglez_final.WiggleZCat as t2 ON t1.name=t2.WiggleZ_Name, tap_upload.upload_table as tup WHERE '+"'t'"+f' = q3c_radial_query(t1.{self.acol_dc},t1.{self.dcol_dc},tup.{self.acol},tup.{self.dcol}, {self.searchrad_deg}) ORDER BY angDist ASC'
                
        results_table = self.client.run_async(adql,
                                              uploads=self.uploads,
                                              timeout=self.timeout).to_table()
        if self.closest_only==True and len(results_table)>1:
            results_table = unique(results_table, self.namecol,
                                   keep='first')
                                   
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


if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.targets)
    qdatacen = searchDC(data=data, acol=args.racol, dcol=args.decol,
                        namecol=args.namecol,
                        searchrad_arcsec=args.search_rad,
                        timeout=args.timeout,
                        closest_only=args.only_find_closest)
    twomass = qdatacen.twomass()
    sixdf = qdatacen.sixdf()
    twodf = qdatacen.twodf()
    gama = qdatacen.gama()
    wigglez = qdatacen.wigglez()
    
    ###write to file
    twomass.write(f'{args.outdir}/x_2MASS.xml', format='votable')
    sixdf.write(f'{args.outdir}/x_6dFGS.xml', format='votable')
    twodf.write(f'{args.outdir}/x_2dFGRS.xml', format='votable')
    gama.write(f'{args.outdir}/x_GAMA.xml', format='votable')
    wigglez.write(f'{args.outdir}/x_WiggleZ.xml', format='votable')
    
