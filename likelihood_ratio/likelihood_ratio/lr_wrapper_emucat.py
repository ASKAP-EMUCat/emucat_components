#!python3

###wrapper for EMUcat likelihood ratio code
###takes VOTable EMUcat files
###converts to work with LoTSS LR code
###runs LR
###converts output back to VOTable for EMUcat


import os
import argparse, subprocess, os, numpy as np, pandas as pd
from astropy.table import Table, join
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from likelihood_ratio_matching import main

###########################################################################
###def functions

def fake_mask(catfile, catformat='votable', acol='ra', dcol='dec',
              outfile='fakemask.fits', overwrite=True):
    
    'fake a mask image (all unmasked) to enable LR code to run using just catalogue info'
    'use AllWISE image parameters, set array values to 1'
    'in lieu of dowloading AllWISE images and mask on the fly'
    'make npx dependent on ra/dec range covered, will work but slow for large regions'
    
    ###obtain ref ra/dec for fake mask from catalogue data
    catdata = Table.read(catfile, format=catformat)

    ra, dec = np.array(catdata[acol]), np.array(catdata[dcol])
    
    ##dec easy, ra needs to account for ra wrapping at 360/0deg
    amin, amax = np.min(ra), np.max(ra)
    dmin, dmax = np.min(dec), np.max(dec)
    
    acen = amin + ((amax-amin)/2)
    dcen = dmin + ((dmax-dmin)/2)
    
    ###account for ra wrapping at 360 - n.b. will break if cone search radius is >90deg!!!
    if amax - amin > 180:
        acen = (360-amax) + ((amin-(360-amax))/2)
    
    ####define number of pixels by wcs world to pix
    ###use lines defining centre of image to determine size (i.e amin and amax versus dcen)
    a0d = SkyCoord(ra=amax, dec=dcen, unit='deg')
    a1d = SkyCoord(ra=amin, dec=dcen, unit='deg')
    d0a = SkyCoord(ra=acen, dec=dmin, unit='deg')
    d1a = SkyCoord(ra=acen, dec=dmax, unit='deg')

    dRA, dDec = a0d.separation(a1d), d0a.separation(d1a)

    ###set npx
    npxa = int((dRA/0.0003819444391411).value)
    npxd = int((dDec/0.0003819444391411).value)

    ###create data array (4095x4095)
    unmasked_array = np.zeros((npxd, npxa))
    
    ##create fits primary HDU
    hdu = fits.PrimaryHDU(unmasked_array)
    
    ###adjust header to privide necessary wcs information
    hdu.header['CTYPE1'] = 'RA---SIN'
    hdu.header['CTYPE2'] = 'DEC--SIN'
    hdu.header['CRVAL1'] = acen
    hdu.header['CRVAL2'] = dcen
    hdu.header['CRPIX1'] = npxa/2
    hdu.header['CRPIX2'] = npxd/2
    hdu.header['CDELT1'] = -0.0003819444391411
    hdu.header['CDELT2'] = 0.0003819444391411
    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CUNIT2'] = 'deg'
    
    ###convert to HDU list and write to file
    hdul = fits.HDUList(hdu)
    
    hdul.writeto(outfile, overwrite=overwrite)
    return


def add_mask_and_stargal_cols(mwfile, fileformat='fits'):
    'add Mask and Stargal columns to multiwavelength data'
    'set all to good for now, update later - stargal can exploit WISE colours'
    
    data = Table.read(mwfile, format=fileformat)

    ###mask 0 = good
    ###sg 1 = good
    maskarray = np.zeros(len(data))
    sgarray = np.ones(len(data))
    
    data['Mask'] = maskarray
    data['stargal'] = sgarray
    
    data.write(mwfile, overwrite=True)
    
    return


def good_data(data, flagcol, flagval, replace=False, param2col='', scalefactor=2,
              warning_col_name='q_warning'):
    ###address issue with flux_int_err == 0
    'either filter out or replace bad data based on column value'
    'assumes data as astropy table'
    
    if replace == True:
        ###arrays to keep easy track of what I'm actually doing here!
        flagcol_array = np.array(data[flagcol])
        param2_array = np.array(data[param2col])
        ##convert bad values
        baddata = (flagcol_array==flagval)
        flagcol_array[baddata] = param2_array[baddata]/scalefactor
        ##e.g. flux_err ==0, if param2 is flux, then this assumes source is scalefactor*sigma
        data[flagcol] = flagcol_array
        ##add warning column
        warning = np.zeros(len(data))
        warning[baddata] = 1
        data[warning_col_name] = warning.astype(int)
    
    else:
        data = data[(data[flagcol]!=flagval)]
        ##add warning col for method consistency (all zeros in this case)
        data[warning_col_name] = np.zeros(len(data)).astype(int)

    return(data)


def VOTable_to_fits(fname_in, fname_out, over_write=True, filterbad=False,
                    flagcol='', flagval='', replace=False, param2col='',
                    flagscale=2):
    ###create fits copy of VOTable from EMUcat for use with LR code
    data = Table.read(fname_in, format='votable')
    for c in data.columns.values():
        if c.dtype == 'object':
            as_str = data[c.name].astype('str')
            data.replace_column(c.name, as_str)

    ##filter/replace bad data (e.g. flux_err == 0)
    if filterbad == True:
        data = good_data(data=data, flagcol=flagcol, flagval=flagval,
                         replace=replace, param2col=param2col,
                         scalefactor=flagscale)

    data.write(fname_out, format='fits', overwrite=over_write)
    return


def output_to_VOTable(fname_in, fname_out, overwrite=True,
                      merge_with_other=False, other_table='',
                      other_format='votable', main_on='', other_on='',
                      needcolumns=[]):
    'take output data and convert to votable for emucat'
    
    data = Table.read(fname_in, format='ascii')
    
    ###if need to merge with other data
    if merge_with_other == True:
        other_data =  Table.read(other_table, format=other_format)
        ###ensure same column name in both tables
        needcolumns = [main_on] + needcolumns ##ensures join key in subset
        if main_on == other_on:
            data = join(data, other_data[needcolumns], keys=main_on, join_type='left')
        else:
            other_data[main_on] = other_data[other_on]
            data = join(data, other_data[needcolumns], keys=main_on, join_type='left')
    
    data.write(fname_out, format='votable', overwrite=overwrite)
    return


def run_lr(fakefile='fakemask.fits', racol='ra', deccol='dec'):
    'run wrapped likelihood ratio code'
    'io format: VOTables'
    
    ##parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--mwcat', type=str, action='store',
                        help='multiwavelength data VOTable to convert to fits for LR code')
    parser.add_argument('--radcat', type=str, action='store',
                        help='radio data VOTable to convert to fits for LR code')
    parser.add_argument('--config', type=str, action='store', default='./lr_config.txt',
                        help='configuration file.')
    args = parser.parse_args()
    
    ###extract info from config file - need outdir, radio ID, flux and err cols
    config_info = pd.read_table(args.config, sep='\s+').replace("'", "", regex=True)
    outdir = config_info[(config_info['parameter']=='outdir')].iloc[0]['value']
    radio_id = config_info[(config_info['parameter']=='radio_id_col')].iloc[0]['value']
    radio_s_col = config_info[(config_info['parameter']=='flux_col')].iloc[0]['value']
    radio_s_err_col = config_info[(config_info['parameter']=='flux_err_col')].iloc[0]['value']

    try:
        os.makedirs(outdir)
    except:
        pass

    ###create fits tables for mw/radcat
    table_data = [args.mwcat, args.radcat]
    fitsnames = [os.path.splitext(fname)[0]+'.fits' for fname in table_data]

    for i in range(len(table_data)):
        if i == 1:
            VOTable_to_fits(fname_in=table_data[i], fname_out=fitsnames[i],
                            filterbad=True, flagcol=radio_s_err_col, flagval=0,
                            replace=True, param2col=radio_s_col, flagscale=2)
        else:
            VOTable_to_fits(fname_in=table_data[i], fname_out=fitsnames[i])
    
    ###add in mask and stargal columns
    add_mask_and_stargal_cols(mwfile=fitsnames[0])

    fake_file = os.path.join(outdir, fakefile)

    ###need to place here what to do about mask image, currently assumes present
    fake_mask(catfile=args.mwcat, acol=racol, dcol=deccol, outfile=fake_file)
    mask = fake_file

    ###run LR code
    main(multiwave_cat=fitsnames[0], radio_cat=fitsnames[1], mask_image=mask, config_file=args.config, overwrite=True,
         snr_cut=5, LR_threshold=0.8)

    ###convert outputs (keep those that end in '_LR_matches.dat') to VOTable
    ###add merging tables to get q_warning

    for filename in os.listdir(outdir):
        if '_LR_matches.dat' in filename:
            newname = os.path.join(outdir, filename.split('.')[0] + '.xml')
            output_to_VOTable(fname_in=os.path.join(outdir, filename), fname_out=newname,
                              overwrite=True, merge_with_other=True, other_table=fitsnames[1],
                              other_format='fits', main_on='radio_ID', other_on=radio_id,
                              needcolumns=['q_warning'])

    ##tidy up - get rid of surplus (add in moving files to target directory)
    '''for filename in file_list:
        if filename not in keep_files:
            os.remove(outdir+filename)

    ###repeat for working directory if outdir != ''
    if outdir!='':
        working_directory_list = os.listdir()
        for filename in working_directory_list:
            if filename not in keep_files:
                os.remove(filename)'''

    return

###########################################################################

if __name__ == "__main__":
    run_lr()

###to do:
##1) upgrade to grab and mosaic WISE images
##2) test with multi-band data
##3) improve stargal and mask flagging in mwcat



