###PLACEHODER - YG to update to working placeholder in next few weeks
###forced photometry script for EMU
###takes list of positions (e.g. from WISE), image, and rms image as inputs
###obtains photometry at those positions using aperture defined by beam
###should output measurements as votable


import numpy as np, pandas as pd
from astropy.io import fits
from astropy.table import Table, join
from astropy.wcs import WCS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.psf import EPSFBuilder, EPSFFitter, extract_stars, DAOGroup, BasicPSFPhotometry
from photutils.background import StdBackgroundRMS
#from photutils import aperture
import argparse


########################################
#time of code check - not important for code to run - remove later
import time
start_time = time.time()
########################################


########################################################################################
########################################################################################
###parameters to load from a config file (at some point)


funit_out = 'mJy'
incat_format = 'fits'
racol = 'ra'
deccol = 'dec'
targetnamecol = 'designation'
outputformat = 'fits'
overwrite = False

###make unit astropy unit
funit_out = u.Unit(funit_out)


########################################################################################
########################################################################################

def emu_image_as_2d(file):
    'takes EMU 4D file and outputs 2D HDU for ease of working with astropy on flat image'
    
    hdu = fits.open(file)
    
    ###extract 2D image array
    imdata = hdu[0].data[0][0]
    
    ###extract header info for 2D
    head2d = hdu[0].header.copy()
    
    hkeys = list(head2d.keys())
    
    crkeys = ['CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CUNIT']
    cr3 = [c + '3' for c in crkeys]
    cr4 = [c + '4' for c in crkeys]
    badkeys = cr3 + cr4 + ['NAXIS3', 'NAXIS4']
    
    for key in hkeys:
        if 'PC03' in key or 'PC04' in key or '_03' in key or '_04' in key:
            badkeys.append(key)
        if key in badkeys:
            del(head2d[key])

    head2d['NAXIS'] = 2

    ###create new 2D hdu
    hdu2d = fits.PrimaryHDU(imdata)
    hdu2d.header = head2d
    hdulist2d = fits.HDUList(hdu2d)
    
    return(hdulist2d)


def divide_by_pixels_area(flux, bmaj, bmin, pxunit=u.arcsec, pxsize=1,
                          noise=False, aperture_correction_factor=1):
    'aperture_sum is the integral of flux density over pixels [mJy*px]'
    'pixel units are in mJy/beam, need to convert mJy*px/beam -> mJy/beam'
    'as aperture == beam resultant flux units should thus be in mJy'
    'i.e. mJy = mJy*px/beam * beam/px => mJy*px/n_px'
    
    ###convert beam sizes to pixel units
    bmaj = bmaj.to(pxunit)
    bmin = bmin.to(pxunit) ###swtches to pxunit but doesn't scale
    bmaj = bmaj/pxsize
    bmin = bmin/pxsize
    
    ###beam area = 2pi/8ln(2) * maj * min ~1.133*maj*min
    beam_area = ((2*np.pi)/(8*np.log(2))) * bmaj * bmin
    
    ##if using for correlated noise take root of beam_area
    if noise==True:
        beam_area = np.sqrt(beam_area.value)*beam_area.unit

    ###correct for pixel area to give flux in mJy/beam
    corrected_flux = flux/beam_area

    ###correct for incomplete aperture if required (likely not if using model psf)
    corrected_flux = aperture_correction_factor*corrected_flux
    
    return(corrected_flux)


def model_gaussian_from_beam(hdu):
    'models 2d gaussian in px from beam info'
    ##beam info
    bmin = hdu[0].header['BMIN']*u.deg
    bmaj = hdu[0].header['BMAJ']*u.deg
    bpa = hdu[0].header['BPA']*u.deg
    
    deg_per_px = (hdu[0].header['CDELT2']*(u.deg/u.pixel))
    
    dx, dy = bmaj/deg_per_px, bmin/deg_per_px
    
    ###convert to sigma from full width half max
    dx = dx/(2*np.sqrt(2*np.log(2)))
    dy = dy/(2*np.sqrt(2*np.log(2)))
    
    ###convert bpa to radians and add pi/2 to get to same ref point as image
    da = bpa.to(u.rad) + (np.pi/2)*u.rad
    
    ##obtain Gaussian Kernal
    gauss = Gaussian2DKernel(dx.value, dy.value, da.value)
    
    return(gauss)


def model_psf_from_beam(hdu, osamp=4, iters=10):
    'model psf from the beam geometry'
    ##define Gaussian
    beam_gauss = np.array(model_gaussian_from_beam(hdu))
    x0, y0 = [beam_gauss.shape[1]/2], [beam_gauss.shape[0]/2]
    beampos = Table()
    beampos['x'] = x0
    beampos['y'] = y0
    beam_data = extract_stars(NDData(beam_gauss), beampos, size=np.min(beam_gauss.shape)-2)
    
    ##create epsf
    fitter = EPSFFitter(fitter=LevMarLSQFitter())
    epsf_build = EPSFBuilder(oversampling=osamp, maxiters=iters, fitter=fitter)
    epsf, beam_model = epsf_build(beam_data)
    return(epsf)



def psf_photo(targetfile, imfile, outunit_S=u.mJy,
              acol='RA', dcol='DEC', namecol='Component_name',
              catformat_in='votable', osamp=4, iters=3,
              bkg_sub=None, n_beams=3):
    'use psf photometry to obtain flux'
    
    ####steps:
    ##1) load data
    hdu_im = emu_image_as_2d(imfile)
    targets = Table.read(targetfile, format=catformat_in)
    imdata = hdu_im[0].data
    header = hdu_im[0].header
    
    
    ##2) convert sky coords to pixel positions
    wcs = WCS(header)
    skypos = SkyCoord(ra=targets[acol], dec=targets[dcol], unit='deg')
    pixel_positions = wcs.world_to_pixel(skypos)
    
    tpos_px = Table()
    tpos_px['Name'] = targets[namecol]
    tpos_px['x_0'] = pixel_positions[0]
    tpos_px['y_0'] = pixel_positions[1]
    
    pixelsize_arcsec = abs(header['CDELT1'])*u.Unit(header['CUNIT1']).to(u.arcsec)
    
    ##3) model psf
    epsf = model_psf_from_beam(hdu=hdu_im, osamp=osamp, iters=iters)
    
    ##4) perform PSF photometry
    bkg = bkg_sub #StdBackgroundRMS() or None which produces better results
    
    ###use N beam widths (in pxls) for grouping and aperture radius for init_guesses
    beam_width_pixels = ((header['BMAJ']*u.deg).to(u.arcsec)).value/pixelsize_arcsec
    pxap = int(np.ceil(n_beams*beam_width_pixels))
    pxap = 10
    if pxap%2 == 0:
        fitsize = pxap+1
    else:
        fitsize = pxap

    ###fix psf positions to input
    epsf.x_0.fixed = True
    epsf.y_0.fixed = True

    photometry = BasicPSFPhotometry(group_maker=DAOGroup(pxap), bkg_estimator=bkg,
                                    psf_model=epsf, fitshape=fitsize,
                                    fitter=LevMarLSQFitter(),
                                    aperture_radius=pxap)

    measurements = photometry(imdata, tpos_px)
    
    ##need to convert flux_fit and uncertainty to real fluxes in correct units
    flux_fit = np.array(measurements['flux_fit'])*u.Jy
    flux_fit = flux_fit.to(outunit_S)*u.arcsec*u.arcsec## == aperture_sum!
    
    flux_err = np.array(measurements['flux_unc'])*u.Jy
    flux_err = flux_err.to(outunit_S)*u.arcsec*u.arcsec
    
    flux = divide_by_pixels_area(flux=flux_fit, bmaj=header['BMAJ']*u.deg,
                                 bmin=header['BMIN']*u.deg, pxsize=pixelsize_arcsec,
                                 aperture_correction_factor=1)
    e_flux = divide_by_pixels_area(flux=flux_err, bmaj=header['BMAJ']*u.deg,
                                   bmin=header['BMIN']*u.deg, pxsize=pixelsize_arcsec,
                                   aperture_correction_factor=1)


    measurements['Flux_ForcedPhoto'] = flux
    measurements['E_Flux_ForcedPhoto'] = e_flux
    
    outcols = ['Name', 'Flux_ForcedPhoto', 'E_Flux_ForcedPhoto']


    ###write to file
    outfilename = targetfile.split('.')[0] + '_EMU_forced_photometry.xml'
    measurements[outcols].write(outfilename, format='votable')
                                 
    return


def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description="files to perform EMU forced photometry on")
    parser.add_argument("targets", help="list of positions [deg] to perform forced photo")
    parser.add_argument("image", help="EMU image to perform photometry on")
    parser.add_argument("--config", action='store', type=str, default='config.txt',
                        help="config file")
        
    args = parser.parse_args()
    return args


########################################################################################
########################################################################################
###main code


if __name__ == "__main__":
    args = parse_args()
    ###add in columns and config file
    config = pd.read_table(args.config ,delim_whitespace=True).replace("'", "",regex=True)
    ra_col = config['value'][np.where(config['parameter']=='ra_col')[0][0]]
    dec_col = config['value'][np.where(config['parameter']=='dec_col')[0][0]]
    name_col = config['value'][np.where(config['parameter']=='name_col')[0][0]]
    
    psf_photo(targetfile=args.targets, imfile=args.image, acol=ra_col,
              dcol=dec_col, namecol=name_col)

