####code to make aperture photometry measurements on radio image
####test on cutouts of a few targets
####make so works with list of positions -- len(list)==1 for single target
#### IN: SkyCoord(s)
#### OUT: aperture sum, guassian flux based on sum, brightest pixel in aperture and modelled guassian flux (include units)
#### define aperture based on beam size


import numpy as np, matplotlib.pyplot as plt, os
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import NDData
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.convolution import Gaussian2DKernel
from photutils import aperture
from photutils.psf import EPSFBuilder, EPSFFitter, extract_stars, DAOGroup, BasicPSFPhotometry
from photutils.background import StdBackgroundRMS
import argparse


#############################################################################
#############################################################################
###params

###test data
tfile = '../data/J212057-550722_emu_cutout.fits'
tname = 'J212057-550722'
tpos = SkyCoord(ra=320.241375, dec=-55.122992, unit='deg')
ts_peak = 1.046*(u.mJy/u.beam)
ts_int = 1.015*u.beam

#tfile = '../data/J212906-561822_emu_cutout.fits'
#tname = 'J212906-561822'
#tpos = SkyCoord(ra=322.277273, dec=-56.306186, unit='deg')
#ts_peak = 0.344*(u.mJy/u.beam)
#ts_int = 0.336*u.beam

#tfile = '../data/J212748-564806_emu_cutout.fits'
#tname = 'J212748-564806'
#tpos = SkyCoord(ra=321.952077, dec=-56.801809, unit='deg')
#ts_peak = 5.144*(u.mJy/u.beam)
#ts_int = 5.026*u.beam

tplist = SkyCoord(ra=[tpos.ra], dec=[tpos.dec])

##header keys for pix size (not always 'CDELT1')
dxkey = 'CDELT1'
xukey = 'CUNIT1'
dykey = 'CDELT2'
yukey = 'CUNIT2'

#############################################################################
#############################################################################
###functionality


def radio_image_as_2d(file):
    'takes 4D image file and outputs 2D HDU for ease of working with astropy on flat image'
    
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
        if 'PC3' in key or 'PC4' in key or '_3' in key or '_4' in key:
            badkeys.append(key)
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


def aperture_sum_to_flux(apsum, a, b, pixel_area, aperture_correction=2):
    'convert aperture sum to flux assuming gausian beam'
    
    beam_area = 2*np.pi*a*b ###a and b are sigma of semi maj/min of beam
    
    px_per_beam = (beam_area/pixel_area) * (1/u.beam)
    
    flux = apsum / px_per_beam
    
    ###correct for incomplete aperture of Gaussian
    ##flux recovered by aperture based on sigma of Gaussian is half of total
    flux_cor = aperture_correction*flux
    
    return(flux_cor)


def aperture_stat(image, header, image_wcs, apertures, stat=np.max):
    'obtain statistic pixel (e.g. max/min/median) value in aperture'
    'stat == np.max (default) is equivalent to peak flux'

    imwcs = WCS(header)
    
    ##create pixel apertures from sky
    pixaps = apertures.to_pixel(WCS(header))
    
    ###create image masks from pixel apertures - returns list of masks
    pixmasks = pixaps.to_mask()
    
    ###obtain maximum pixel value in cutout defined by mask
    statpx = np.array([stat(i.cutout(image)) for i in pixmasks])

    return(statpx)



def aperture_flux(targets, image_hdu, outunit=u.mJy, apmethod='exact',
                  nsubpx=3):
    'get flux, error, and max pixel value in aperture'
    ####errors require rms image, need to get hold of for testing
    
    image, header = image_hdu[0].data, image_hdu[0].header
    imwcs = WCS(header)
    
    ##pixel size
    pxsizex = abs(header[dxkey])*u.Unit(header[xukey])
    pxsizey = abs(header[dykey])*u.Unit(header[yukey])
    pixarea = pxsizex*pxsizey
    pxunit = u.Unit(header['BUNIT'])
    
    ##beam geometry (semi major/minor are FWHM)
    beam = {'major': (header['BMAJ']/2)*u.deg,
            'minor': (header['BMIN']/2)*u.deg,
            'posangle': header['BPA']*u.deg}

    ###get sigma from fwhm for semi major/minor
    beam_a = beam['major']/(2*np.sqrt(2*np.log(2)))
    beam_b = beam['minor']/(2*np.sqrt(2*np.log(2)))

    ###create apertures
    target_apertures = aperture.SkyEllipticalAperture(positions=targets,
                                                      a=beam_a,
                                                      b=beam_b,
                                                      theta=beam['posangle'])
    
    ###get aperture sum
    apphoto = aperture.aperture_photometry(data=image,
                                           apertures=target_apertures,
                                           wcs=imwcs,
                                           method=apmethod,
                                           subpixels=nsubpx)
        
    ###aperture sum units (px units / beam):
    unitcols = ['aperture_sum']
    for col in unitcols:
            apphoto[col].unit = pxunit

    ###aperture flux = divide by area (px/beam) and (assuming Guassian beam) x2.
    apphoto['ap_flux'] = aperture_sum_to_flux(apphoto['aperture_sum'],
                                              a=beam_a,
                                              b=beam_b,
                                              pixel_area=pixarea)

    ###flux error -- need RMS image (take median of rms?)

    ###brightest pixel value
    apphoto['ap_max'] = aperture_stat(image=image,
                                      header=header,
                                      image_wcs=imwcs,
                                      apertures=target_apertures,
                                      stat=np.nanmax)*pxunit

    ###convert to output units
    newunits = {'aperture_sum': outunit/u.beam,
                'ap_flux': outunit,
                'ap_max': outunit/u.beam}

    for key in list(newunits.keys()):
        apphoto[key] = apphoto[key].to(newunits[key])

    aperture_results = {'aperture_flux': apphoto['ap_flux'],
                        'aperture_max': apphoto['ap_max']}

    return(aperture_results)




####model psf as gaussian and use psf photometry

def model_gaussian_from_beam(image_hdu):
    'models 2d gaussian in px from beam info'
    ##beam info
    bmin = image_hdu[0].header['BMIN']*u.deg
    bmaj = image_hdu[0].header['BMAJ']*u.deg
    bpa = image_hdu[0].header['BPA']*u.deg
    
    ###convert to sigma from FWHM
    bmaj = bmaj/(2*np.sqrt(2*np.log(2)))
    bmin = bmin/(2*np.sqrt(2*np.log(2)))
    
    deg_per_px = (abs(image_hdu[0].header[dxkey])*(u.deg/u.pixel))
    dx, dy = bmaj/deg_per_px, bmin/deg_per_px

    ###convert bpa to radians and add pi/2 to get to same ref point as image
    da = bpa.to(u.rad) + (np.pi/2)*u.rad
    
    ##obtain Gaussian Kernal -- linear interp seems better than center/oversample
    gauss = Gaussian2DKernel(dx.value, dy.value, da.value, mode='linear_interp')
    
    return(gauss)


def model_psf_from_beam(image_hdu, osamp=4, iters=10):
    'model psf from the beam geometry'
    ##define Gaussian
    beam_gauss = np.array(model_gaussian_from_beam(image_hdu))
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


def psf_flux(targets, image_hdu,
             outunit=u.mJy,
             osamp=4, iters=3,
             bkg_sub=None, n_beams=3):
    
    image, header = image_hdu[0].data, image_hdu[0].header
    
    ##pixel size
    pxsizex = abs(header[dxkey])*u.Unit(header[xukey])
    pxsizey = abs(header[dykey])*u.Unit(header[yukey])
    pixarea = pxsizex*pxsizey
    pxunit = u.Unit(header['BUNIT'])
    
    ##beam geometry (semi major/minor are FWHM)
    beam = {'major': (header['BMAJ'])*u.deg,
            'minor': (header['BMIN'])*u.deg,
            'angle': header['BPA']*u.deg}
    
    ###get sigma from fwhm for major/minor
    beam_a = beam['major']/(2*np.sqrt(2*np.log(2)))
    beam_b = beam['minor']/(2*np.sqrt(2*np.log(2)))


    ##2) convert sky coords to pixel positions
    wcs = WCS(header)
    pixel_positions = wcs.world_to_pixel(targets)

    
    tpos_px = Table()
    tpos_px['x_0'] = pixel_positions[0]
    tpos_px['y_0'] = pixel_positions[1]
    ###would using max pixel value in aperture as initial flux guess improve results?

    ##3) model psf
    epsf = model_psf_from_beam(image_hdu=image_hdu, osamp=osamp, iters=iters)
    
    ##4) perform PSF photometry
    bkg = bkg_sub #StdBackgroundRMS() or None which produces better results
    
    ###use N beam widths (in pxls) for grouping and aperture radius for init_guesses
    beam_width_pixels = (beam_a/pxsizex)
    pxap = int(np.ceil(n_beams*beam_width_pixels))
#    pxap = 10
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

    measurements = photometry(image, tpos_px)
    
    ##need to convert flux_fit and uncertainty to real fluxes in correct units
    ffit = np.array(measurements['flux_fit'])*u.Jy
    ffit = ffit.to(outunit)*(1/u.beam)## == aperture_sum!

    flux = aperture_sum_to_flux(apsum=ffit, a=beam_a, b=beam_b,
                                pixel_area=pixarea, aperture_correction=1)

    ffit_err = np.array(measurements['flux_unc'])*u.Jy
    ffit_err = ffit_err.to(outunit)*(1/u.beam)## == aperture_sum!

    flux_err = aperture_sum_to_flux(apsum=ffit_err, a=beam_a, b=beam_b,
                                    pixel_area=pixarea, aperture_correction=1)

    psf_results = {'psf_flux': flux, 'psf_flux_err': flux_err}

    return(psf_results)



def get_photometry(target_list, image_file, load_image_data=radio_image_as_2d,
                   output_units=u.mJy, apmethod='subpixel', nsubpx=3,
                   namecol='name', acol='ra', dcol='dec', include_psf_photo=True):
    'get aperture and psf photometry and output data'
    
    ###make sure columns are in data -- dont want this to break AFTER the hard work is done
    needcols = [namecol, acol, dcol]
    for col in needcols:
        if col not in target_list.colnames:
            print(col + ' not in target list colnames, aborting')
            return
    
    print('loading image')
            
    
    ###load image data
    image_hdu = load_image_data(image_file)
    
    ###create skycoords from target list -- requires units in input table
    targets = SkyCoord(ra=target_list[acol], dec=target_list[dcol])

    print('getting aperture photometry')
    ###get photometric measurements
    apphot = aperture_flux(targets=targets, image_hdu=image_hdu,
                           outunit=output_units, apmethod=apmethod,
                           nsubpx=nsubpx)

    if include_psf_photo==True:
        print('getting psf photometry -- this may take a while')
        psfphot = psf_flux(targets=targets, image_hdu=image_hdu,
                           outunit=output_units)
    else:
        psfphot = {}

    print('measurements made, writing table')
    ###combine with target name in dict
    photometry = {**{namecol: target_list[namecol]}, **apphot, **psfphot}
    
    ###write astropy table
    phottab = Table(photometry)
    
    return(phottab)





#############################################################################
#############################################################################
###test
#####2D Gaus
#### V = 2pi*A*sx*sy
#### A = V / 2pi*sx*sy
#### sx; sy = a/2.355, b/2.355

outunit=u.mJy ###units to output flux in.



targets = tplist

hdu = fits.open(tfile)


image, header = hdu[0].data, hdu[0].header
imwcs = WCS(header)

##pixel size
pxsizex = abs(header[dxkey])*u.Unit(header[xukey])
pxsizey = abs(header[dykey])*u.Unit(header[yukey])
pixarea = pxsizex*pxsizey
pxunit = u.Unit(header['BUNIT'])

##beam geometry (semi major/minor are FWHM)
beam = {'major': (header['BMAJ']/2)*u.deg,
        'minor': (header['BMIN']/2)*u.deg,
        'posangle': header['BPA']*u.deg}

###get sigma from fwhm for semi major/minor
beam_a = beam['major']/(2*np.sqrt(2*np.log(2)))
beam_b = beam['minor']/(2*np.sqrt(2*np.log(2)))


#imload = fits.open
#imfile = tfile
imload = radio_image_as_2d
imfile = ('../../../../../survey_data/EMU/Pilot_survey/images/t0/'
          + 'image.i.SB9442.cont.taylor.0.restored.fits')
targettab = Table({'name': [tname], 'ra': tplist.ra, 'dec': tplist.dec})

targettab = Table.read('../data/test_200_point_fromSB9442.fits')
targettab = Table.read('../data/SB9442_central_point_sources.fits')


tcat = SkyCoord(ra=targettab['ra_deg_cont'], dec=targettab['dec_deg_cont'])


hd2 = radio_image_as_2d(imfile)
#test = aperture_flux(tcat, hd2)

target_photo = get_photometry(target_list=targettab, image_file=imfile,
                              load_image_data=imload, output_units=u.mJy,
                              namecol='component_name', acol='ra_deg_cont',
                              dcol='dec_deg_cont', include_psf_photo=False)

#test = get_photometry(target_list=targettab, image_file=imfile, load_image_data=imload,
#                      output_units=outunit)
#print(test)

#tflux = test['aperture_flux'][0]
#tflux = test['psf_flux'][0]

#intrat = (tflux/ts_int).value
#peakrat = (tflux/ts_peak).value
#
#print(' ')
#print('total: ', np.round(100*intrat, 2), '%')
#print('peak: ', np.round(100*peakrat, 2), '%')
#print(' ')



#test = aperture_flux(targets=tplist, image_hdu=hdu, outunit=u.mJy)
#tint = test['ap_flux'][0]
#tpeak = test['ap_max'][0]
#
#
#intrat = (tint/ts_int).value
#peakrat = (tint/ts_peak).value
#
#print('total: ', np.round(100*intrat, 2))
#print('peak: ', np.round(100*peakrat, 2))


#############################################################################
#############################################################################
###main



