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
from astropy.stats import SigmaClip
from photutils import aperture
from photutils.psf import EPSFBuilder, EPSFFitter, extract_stars, DAOGroup, BasicPSFPhotometry
from photutils.background import StdBackgroundRMS, Background2D, MedianBackground
import argparse

#############################################################################
#############################################################################
###params

###time not important other than to show progress - background subtraction and psf flux can be slow for large images (on my laptop at any rate!)
########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################


###could make the below params a config file? probably fine without within emucat
###target file parameters
acol = 'ra'
dcol = 'dec'
namecol = 'designation'
tposunits = 'deg'
outname = 'forced_emu_photometry.xml'
outformat = 'votable'

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



class ForcedPhoto:
    'Aperture photometry functionality'
    
    def __init__(self, image_file, rms_file,
                 load_method=radio_image_as_2d,
                 outunit=u.mJy,
                 dxkey='CDELT1',
                 dykey='CDELT2',
                 xukey='CUNIT1',
                 yukey='CUNIT2'):
        
        ###load data
        self.hdu_im = load_method(image_file)
        self.hdu_rms = load_method(rms_file)
        
        ###get images, headers, beam and pixel parameters
        self.image = self.hdu_im[0].data
        self.rms = self.hdu_rms[0].data
        self.head = self.hdu_im[0].header
        self.wcs = WCS(self.head)
        self.beam = {'semi_major': (self.head['BMAJ']/2)*u.deg,
                     'semi_minor': (self.head['BMIN']/2)*u.deg,
                     'posangle': self.head['BPA']*u.deg}
        
        self.pxsizex = abs(self.head[dxkey])*u.Unit(self.head[xukey])
        self.pxsizey = abs(self.head[dykey])*u.Unit(self.head[yukey])
        self.pxarea = self.pxsizex*self.pxsizey
        self.pxunit = u.Unit(self.head['BUNIT'])
        self.outunit = outunit
        
        ###subtract background
        self.bkg = self.estimate_background(image=self.image)
        self.image = self.image - self.bkg

    
    
    def estimate_background(self, image, nsigma=3, box_size=(15, 15), filter_size=(3, 3)):
        'create background image array to subtract from main image'
        t0 = time.time()
        
        sigma_clip = SigmaClip(sigma=nsigma)
        bkg_estimator = MedianBackground()
        
        bkg = Background2D(image, box_size=box_size, filter_size=filter_size,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            
        bkg_array = bkg.background
                           
        print('bkg est time: ', np.round(time.time()-t0, 2), 's')
                           
        return(bkg_array)


    def aperture_sum_to_flux(self, apsum, a, b, pixel_area, aperture_correction=2):
        'convert aperture sum to flux assuming gausian beam'
    
        beam_area = 2*np.pi*a*b ###a and b are sigma of semi maj/min of beam
    
        px_per_beam = (beam_area/pixel_area) * (1/u.beam)
    
        flux = apsum / px_per_beam
    
        ###correct for incomplete aperture of Gaussian
        ##flux recovered by aperture based on sigma of Gaussian is half of total
        flux_cor = aperture_correction*flux
    
        return(flux_cor)


    def aperture_stat(self, image, image_wcs, apertures, stat=np.max):
        'obtain statistic pixel (e.g. max/min/median) value in aperture'
        'stat == np.max (default) is ≈equivalent to peak flux'
    
        ##create pixel apertures from sky
        pixaps = apertures.to_pixel(image_wcs)
    
        ###create image masks from pixel apertures - returns list of masks
        pixmasks = pixaps.to_mask()
    
        ###obtain maximum pixel value in cutout defined by mask
        statpx = np.array([stat(i.cutout(image)) for i in pixmasks])
    
        return(statpx)


    def aperture_measurements(self, targets, apmethod='exact',
                              nsubpx=3):
        'get flux, error, and max pixel value in aperture'
        t0 = time.time()

        ###get sigma from fwhm for semi major/minor
        beam_a = self.beam['semi_major']/(2*np.sqrt(2*np.log(2)))
        beam_b = self.beam['semi_minor']/(2*np.sqrt(2*np.log(2)))

        ###create apertures
        target_apertures = aperture.SkyEllipticalAperture(positions=targets,
                                                          a=beam_a,
                                                          b=beam_b,
                                                          theta=self.beam['posangle'])
                                                          
        ###get aperture sum
        apphoto = aperture.aperture_photometry(data=self.image,
                                               error=self.rms,
                                               apertures=target_apertures,
                                               wcs=self.wcs,
                                               method=apmethod,
                                               subpixels=nsubpx)
        
        ###aperture sum units (px units / beam):
        unitcols = ['aperture_sum']
        for col in unitcols:
            apphoto[col].unit = self.pxunit

        ###aperture flux = divide by area (px/beam) and (assuming Guassian beam) x2.
        apphoto['ap_flux'] = self.aperture_sum_to_flux(apphoto['aperture_sum'],
                                                       a=beam_a,
                                                       b=beam_b,
                                                       pixel_area=self.pxarea)

        ###flux error -- need RMS image (return flux of ap in rms im and median rms in ap)
        apphoto['aperture_sum_err'].unit = apphoto['aperture_sum'].unit
        apphoto['median_ap_rms'] = self.aperture_stat(image=self.rms,
                                                      image_wcs=self.wcs,
                                                      apertures=target_apertures,
                                                      stat=np.nanmedian)*self.pxunit
        apphoto['ap_flux_err'] = self.aperture_sum_to_flux(apphoto['aperture_sum_err'],
                                                           a=beam_a,
                                                           b=beam_b,
                                                           pixel_area=self.pxarea)

        ###brightest pixel value
        apphoto['ap_max'] = self.aperture_stat(image=self.image,
                                               image_wcs=self.wcs,
                                               apertures=target_apertures,
                                               stat=np.nanmax)*self.pxunit

        ###convert to output units
        newunits = {'aperture_sum': self.outunit/u.beam,
                    'ap_flux': self.outunit,
                    'ap_flux_err': self.outunit,
                    'ap_max': self.outunit/u.beam,
                    'median_ap_rms': self.outunit/u.beam}

        for key in list(newunits.keys()):
            apphoto[key] = apphoto[key].to(newunits[key])
    
        aperture_results = {'aperture_flux': apphoto['ap_flux'],
                            'aperture_flux_err': apphoto['ap_flux_err'],
                            'aperture_max': apphoto['ap_max'],
                            'aperture_rms_median': apphoto['median_ap_rms']}
                                
        print('apphoto time: ', np.round(time.time()-t0, 2), 's')
    
        return(aperture_results)
                            

    def model_gaussian_from_beam(self):
        'models 2d gaussian in px from beam info'
        ##beam info
        bmaj = (2*self.beam['semi_major'])/(2*np.sqrt(2*np.log(2)))
        bmin = (2*self.beam['semi_minor'])/(2*np.sqrt(2*np.log(2)))

        deg_per_px = self.pxsizex/u.pixel #(abs(image_hdu[0].header[dxkey])*(u.deg/u.pixel))
        dx, dy = bmaj/deg_per_px, bmin/deg_per_px
    
        ###convert bpa to radians and add pi/2 to get to same ref point as image
        da = self.beam['posangle'].to(u.rad) + (np.pi/2)*u.rad
    
        ##obtain Gaussian Kernal -- linear interp seems better than center/oversample
        gauss = Gaussian2DKernel(dx.value, dy.value, da.value, mode='linear_interp')
    
        return(gauss)


    def model_psf_from_beam(self, osamp=4, iters=10):
        'model psf from the beam geometry'
        ##define Gaussian
        beam_gauss = np.array(self.model_gaussian_from_beam())
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


    def psf_measurements(self, targets,
                         osamp=4, iters=3,
                         n_beams=3):
        'obtain flux by modelling beam psf as gaussian'
        t0 = time.time()

        ###get sigma from fwhm for major/minor
        beam_a = (2*self.beam['semi_major'])/(2*np.sqrt(2*np.log(2)))
        beam_b = (2*self.beam['semi_minor'])/(2*np.sqrt(2*np.log(2)))


        ##2) convert sky coords to pixel positions
        pixel_positions = self.wcs.world_to_pixel(targets)
    
        tpos_px = Table()
        tpos_px['x_0'] = pixel_positions[0]
        tpos_px['y_0'] = pixel_positions[1]
        ###would using max pixel value in aperture as initial flux guess improve results?
    
        ##3) model psf
        epsf = self.model_psf_from_beam(osamp=osamp, iters=iters)

        ##4) set up pixel apertures and perform photometry
        ###use N beam widths (in pxls) for grouping and aperture radius for init_guesses
        beam_width_pixels = (beam_a/self.pxsizex)
        pxap = int(np.ceil(n_beams*beam_width_pixels))
        if pxap%2 == 0:
            fitsize = pxap+1
        else:
            fitsize = pxap

        ###fix psf positions to input
        epsf.x_0.fixed = True
        epsf.y_0.fixed = True

        photometry = BasicPSFPhotometry(group_maker=DAOGroup(pxap), bkg_estimator=None,
                                        psf_model=epsf, fitshape=fitsize,
                                        fitter=LevMarLSQFitter(),
                                        aperture_radius=pxap)

        measurements = photometry(self.image, tpos_px)

        ##need to convert flux_fit and uncertainty to real fluxes in correct units
        ffit = np.array(measurements['flux_fit'])*u.Jy
        ffit = ffit.to(self.outunit)*(1/u.beam)## == aperture_sum!
    
        flux = self.aperture_sum_to_flux(apsum=ffit,
                                         a=beam_a,
                                         b=beam_b,
                                         pixel_area=self.pxarea,
                                         aperture_correction=1)
        
        ffit_err = np.array(measurements['flux_unc'])*u.Jy
        ffit_err = ffit_err.to(self.outunit)*(1/u.beam)## == aperture_sum!
                                
        flux_err = self.aperture_sum_to_flux(apsum=ffit_err,
                                             a=beam_a,
                                             b=beam_b,
                                             pixel_area=self.pxarea,
                                             aperture_correction=1)
                                
        psf_results = {'psf_flux': flux, 'psf_flux_err': flux_err}
                                
        print('psf time: ', np.round(time.time()-t0, 2), 's')

        return(psf_results)


    def photoTable(self, targets, include_psf=True):
        'output astropy table with results, option to not do psf photo as this is slower'
        
        aphot = self.aperture_measurements(targets=targets)
        if include_psf == True:
            pphot = self.psf_measurements(targets=targets)
        else:
            pphot = {}

        return(Table({**aphot, **pphot}))


def parse_args():
    "parse input args, i.e. target list, image, and rms file names"
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("targets", help="file with list of sky coords to perform forced photometry on")
    parser.add_argument("image", help="image file")
    parser.add_argument("noise", help="noise file")
                        
    args = parser.parse_args()


    return args


#############################################################################
#############################################################################
###main code -- make command line


if __name__ == '__main__':
    ##parse command line arguments
    args = parse_args()
    
    ###load data, grab target names and positions
    targets = Table.read(args.targets)
    targetids = targets[namecol]
    targetpos = SkyCoord(ra=targets[acol], dec=targets[dcol], unit=tposunits)
    
    ##perform photometry
    photo = ForcedPhoto(image_file=args.image, rms_file=args.noise, load_method=fits.open)
    phototab = photo.photoTable(targets=targetpos, include_psf=True)
    
    ##finalise table and write to file
    phototab.add_column(col=targetids, index=0) ##adds ID column
    phototab.write(outname, format=outformat)





