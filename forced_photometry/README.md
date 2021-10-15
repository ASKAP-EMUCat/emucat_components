**emu_forced_photo.py readme**

command line code for performing forced photometry on EMU images for a given set of positions (e.g. WISE sources). requires file of target positions, EMU image and EMU rms image. Run as:

*python3 emu_forced_photo.py target_list_file emu_image_file emu_rms_image_file*

the directory *test_data* contains a small list of AllWISE targets and 1 square degree image and rms files on which to test the code. To use these use the following command:

*python3 emu_forced_photo.py test_data/sb9442_cen_1degsq_AWISE_AGN.xml test_data/sb9442_central_1sqdeg_t0.fits test_data/sb9442_central_1sqdeg_t0rms.fits*

The optional command --model_psf can also be used to disable the slower photutils psf modelling used to obtain flux measurements. By default this is set to True, but if set to False the code will run substantially faster and the output file will not contain the *psf_flux* and *psf_flux_err* columns.
The file of targets should include columns for the name/id, RA and Dec of the target. The default assumption is that these columns are called 'designation', 'ra' and 'dec', and that 'ra' and 'dec' are in decimal degrees -- these parameters can be changed at the top of the code if necessary. The output is a votable format file containing the names/id of the sources and basic photometric measurements.


**Output Data**

parameter | description
----------|------------
name | name/id of the object (note this column name will be the same as for the input target file)
aperture_flux | flux within an aperture based on beam geometry [mJy]
aperture_flux_err | 1 sigma error on flux within an aperture based on beam geometry [mJy]
aperture_max | brightest pixel value within the aperture [mJy / beam]
aperture_rms_median | median rms value within the aperture [mJy / beam]
psf_flux | point source flux using psf modelled on beam geometry [mJy]
psf_flux_err | 1 sigma error on point source flux using psf modelled on beam geometry [mJy]


**Code Dependencies**
This code was developed using python v3.7.2 along with the following packages (version used in development):
* argparse (1.1)
* astropy (4.3.1)
* matplotlib (2.2.2)
* numpy (1.17.4)
* photutils (0.7)



