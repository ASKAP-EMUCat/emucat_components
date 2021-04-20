**emu_forced_photo.py readme**

command line code for performing forced photometry on EMU images for a given set of positions (e.g. WISE sources). requires file of target positions, EMU image and EMU rms image. Run as:

*python3 emu_forced_photo.py target_list_file emu_image_file emu_rms_image_file*

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





