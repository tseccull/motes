#! /home/tseccull/anaconda2/envs/anaconda3/bin/python

###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import astropy.io.fits as fits
import datetime
import glob
import motes.common as common
import motes.harvester as harvester
import motes.startup as startup
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

###############################################################################
# FUNCTIONS //////////////////////////////////////////////////////////////////#
###############################################################################
# ////////////////////////////////////////////////////////////////////////////#
# HANDLE COSMIC RAYS AND BAD PIXELS

def cr_handling(pars, axdict, lext, hext, fdict, hpars):
    sys.stdout.write(' >>> Separating spectrum and background regions for cosmic ray detection.\n')

    # DIAGNOSTICS - Plot limits defined for CR masking process.
    if pars['-DIAG_PLOT_CRMASK']:
        daxis = np.linspace(axdict['wavstart'], axdict['wavstart'] + len(axdict['waxis']), num=len(axdict['waxis']))
        loline = [lext + axdict['imgstart']] * len(axdict['waxis'])
        hiline = [hext + axdict['imgstart']] * len(axdict['waxis'])
        dlines = np.array([daxis, loline, daxis, hiline])
        common.show_img(fdict['data'], axdict, hpars, dlines, '2D Spectrum Mask Regions - Pre Masking')


    sys.stdout.write(' >>> Masking cosmic rays via sigma-clip:\n')
    sys.stdout.write('     Sigma multiplier = ' + str(pars['-CR_CLIP_STD_MULT']) + '\n')
    sys.stdout.write('     Number of sigma-clip iterations = ' + str(pars['-CR_CLIP_ITER']) + '\n')
    sys.stdout.write(' >>> Masking bad pixels.\n')
    fdict['cmask'] = np.ones(np.shape(fdict['data']))
    for it in range(int(pars['-CR_CLIP_ITER'])):
        crmasque = common.mask_cosmic_rays(fdict['data'], lext, hext, multiplier=pars['-CR_CLIP_STD_MULT'])
        fdict['qual'] *= crmasque
        fdict['data'] *= fdict['qual']
        fdict['errs'] *= fdict['qual']
        fdict['cmask'] *= crmasque
    sys.stdout.write(' >>> Cosmic rays and bad pixels masked.\n')

    # DIAGNOSTICS - Plot 2D spectrum with CRs masked.
    if pars['-DIAG_PLOT_CRMASK']:
        common.show_img(fdict['data'], axdict, hpars, dlines, '2D Spectrum Mask Regions - Post Masking')

    # If this is a GMOS spectrum restore zeros where NaNs were put in in the chip gaps.
    # Has no effect on spectra from other instruments.
    nanlocs = np.where(np.isnan(fdict['data'])==True)
    fdict['data'][nanlocs] == 0
    fdict['errs'][nanlocs] == 0
    fdict['qual'][nanlocs] == 0

    # Determine the location of boundaries between dispersion bins in the 2D spectrum that will be fitted with Moffat
    # functions to precisely localise the spectrum to allow sky subtraction.
    # If cosmic rays have been masked, and the user wants them to be replaced with appropriate flux values the
    # process is performed while determining the location of the dispersion bins. See the docs or comments in the
    # common.get_bins function for more information.
    # Replacement is mostly, but not entirely effective. It may cause particular problems when used on fainter
    # spectra, or those with very bright sky lines. Some fiddling with the CR sigma-clipping multipliers may be
    # required in the sky subtraction and CR masking functions. Use with care.
    
    # Determine the boundaries of the localisation bins on the dispersion axis
    binpars, fdict = common.get_bins(fdict, int(np.floor(lext)), int(np.ceil(hext)), axdict['dispaxislen'], pars, sky=True)
    
    if pars['-REPLACE_CRBP']:
        sys.stdout.write(' >>> Cosmic rays and bad pixels replaced.\n')

        # DIAGNOSTICS - Show 2D spectrum with bad pixels and CRs replaced.
        if pars['-DIAG_PLOT_CRREPLACE']:
            common.show_img(fdict['data'], axdict, hpars, dlines, '2D Spectrum with Cosmic Rays and Bad Pixels Replaced')

    # If replacement of cosmic rays is not requested, just get the spectrum localisation bins from the 2D spectrum
    # that has had it's CRs masked.
    else:
        sys.stdout.write(' >>> Cosmic ray and bad pixel replacement not performed.\n')

    common.get_bins_output(binpars, pars, lext, hext, fdict['data'], hpars, axdict)

    return binpars, fdict


# ////////////////////////////////////////////////////////////////////////////#
# RUN MOTES
def motes():
    # RUN STARTUP FUNCTIONS
    params = startup.read_parfile()    # Import parameters from file to dict
    intreg = startup.read_regions()    # Search for, and read in, reg.txt

    # OPEN AND PROCESS EACH 2D SPECTRUM IN TURN
    for i, file_2D in enumerate(sorted(glob.glob('*.fits'))):

        sys.stdout.write(('/' * (70-len(file_2D[:70]))) + ' ' + file_2D[:70] + '\n')
        sys.stdout.write(' >>> Beginning MOTES Processing\n')

        # GATHER HEADER INFO AND DATA FROM THE 2D IMAGE FILE
        sys.stdout.write(' >>> Gathering image frames and header data from input file.\n')
        headparams, framedict, axesdict, imghead = harvester.data_harvest(i, file_2D, intreg, params)
        framedict['ogdata'] = framedict['data']
        framedict['ogerrs'] = framedict['errs']
        sys.stdout.write(' >>> Gathering image frames and header data from input file completed.\n')

        # Perform initial least-squares Moffat fitting over the entire 2D image collapsed along the dispersion axis with a median.
        sys.stdout.write(' >>> Fitting Moffat profile to median spatial profile of entire spectrum. ')
        sys.stdout.flush()
        datadispcollapse = np.nanmedian(framedict['data'], axis=1)
        datascale = 10**np.abs(np.floor(np.log10(np.median(datadispcollapse))))
        moffparams = common.moffat_least_squares(axesdict['saxis'], datadispcollapse*datascale, headparams['seeing'], headparams['pixresolution'])
        headparams['seeing'] = 2*moffparams[2]*np.sqrt((2**(1/moffparams[3]))-1)

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle.
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled down again after the fitting.
        moffparams[0] /= datascale
        moffparams[4] /= datascale
        moffparams[5] /= datascale

        sys.stdout.write('DONE.\n')
        sys.stdout.write(' >>> FWHM of median spatial profile is ' + str(round(headparams['seeing'], 2)) + ' Pixels, or ' + str(round(headparams['seeing']*headparams['pixresolution'], 2)) + '"\n')

        # Use the parameters of the Moffat profile fitted to the median spatial profile of the entire spectrum to determine
        # spatial limits that are used to bound the region of the spectrum used by the common.get_bins function to
        # to bin the 2D spectrum while taking account of its S/N. These spatial limits may also be used during the
        # process of masking, removing, and replacing cosmic rays and bad pixels.
        lowext, highext, fwhm, cent = common.extraction_limits(moffparams, params['-CR_FWHM_MULTIPLIER'], axesdict)
        sys.stdout.write(' >>> Spectrum localised to aperture in range of spatial pixel rows ' + str(int(lowext+axesdict['imgstart'])) + '-' + str(int(highext+axesdict['imgstart'])) + '\n')

        # DIAGNOSTICS -  Plot fitted Moffat profile over collapsed 2D spectrum and print the parameters of the fitted Moffat profile.
        if params['-DIAG_PLOT_COLLAPSED_2D_SPEC']:
            common.printmoffparams(moffparams, axesdict['imgstart'], datascale)
            common.plot_fitted_spatial_profile(axesdict['saxis'], datadispcollapse, axesdict['hrsaxis'], moffparams, axesdict['imgstart'], headparams)


        # Mask and remove cosmic rays and bad pixels and, if selected, replace them in the image.
        if params['-MASK_CR']:
            binparams, framedict = cr_handling(params, axesdict, lowext, highext, framedict, headparams)

        # No cosmic ray masking requested, just get localisation bins
        else:
            sys.stdout.write(' >>> Cosmic ray masking deactivated and not performed.\n')

            # If this is a GMOS spectrum restore zeros where NaNs were put in in the chip gaps.
            # Has no effect on spectra from other instruments.
            wherenans = np.where(np.isnan(framedict['data'])==True)
            framedict['data'][wherenans] == 0
            framedict['qual'][wherenans] == 0

            framedict['cmask'] = np.ones(np.shape(framedict['data']))
            binparams, framedict = common.get_bins(framedict, int(np.floor(lowext)), int(np.ceil(highext)), axesdict['dispaxislen'], params, sky=True)
                                                                              
            common.get_bins_output(binparams, params, axesdict['saxis'], lowext, highext, axesdict['wavstart'], axesdict['imgstart'], framedict['data'], axesdict['waxis'], headparams)
            sys.stdout.write(' >>> Bad pixels replaced.\n')


        #SUBTRACT THE SKY SPECTRUM IF REQUESTED BY THE USER.
        if params['-SUBTRACT_SKY']:
            framedict, skybinpars, skyextlims = skyloc(framedict, axesdict, datascale, headparams, binparams, params)

        binparams, framedict = common.get_bins(framedict, int(np.floor(lowext)), int(np.ceil(highext)), axesdict['dispaxislen'], params)

        common.get_bins_output(binparams, params, lowext, highext, framedict['data'], headparams, axesdict)

        # For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
        # moffat profile to the median data and then use the parameters of the fitted Moffat function to localise the 2D
        # spectrum.
        sys.stdout.write(' >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n')

        extbin = []
        extractionlimits = []
        for bin in binparams:
			# Take the median spatial profile of the dispersion 
		    # bin, and leave out pixel columns in the chip gaps if this is a GMOS spectrum.
            binimg = framedict['data'][:, bin[0]:bin[1]]
            chipgap = np.where(np.median(binimg, axis=0)!=1)
            bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)

            # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median spatial profile and
            # return its parameters.
            binmoffparams = common.moffat_least_squares(axesdict['saxis'], bindata*datascale, headparams['seeing'], headparams['pixresolution'])

            binmoffparams[0] /= datascale
            binmoffparams[4] /= datascale
            binmoffparams[5] /= datascale

            # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
            # previously fitted to it.
            LowExt, HighExt, fwhm, centre = common.extraction_limits(binmoffparams, params['-FWHM_MULTIPLIER'], axesdict, highres=True)

            extractionlimits.append([(bin[0] + bin[1]) * 0.5, LowExt, HighExt, centre])

            # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
            # locations so they can be saved as metadata along with the extracted spectrum.
            binmoffparams.append(bin[0]+axesdict['wavstart'])
            binmoffparams.append(bin[1]+axesdict['wavstart'])
            extbin.append(binmoffparams)

            # DIAGNOSTICS - Plot computed moffat profile over data for each bin
            if params['-DIAG_PLOT_MOFFAT']:
                common.plot_fitted_spatial_profile(axesdict['saxis'], bindata, axesdict['hrsaxis'], binmoffparams, axesdict['imgstart'], headparams)

        binpars = np.array(extbin)
        sys.stdout.write('     Fitting complete.\n')

        sys.stdout.write(' >>> Drawing extraction aperture limits. ')
        sys.stdout.flush()
        extractionlimits = np.array(extractionlimits).T

        # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
        if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
            drawlines = [extractionlimits[0] + axesdict['wavstart'], (extractionlimits[1] * 0.02) + axesdict['imgstart'],
                         extractionlimits[0] + axesdict['wavstart'], (extractionlimits[2] * 0.02) + axesdict['imgstart']]
            
            common.show_img(framedict['data'], axesdict, headparams, drawlines, '2D Spectrum Overplotted with Extraction Limits')

        # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
        # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.
        finalextractionlims = common.interpolate_extraction_lims(extractionlimits, axesdict['dispaxislen'], params['-INTERP_KIND'])
        sys.stdout.write('DONE.\n')

        # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis.
        if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
            drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (finalextractionlims[0]*0.02) + axesdict['imgstart'],
                         np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (finalextractionlims[1]*0.02) + axesdict['imgstart']]
            
            common.show_img(framedict['data'], axesdict, headparams, drawlines, '2D Spectrum Overplotted with Full Extraction Limits')

        # Extract the spectrum from a supersampled version of the 2D image using the extraction limits.
        sys.stdout.write(' >>> Extracting 1D spectrum. ')
        sys.stdout.flush()
        supersampdata2D = (np.repeat(framedict['data'], 50, axis=0) * 0.02).T
        supersamperrs2D = (np.repeat(framedict['errs'], 50, axis=0) * 0.02).T
        data1D = []
        errs1D = []

        for k in range(axesdict['dispaxislen']):
            data1D.append(np.sum(supersampdata2D[k][int(finalextractionlims[0][k]):int(finalextractionlims[1][k])]))
            errs1D.append(np.sqrt(np.sum(np.power(supersamperrs2D[k][int(finalextractionlims[0][k]):int(finalextractionlims[1][k])], 2))))
        data1D = np.array(data1D)
        errs1D = np.array(errs1D)
        finalextractionlims = np.array(finalextractionlims) * 0.02
        sys.stdout.write('DONE.\n')

        # DIAGNOSTICS - Plot extracted spectrum.
        if params['-PLOT_EXTRACTED_SPECTRUM']:
            # DIAGNOSTICS, EXTRACTED SPECTRUM
            plt.plot(axesdict['waxis'], data1D, color='k')
            plt.grid(alpha=0.5, linestyle='dotted')
            plt.title('Extracted 1D Spectrum')
            plt.ylabel('Flux, ' + headparams['fluxunit'])
            plt.xlabel('Wavelength, ' + headparams['wavunit'])
            plt.show()

        # Save the extracted spectrum to a new .fits file.
        if params['-SAVE']:
            sys.stdout.write(' >>> Saving 1D spectrum and metadata.\n')
            sys.stdout.flush()
            save_fits(axesdict, headparams, data1D, errs1D, imghead, params, file_2D, moffparams, framedict, binpars, finalextractionlims, skybinpars, skyextlims)

        sys.stdout.write(' >>> Extraction of ' + file_2D + ' completed.\n')
    sys.stdout.write(' >>> MOTES Processing Complete.\n\n')
    
    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function saves the extracted spectrum and intermediate products in a single FITS file.
def save_fits(axdict, hparams, flux, errs, head, pars, filename, moffpars, fdict, bpars, extractionlims, sbpars, skyextractionlims):

    head['MOTES']='######## Extracted 1D Spectrum Metadata ########'
    head.add_blank('', before='MOTES')
    head['HIERARCH UTC EXT DATE'] = datetime.datetime.utcnow().strftime(
        '%Y-%m-%dT%H:%M:%S'), 'file creation date'

    head['HIERARCH SPATPIXL'] = axdict['imgstart'], 'lower limit of spatial axis, pix'
    head['HIERARCH SPATPIXH'] = axdict['imgend'], 'upper limit of spatial axis, pix'
    head['HIERARCH DISPPIXL'] = axdict['wavstart'], 'lower limit of dispersion axis, pix'
    head['HIERARCH DISPPIXH'] = axdict['wavend'], 'upper limit of dispersion axis, pix'
    head['HIERARCH WAVL'] = np.floor(axdict['waxis'][0]), 'lower limit of wav range, ' + hparams['wavunit']
    head['HIERARCH WAVH'] = np.ceil(axdict['waxis'][-1]), 'upper limit of wav range, ' + hparams['wavunit']

    head['HIERARCH MOFF A'] = round(moffpars[0], 5), 'moffat profile amplitude'
    head.add_blank('Parameters fit to the median spatial profile of the spectrum',
                   before='HIERARCH MOFF A')
    head['HIERARCH MOFF C'] = round(moffpars[1] + axdict['imgstart'], 5), 'moffat profile centre'
    head['HIERARCH MOFF ALPHA'] = round(moffpars[2], 5), 'moffat profile alpha value'
    head['HIERARCH MOFF BETA'] = round(moffpars[3], 5), 'moffat profile beta value'
    head['HIERARCH MOFF BACK'] = round(moffpars[4], 5), 'moffat profile background level'
    head['HIERARCH MOFF GRAD'] = round(moffpars[5], 5), 'moffat profile background slope'
    head['HIERARCH IQ'] = round(hparams['seeing']*hparams['pixresolution'], 2), 'IQ measured from median profile, "'
    
    head['HIERARCH SNR BIN LIMIT'] = pars['-SNR_BIN_LIM'], 'maximum SNR per bin'
    head.add_blank('Dispersion Binning and Spectrum Extraction', before='HIERARCH SNR BIN LIMIT')
    head['HIERARCH COL BIN LIMIT'] = int(pars['-COL_BIN_LIM']), 'minimum number of columns per bin'
    head['HIERARCH FWHM MULTIPLIER'] = pars['-FWHM_MULTIPLIER'], 'FWHM used to define the extraction limits'
    head['HIERARCH INTERP KIND'] = pars['-INTERP_KIND'], 'interpolation mode used'

    if pars['-SUBTRACT_SKY']:
        head['HIERARCH SKYSUB FWHM MULT'] = pars['-BG_FWHM_MULTIPLIER'], 'FWHM multiplier for defining background'
        head.add_blank('Sky Subtraction', before='HIERARCH SKYSUB FWHM MULT')
        head['HIERARCH SKYSUB SNR BIN LIM'] = pars['-SKY_SNR_BIN_LIM'], 'max SNR per bin for sky subtraction'
        skymodhdu = fits.ImageHDU(fdict['skymod'])
        skymodhdu.header['EXTNAME'] = '2D_SKY'
        skybinhdu = fits.ImageHDU(sbpars)
        skybinhdu.header['EXTNAME'] = 'SKY_BIN_PARS'
        skyextractionlims = fits.ImageHDU(skyextractionlims)
        skyextractionlims.header['EXTNAME'] = 'SKY_EXT_LIMS'

    if pars['-MASK_CR']:
        head['HIERARCH CR CLIP ITER'] = int(pars['-CR_CLIP_ITER']), 'cosmic ray sigma-clipping iterations'
        head.add_blank('Cosmic Ray Masking and Bad Pixel Replacement', before='HIERARCH CR CLIP ITER')
        head['HIERARCH CR SIGMA MULT'] = pars['-CR_CLIP_STD_MULT'], 'cosmic ray sigma-clipping multiplier'
        head['HIERARCH CR FWHM MULT'] = pars['-CR_FWHM_MULTIPLIER'], 'FWHM multiplier for defining CR clip regions'
        if pars['-REPLACE_CRBP']:
            head['HIERARCH CR REPLACED?'] = 'YES'
        else:
            head['HIERARCH CR REPLACED?'] = 'NO'
        crmaskhdu = fits.ImageHDU(fdict['cmask'])
        crmaskhdu.header['EXTNAME'] = '2D_CR_MASK'
    
    head['HIERARCH EXTRACTED HDU ROW 0'] = 'Wavelength Axis, ' + hparams['wavunit']
    head.add_blank('Data Saved in the Extracted Spectrum HDU', before='HIERARCH EXTRACTED HDU ROW 0')
    head['HIERARCH EXTRACTED HDU ROW 1'] = 'Flux, ' + hparams['fluxunit']
    head['HIERARCH EXTRACTED HDU ROW 2'] = 'Flux Uncertainty, ' + hparams['fluxunit']
    head['EXTNAME'] = '1D_SPEC'
    
    fluxhdu = fits.PrimaryHDU([axdict['waxis'], flux, errs], header=head)
    spec2Dhdu = fits.ImageHDU(fdict['ogdata'])
    spec2Dhdu.header['EXTNAME'] = 'ORIG_2D_SPEC'
    errs2Dhdu = fits.ImageHDU(fdict['ogerrs'])
    errs2Dhdu.header['EXTNAME'] = 'ORIG_2D_ERRS'
    qual2Dhdu = fits.ImageHDU(fdict['ogqual'])
    qual2Dhdu.header['EXTNAME'] = 'ORIG_2D_QUAL'
    binhdu = fits.ImageHDU(bpars)
    binhdu.header['EXTNAME'] = 'EXT_BIN_PARS'
    extractionlims = fits.ImageHDU(extractionlims)
    extractionlims.header['EXTNAME'] = 'EXT_LIMS'
    hdu_list = [fluxhdu, spec2Dhdu, errs2Dhdu, qual2Dhdu, binhdu, extractionlims]
    
    if pars['-SUBTRACT_SKY']:
        hdu_list.append(skymodhdu)
        hdu_list.append(skybinhdu)
        hdu_list.append(skyextractionlims)
    if pars['-MASK_CR']:
        hdu_list.append(crmaskhdu)

    hdulist = fits.HDUList(hdu_list)
    filenamelist = filename.split('_')
    hdulist.writeto('m' + '_'.join(filenamelist[0:-1]) + '_' + filenamelist[-1])
    hdulist.close()

    sys.stdout.write(' >>> Spectrum extracted and saved:\n')
    sys.stdout.write('     ' + os.getcwd() + '/' + 'm' + '_'.join(filenamelist[0:-1]) + '_' + filenamelist[-1] +'\n')
    return None
    
    
# ////////////////////////////////////////////////////////////////////////////#
# SKY LOCALISATION AND SUBTRACTION
def skyloc(framedict, axesdict, datascale, headparams, binparams, params):
	# For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
    # moffat profile to the median data and then use the parameters of the fitted Moffat function to localise the 2D spectrum.
    sys.stdout.write(' >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n')

    skybin = []
    extractionlimits = []
    for bin in binparams:
		# Take the median spatial profile of the dispersion 
		# bin, and leave out pixel columns in the chip gaps if this is a GMOS spectrum.
        binimg = framedict['data'][:, bin[0]:bin[1]]
        chipgap = np.where(np.median(binimg, axis=0)!=1)
        bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)

        # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median spatial profile and
        # return its parameters.
        binmoffparams = common.moffat_least_squares(axesdict['saxis'], bindata*datascale, headparams['seeing'], headparams['pixresolution'])

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle (see *1E18 above).
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled down again after the fitting.
        binmoffparams[0] /= datascale
        binmoffparams[4] /= datascale
        binmoffparams[5] /= datascale

        # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
        # previously fitted to it.
        LowExt, HighExt, fwhm, centre = common.extraction_limits(binmoffparams, params['-BG_FWHM_MULTIPLIER'], axesdict, highres=True)

        extractionlimits.append([(bin[0] + bin[1]) * 0.5, LowExt, HighExt, centre])

        # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
        # locations so they can be saved as metadata along with the extracted spectrum.
        binmoffparams.append(bin[0]+axesdict['wavstart'])
        binmoffparams.append(bin[1]+axesdict['wavstart'])
        skybin.append(binmoffparams)

        # DIAGNOSTICS - Plot computed moffat profile over data for each bin
        if params['-DIAG_PLOT_MOFFAT']:
            common.plot_fitted_spatial_profile(axesdict['saxis'], bindata, axesdict['hrsaxis'], binmoffparams, axesdict['imgstart'], headparams)

    skybin = np.array(skybin)

    sys.stdout.write('     Fitting complete.\n')

    sys.stdout.write(' >>> Drawing target/sky boundaries. ')
    sys.stdout.flush()
    extractionlimits = np.array(extractionlimits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [extractionlimits[0] + axesdict['wavstart'], (extractionlimits[1] * 0.02) + axesdict['imgstart'],
                     extractionlimits[0] + axesdict['wavstart'], (extractionlimits[2] * 0.02) + axesdict['imgstart']]
        
        common.show_img(framedict['data'], axesdict, headparams, drawlines, '2D Spectrum Overplotted with Target/Sky Boundaries')

    # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
    # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.
    skyextractionlims = common.interpolate_extraction_lims(extractionlimits, axesdict['dispaxislen'], params['-INTERP_KIND'])

    skyextractionlims[0] *= 0.02
    skyextractionlims[1] *= 0.02

    sys.stdout.write('DONE.\n')

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[0]) + axesdict['imgstart'],
                     np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[1]) + axesdict['imgstart']]
        
        common.show_img(framedict['data'], axesdict, headparams, drawlines, '2D Spectrum Overplotted with Full Target/Sky Boundaries')

    sys.stdout.write(' >>> Subtracting sky. ')
    sys.stdout.flush()

    framedict = common.subtract_sky(skyextractionlims[0], skyextractionlims[1], framedict, axesdict, crmask_multiplier=1.)

    sys.stdout.write('DONE.\n')
    sys.stdout.flush()

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[0]) + axesdict['imgstart'],
                     np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[1]) + axesdict['imgstart']]
        
        common.show_img(framedict['data'], axesdict, headparams, drawlines, 'Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries')
	
    return framedict, skybin, skyextractionlims
	

###############################################################################
# MAIN ///////////////////////////////////////////////////////////////////////#
###############################################################################

if __name__ == '__main__':
    motes()
