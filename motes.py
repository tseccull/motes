###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import astropy.io.fits as fits
import datetime
import glob
import motes.common as common
import motes.harvester as harvester
import motes.startup as startup
import numpy as np
import os
import sys

# Force TkAgg backend plots
import matplotlib
matplotlib.use('TkAgg', force=True)
import matplotlib.pyplot as plt

###############################################################################
# FUNCTIONS //////////////////////////////////////////////////////////////////#
###############################################################################
# ////////////////////////////////////////////////////////////////////////////#
# HANDLE COSMIC RAYS AND BAD PIXELS


def cr_handling(pars, axdict, lext, hext, fdict, hpars):
    sys.stdout.write(' >>> Separating spectrum and background regions for cosmic ray detection.\n')

    # DIAGNOSTICS - Plot limits defined for CR masking process.
    if pars['-DIAG_PLOT_CRMASK']:
        daxis = np.linspace(axdict['wavstart'], axdict['wavstart'] +
                            len(axdict['waxis']), num=len(axdict['waxis']))
        loline = [lext + axdict['imgstart']] * len(axdict['waxis'])
        hiline = [hext + axdict['imgstart']] * len(axdict['waxis'])
        dlines = np.array([daxis, loline, daxis, hiline])
        common.show_img(fdict['data'], axdict, hpars, dlines,
                        '2D Spectrum Mask Regions - Pre Masking')

    sys.stdout.write(' >>> Masking cosmic rays via sigma-clip:\n')
    sys.stdout.write('     Sigma multiplier = ' + str(pars['-CR_CLIP_STD_MULT']) + '\n')
    sys.stdout.write('     Number of sigma-clip iterations = ' + str(pars['-CR_CLIP_ITER']) + '\n')
    sys.stdout.write(' >>> Masking bad pixels.\n')
    fdict['cmask'] = np.ones(np.shape(fdict['data']))
    for it in range(int(pars['-CR_CLIP_ITER'])):
        crmasque = common.mask_cosmic_rays(
            fdict['data'], lext, hext, multiplier=pars['-CR_CLIP_STD_MULT'])
        fdict['qual'] *= crmasque
        fdict['data'] *= fdict['qual']
        fdict['errs'] *= fdict['qual']
        fdict['cmask'] *= crmasque
    sys.stdout.write(' >>> Cosmic rays and bad pixels masked.\n')

    # DIAGNOSTICS - Plot 2D spectrum with CRs masked.
    if pars['-DIAG_PLOT_CRMASK']:
        common.show_img(fdict['data'], axdict, hpars, dlines,
                        '2D Spectrum Mask Regions - Post Masking')

    # Determine the location of boundaries between dispersion bins in the 2D spectrum that will be fitted with Moffat
    # functions to precisely localise the spectrum to allow sky subtraction.
    # If cosmic rays have been masked, and the user wants them to be replaced with appropriate flux values the
    # process is performed while determining the location of the dispersion bins. See the docs or comments in the
    # common.get_bins function for more information.
    # Replacement is mostly, but not entirely effective. It may cause particular problems when used on fainter
    # spectra, or those with very bright sky lines. Some fiddling with the CR sigma-clipping multipliers may be
    # required in the sky subtraction and CR masking functions. Use with care.

    # Determine the boundaries of the localisation bins on the dispersion axis
    binpars, fdict = common.get_bins(
        fdict, int(
            np.floor(lext)), int(
            np.ceil(hext)), axdict['dispaxislen'], pars, sky=True)

    if pars['-REPLACE_CRBP']:
        sys.stdout.write(' >>> Cosmic rays and bad pixels replaced.\n')

        # DIAGNOSTICS - Show 2D spectrum with bad pixels and CRs replaced.
        if pars['-DIAG_PLOT_CRREPLACE']:
            common.show_img(fdict['data'], axdict, hpars, dlines,
                            '2D Spectrum with Cosmic Rays and Bad Pixels Replaced')

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

        sys.stdout.write(('/' * (70 - len(file_2D[:70]))) + ' ' + file_2D[:70] + '\n')
        sys.stdout.write(' >>> Beginning MOTES Processing\n')

        # GATHER HEADER INFO AND DATA FROM THE 2D IMAGE FILE
        sys.stdout.write(' >>> Gathering image frames and header data from input file.\n')
        headparams, framedict, axesdict, imghead = harvester.data_harvest(
            i, file_2D, intreg, params)

        # Save the original data
        framedict['ogdata'] = framedict['data']
        framedict['ogerrs'] = framedict['errs']
        sys.stdout.write(' >>> Gathering image frames and header data from input file completed.\n')
        sys.stdout.write(
            ' >>> Fitting Moffat profile to median spatial profile of entire spectrum. ')
        sys.stdout.flush()

        # Perform an initial least-squares Moffat fitting over the entire 2D image
        # collapsed along the dispersion axis with a median.
        datadispcollapse = np.nanmedian(framedict['data'], axis=1)
        # For a more accurate fit of the spatial profile, introduce the Median
        # Absolute Deviation (MAD) as the measure of spread of each spatial column flux distribution
        # in the spectrogram.
        errdispcollapse = np.empty((0))
        for row in range(framedict['data'].shape[0]):
            mad = np.nanmedian(np.abs(framedict['data'][row, :] -
                                      np.nanmedian(framedict['data'][row, :])))
            errdispcollapse = np.append(errdispcollapse, mad)

        # In NOODING mode of the X-Shooter, the scale can include negatives use
        # np.abs() to the 86th percentile of the dataiscollapsed to ensure its positive.
        # In version 0.2.0 np.percentile() is replaced with np.nanmedian().
        # However, specifically in NODDING spectrograms there is a possibility of
        # a median being closer zero. To ensure that datascale is large, using 86%
        # of the dataiscollapse yields better results.
        datascale = 10**np.abs(np.floor(np.log10(np.abs(np.percentile(datadispcollapse, 86)))))
        print('Datascale: {0:.1E}'.format(datascale))

        # Perform the fitting here over a new routine
        moffparams, _ = common.moffat_weighted_curvefit(
            axesdict['saxis'],
            datadispcollapse * datascale,
            headparams['seeing'],
            errdispcollapse * datascale, params['-ABBA'])

        # Obtain the (new, measured?) seeing from the profile parameters.
        headparams['seeing'] = 2 * moffparams[2] * np.sqrt((2**(1 / moffparams[3])) - 1)

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle.
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled
        # down again after the fitting.
        moffparams[0] /= datascale
        moffparams[4] /= datascale
        moffparams[5] /= datascale

        sys.stdout.write('DONE.\n')
        sys.stdout.write(' >>> FWHM of median spatial profile is ' +
                         str(round(headparams['seeing'], 2)) +
                         ' Pixels, or ' +
                         str(round(headparams['seeing'] *
                                   headparams['pixresolution'], 2)) +
                         '"\n')

        # DIAGNOSTICS -  Plot fitted Moffat profile over collapsed 2D spectrum and
        # print the parameters of the fitted Moffat profile.
        if params['-DIAG_PLOT_COLLAPSED_2D_SPEC']:
            common.printmoffparams(moffparams, axesdict['imgstart'], datascale)
            common.plot_fitted_spatial_profile(
                axesdict['saxis'],
                datadispcollapse,
                axesdict['hrsaxis'],
                moffparams,
                axesdict['imgstart'],
                headparams)

        # Use the parameters of the Moffat profile fitted to the median spatial profile of the entire spectrum to determine
        # spatial limits that are used to bound the region of the spectrum used by the common.get_bins() function to
        # to bin the 2D spectrum while taking account of its S/N.

        # These spatial limits may also be used during the process of masking,
        # removing, and replacing cosmic rays and bad pixels.
        lowext, highext, fwhm, cent = common.extraction_limits(
            moffparams, params['-CR_FWHM_MULTIPLIER'], axesdict)
        sys.stdout.write(' >>> Spectrum localised to aperture in range of spatial pixel rows ' + str(
            int(lowext + axesdict['imgstart'])) + '-' + str(int(highext + axesdict['imgstart'])) + '\n')

        # Save the (params['-CR_FWHM_MULTIPLIER']) * FWHM limits of median profile, of the whole collapsed spectra for
        # re-estimating the bin sizes after sky subtraction. Otherwise, the last
        # bin of the sky subtraction bin loop is used in S/N estimations (in
        # common.get_bins()) in later stages of the code.
        lowext_median = lowext
        highext_median = highext

        """
        ///////// Cosmic Ray/Bad Pixel Masking /////////
        """

        # IF (for COSMICS) ######################################################

        # Mask and remove cosmic rays and bad pixels and, if selected, replace them in the image.
        if params['-MASK_CR']:
            binparams, framedict = cr_handling(
                params, axesdict, lowext_median, highext_median, framedict, headparams)

        # No cosmic ray masking requested, just get localisation bins
        else:
            sys.stdout.write(' >>> Cosmic ray masking deactivated and not performed.\n')
            framedict['cmask'] = np.ones(np.shape(framedict['data']))
            binparams, framedict = common.get_bins(
                framedict, int(
                    np.floor(lowext_median)), int(
                    np.ceil(highext_median)), axesdict['dispaxislen'], params, sky=True)
            # Print and plot common.get_bins
            common.get_bins_output(
                binparams,
                params,
                axesdict['saxis'],
                lowext_median,
                highext_median,
                axesdict['wavstart'],
                axesdict['imgstart'],
                framedict['data'],
                axesdict['waxis'],
                headparams)
            sys.stdout.write(' >>> Bad pixels replaced.\n')

        #############################
        # Note: lowext_median and highext_median were used here for common.get_bins() and common.mask_cosmic_rays(). For now, these limits are determined by the user via (params['-CR_FWHM_MULTIPLIER']) over the whole spectrum's median spatial profile.
        #############################

        """
        ///////// Sky Subtraction /////////
        """

        # SUBTRACT THE SKY SPECTRUM IF REQUESTED BY THE USER.
        if params['-SUBTRACT_SKY']:
            framedict, skybinpars, skyextlims = skyloc(
                framedict, axesdict, datascale, headparams, binparams, params)

        # Last use of common.get_bins() to update the bins and the data before extraction.
        # Change lowext (first bin result) to lowext_median from the first median collapsed spectrum.
        # Otherwise, it uses the last fitted bin in the sky loop) which is not be
        # correct - probably deprecated, now this is separate in skyloc()

        binparams, framedict = common.get_bins(
            framedict, int(
                np.floor(lowext_median)), int(
                np.ceil(highext_median)), axesdict['dispaxislen'], params)
        # Print and plot common.get_bins
        common.get_bins_output(
            binparams,
            params,
            lowext_median,
            highext_median,
            framedict['data'],
            headparams,
            axesdict)

        """
        ///////// Extraction /////////
        """

        # For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
        # moffat profile to the median data and then use the parameters of the fitted Moffat function to localise the 2D
        # spectrum.
        sys.stdout.write(' >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n')

        extbin = []
        extractionlimits_optimal = []
        p_value_extraction = np.empty((0))

        for bin in binparams:

            # Filter out any nans or inf and ensure the errors are positive
            framedict['data'], framedict['errs'] = common.filter_data(
                framedict['data'], np.abs(framedict['errs']))

            # Leave out pixel columns in the chip gaps if this is a GMOS spectrum.
            binimg = framedict['data'][:, bin[0]:bin[1]]
            chipgap = np.where(np.median(binimg, axis=0) != 1)

            # Use a weighted Levenberg-Marquardt Least Squares method to fit a Moffat function to the (median) spatial profile and
            # return its parameters. Determine the median absolute deviation as a measure of spread.
            bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)
            binerrs = np.empty((0))
            for row in range(bindata.shape[0]):
                mad = np.nanmedian(np.abs(binimg[row, chipgap[0]] -
                                          np.nanmedian(binimg[row, chipgap[0]])))
                binerrs = np.append(binerrs, mad)
            binmoffparams, pval = common.moffat_weighted_curvefit(
                axesdict['saxis'],
                bindata * datascale,
                headparams['seeing'],
                binerrs * datascale, params['-ABBA'])
            p_value_extraction = np.append(p_value_extraction, pval)

            '''
            # Take the median spatial profile of the dispersion
            # bin, and leave out pixel columns in the chip gaps if this is a GMOS spectrum.
            binimg = framedict['data'][:, bin[0]:bin[1]]
            chipgap = np.where(np.median(binimg, axis=0) != 1)
            bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)

            # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median spatial profile and
            # return its parameters.
            binmoffparams = common.moffat_least_squares(
                axesdict['saxis'],
                bindata * datascale,
                headparams['seeing'],
                headparams['pixresolution'])
            '''

            # Again, descalling variables
            binmoffparams[0] /= datascale
            binmoffparams[4] /= datascale
            binmoffparams[5] /= datascale

            # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
            # previously fitted to it. NOTE that highres=True
            lowext_optimal, highext_optimal, fwhm_optimal, centre_optimal = common.extraction_limits(
                binmoffparams, params['-FWHM_MULTIPLIER'], axesdict, highres=True)
            extractionlimits_optimal.append(
                [(bin[0] + bin[1]) * 0.5, lowext_optimal, highext_optimal, centre_optimal])

            # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
            # locations so they can be saved as metadata along with the extracted spectrum.
            binmoffparams.append(bin[0] + axesdict['wavstart'])
            binmoffparams.append(bin[1] + axesdict['wavstart'])
            extbin.append(binmoffparams)

            # DIAGNOSTICS - Plot computed moffat profile over data for each bin
            if params['-DIAG_PLOT_MOFFAT']:
                common.plot_fitted_spatial_profile(
                    axesdict['saxis'],
                    bindata,
                    axesdict['hrsaxis'],
                    binmoffparams,
                    axesdict['imgstart'],
                    headparams)

        # Renamed to bin_save for .FITS file saving
        bin_save = np.array(extbin)
        sys.stdout.write('     Fitting complete.\n')
        sys.stdout.write(' >>> Drawing extraction aperture limits. ')
        sys.stdout.flush()

        # Transpose extraction limits
        extractionlimits_optimal = np.array(extractionlimits_optimal).T

        '''
        # NEW DIAGNOSTIC - Plot the p_values in extraction bins
        plt.plot(p_value_extraction, marker='.', linestyle=' ', color='k')
        plt.title('p-value')
        plt.xlabel('Bin')
        plt.ylabel('p-value')
        plt.show()
        print '\n Fitting Median p-Value: ' + str(np.nanmedian(p_value_extraction))
        '''

        # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
        if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
            drawlines = [extractionlimits_optimal[0] + axesdict['wavstart'], (extractionlimits_optimal[1] * 0.02) + axesdict['imgstart'],
                         extractionlimits_optimal[0] + axesdict['wavstart'], (extractionlimits_optimal[2] * 0.02) + axesdict['imgstart']]
            common.show_img(
                framedict['data'],
                axesdict,
                headparams,
                drawlines,
                '2D Spectrum Overplotted with Extraction Limits')

        # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
        # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.
        # This also includes extrapolated regions!
        finalextractionlims_optimal = common.interpolate_extraction_lims(
            extractionlimits_optimal, axesdict['dispaxislen'], params['-INTERP_KIND'])
        sys.stdout.write('DONE.\n')

        # DIAGNOSTICS - Plot the final extraction limits including the
        # extrapolated sections at the ends of the wavelength axis. Ensure the the
        # final_extractionlims are converted back to original spatial resoltion.

        # For the drawlines variable, np.floor(finalextractionlims_optimal[0]) and np.ceil(finalextractionlims_optimal[1]).
        # This is in line with the current version of optimal extraction that uses
        # non-supersampled data.

        if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
            drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (np.floor(finalextractionlims_optimal[0] * 0.02)) + axesdict['imgstart'],
                         np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (np.ceil(finalextractionlims_optimal[1] * 0.02)) + axesdict['imgstart']]
            common.show_img(
                framedict['data'],
                axesdict,
                headparams,
                drawlines,
                '2D Spectrum Overplotted with Full Extraction Limits')

        """
        ///////// Supersampling Spatial Limits /////////
        """

        # TODO Finish implementing in 0.4.0_dev?
        # Supersampling below is not quite right and in fact breaks current edition of
        # optimal extraction but failing to plot the bins.
        '''
        sample_numer = 50
        # Supersample data in the spatial direction (i.e. rows)
        # Data itself is simple, supersampling in spatial direction and ensuring
        # that the total flux is not biased by having a multiplier reciprocal to
        # sample number.
        supersamp_data2D = (np.repeat(framedict['data'], sample_numer, axis=0) * (1 / sample_numer))

        # The original error in one pixel needs to be split into subpixels those propagated sum (of errors) returns the original error.
        # Therefore, one should reverse standard error propagation (assuming Gaussian errors as it is the case with ESO instrument results) to obtain a non-biased total error (our non-supersampled original error)
        # In order to do that, one therefore could find a scalar factor to divide
        # the entire non-supersampled error frame by (such factor will depend on
        # the sample size) then supersample without a further multiplicative
        # factor.
        corrected_error_frame = np.divide(framedict['errs'], np.sqrt(sample_numer))
        supersamp_errs2D = (np.repeat(corrected_error_frame, sample_numer, axis=0))
        '''

        """
        ///////// Optimal Extraction /////////
        """

        # Convert to normal resolution prior to feeding it to feeding it optimal extraction routine.
        finalextractionlims_optimal[0] *= 0.02
        finalextractionlims_optimal[1] *= 0.02

        sys.stdout.write(' >>> Optimally Extracting 1D spectrum. ')
        sys.stdout.flush()

        # Run the optimal extraction module.
        optdata1D, opterrs1D, stddata1D, stderrs1D = common.optimal_extraction(
            framedict['data'], framedict['errs'], extractionlimits=finalextractionlims_optimal, binparameters=binparams, datascale=datascale, seeing=headparams['seeing'])

        """
        ///////// Plot /////////
        """

        # DIAGNOSTICS - Plot extracted spectrum.
        if params['-PLOT_EXTRACTED_SPECTRUM']:
            def plot_step_spec(x, y, err, lb, colour):
                plt.step(x,  # add units to main plot
                         y,
                         marker='.',
                         linestyle=' ',
                         markersize=1,
                         where='mid',
                         linewidth=0.5,
                         label=lb,
                         color=colour)
                plt.errorbar(x,
                             y,  # Units are not supported for error bars
                             err,
                             linestyle=' ',
                             ecolor=colour,
                             capsize=0.5,
                             zorder=1)

            fig = plt.figure(figsize=(20, 10))
            ax = fig.add_subplot(111)
            plot_step_spec(axesdict['waxis'], optdata1D, opterrs1D, 'Optimal', 'blue')
            plot_step_spec(axesdict['waxis'], stddata1D, stderrs1D, 'Standard', 'green')
            plt.grid(alpha=0.5, linestyle='dotted')
            plt.title('Extracted Optimal 1D Spectrum')
            plt.ylabel('Flux, ' + headparams['fluxunit'])
            plt.xlabel('Wavelength, ' + headparams['wavunit'])
            plt.show()

        """
        ///////// Saving /////////
        """
        # Add predefined values for sbpars and skyextractionlims to avoid local
        # runtime errors if sky was not subtracted. This avoids issues when saving
        # (and loading) the .fits file.
        if not params['-SUBTRACT_SKY']:
            skybinpars = 0.0
            skyextlims = 0.0

        # Save the extracted spectrum to a new .fits file.
        if params['-SAVE']:
            sys.stdout.write(' >>> Saving 1D spectrum and metadata.\n')
            sys.stdout.flush()
            save_fits(
                axesdict,
                headparams,
                optdata1D,
                opterrs1D,
                imghead,
                params,
                file_2D,
                moffparams,
                framedict,
                bin_save,
                finalextractionlims_optimal,
                skybinpars,
                skyextlims)

        sys.stdout.write(' >>> Extraction of ' + file_2D + ' completed.\n')
    sys.stdout.write(' >>> MOTES Processing Complete.\n\n')

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function saves the extracted spectrum and intermediate products in a single FITS file.
def save_fits(axdict, hparams, flux, errs, head, pars, filename, moffpars,
              fdict, bpars, extractionlims, sbpars, skyextractionlims):

    head['MOTES'] = '######## Extracted 1D Spectrum Metadata ########'
    head.add_blank('', before='MOTES')
    head['HIERARCH UTC EXT DATE'] = datetime.datetime.utcnow().strftime(
        '%Y-%m-%dT%H:%M:%S'), 'file creation date'

    head['HIERARCH SPATPIXL'] = axdict['imgstart'], 'lower limit of spatial axis, pix'
    head['HIERARCH SPATPIXH'] = axdict['imgend'], 'upper limit of spatial axis, pix'
    head['HIERARCH DISPPIXL'] = axdict['wavstart'], 'lower limit of dispersion axis, pix'
    head['HIERARCH DISPPIXH'] = axdict['wavend'], 'upper limit of dispersion axis, pix'
    head['HIERARCH WAVL'] = np.floor(
        axdict['waxis'][0]), 'lower limit of wav range, ' + hparams['wavunit']
    head['HIERARCH WAVH'] = np.ceil(
        axdict['waxis'][-1]), 'upper limit of wav range, ' + hparams['wavunit']

    head['HIERARCH MOFF A'] = round(moffpars[0], 5), 'moffat profile amplitude'
    head.add_blank('Parameters fit to the median spatial profile of the spectrum',
                   before='HIERARCH MOFF A')
    head['HIERARCH MOFF C'] = round(moffpars[1] + axdict['imgstart'], 5), 'moffat profile centre'
    head['HIERARCH MOFF ALPHA'] = round(moffpars[2], 5), 'moffat profile alpha value'
    head['HIERARCH MOFF BETA'] = round(moffpars[3], 5), 'moffat profile beta value'
    head['HIERARCH MOFF BACK'] = round(moffpars[4], 5), 'moffat profile background level'
    head['HIERARCH MOFF GRAD'] = round(moffpars[5], 5), 'moffat profile background slope'
    head['HIERARCH IQ'] = round(hparams['seeing'] *
                                hparams['pixresolution'], 2), 'IQ measured from median profile, "'

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
        head['HIERARCH CR CLIP ITER'] = int(
            pars['-CR_CLIP_ITER']), 'cosmic ray sigma-clipping iterations'
        head.add_blank(
            'Cosmic Ray Masking and Bad Pixel Replacement',
            before='HIERARCH CR CLIP ITER')
        head['HIERARCH CR SIGMA MULT'] = pars['-CR_CLIP_STD_MULT'], 'cosmic ray sigma-clipping multiplier'
        head['HIERARCH CR FWHM MULT'] = pars['-CR_FWHM_MULTIPLIER'], 'FWHM multiplier for defining CR clip regions'
        if pars['-REPLACE_CRBP']:
            head['HIERARCH CR REPLACED?'] = 'YES'
        else:
            head['HIERARCH CR REPLACED?'] = 'NO'
        crmaskhdu = fits.ImageHDU(fdict['cmask'])
        crmaskhdu.header['EXTNAME'] = '2D_CR_MASK'

    head['HIERARCH EXTRACTED HDU ROW 0'] = 'Wavelength Axis, ' + hparams['wavunit']
    head.add_blank('Data Saved in the Extracted Spectrum HDU',
                   before='HIERARCH EXTRACTED HDU ROW 0')
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
    sys.stdout.write('     ' + os.getcwd() + '/' + 'm' +
                     '_'.join(filenamelist[0:-1]) + '_' + filenamelist[-1] + '\n')
    return None


# ////////////////////////////////////////////////////////////////////////////#
# SKY LOCALISATION AND SUBTRACTION

def skyloc(framedict, axesdict, datascale, headparams, binparams, params):

    # For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
    # moffat profile to the median data and then use the parameters of the
    # fitted Moffat function to localise the 2D spectrum.
    sys.stdout.write(' >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n')

    skybin = []
    extractionlimits_sky = []
    p_value_sky = np.empty((0))

    for bin in binparams:

        # Filter out any nans or inf and ensure the errors are positive
        framedict['data'], framedict['errs'] = common.filter_data(
            framedict['data'], np.abs(framedict['errs']))

        # Leave out pixel columns in the chip gaps if this is a GMOS spectrum.
        binimg = framedict['data'][:, bin[0]:bin[1]]
        chipgap = np.where(np.median(binimg, axis=0) != 1)

        # Use a weighted Levenberg-Marquardt Least Squares method to fit a Moffat function to the (median) spatial profile and
        # return its parameters. Determine the median absolute deviation as a measure of spread.
        bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)
        binerrs = np.empty((0))
        for row in range(bindata.shape[0]):
            mad = np.nanmedian(np.abs(binimg[row, chipgap[0]] -
                                      np.nanmedian(binimg[row, chipgap[0]])))
            binerrs = np.append(binerrs, mad)
        binmoffparams, pval = common.moffat_weighted_curvefit(
            axesdict['saxis'],
            bindata * datascale,
            headparams['seeing'],
            binerrs * datascale, params['-ABBA'])
        p_value_sky = np.append(p_value_sky, pval)

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle (see *1E18 above).
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled
        # down again after the fitting.
        binmoffparams[0] /= datascale
        binmoffparams[4] /= datascale
        binmoffparams[5] /= datascale

        # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
        # previously fitted to it.
        # Previosly, the use of common.extraction_limits overwrote lowext and highext.
        # Rename them to lowext_sky / highext_sky instead.
        lowext_sky, highext_sky, fwhm_sky, centre_sky = common.extraction_limits(
            binmoffparams, params['-BG_FWHM_MULTIPLIER'], axesdict, highres=True)
        extractionlimits_sky.append(
            [(bin[0] + bin[1]) * 0.5, lowext_sky, highext_sky, centre_sky])

        # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
        # locations so they can be saved as metadata along with the extracted spectrum.
        binmoffparams.append(bin[0] + axesdict['wavstart'])
        binmoffparams.append(bin[1] + axesdict['wavstart'])
        skybin.append(binmoffparams)

        # DIAGNOSTICS - Plot computed moffat profile over data for each bin
        if params['-DIAG_PLOT_MOFFAT']:
            common.plot_fitted_spatial_profile(
                axesdict['saxis'],
                bindata,
                axesdict['hrsaxis'],
                binmoffparams,
                axesdict['imgstart'],
                headparams)

    # Bin fitting loop finished.
    # Convert the skybin into an array to save as metadata(?) and transpose
    # the 'extractionlimits_sky' array
    sys.stdout.write('     Fitting complete.\n')
    skybin = np.array(skybin)
    extractionlimits_sky = np.array(extractionlimits_sky).T
    sys.stdout.write(' >>> Drawing target/sky boundaries. ')
    sys.stdout.flush()

    '''
    # NEW DIAGNOSTIC - Plot the p_values in sky bins
    plt.plot(p_value_sky, marker='.', linestyle=' ', color='k')
    plt.title('p-value')
    plt.xlabel('Bin')
    plt.ylabel('p-value')
    plt.show()
    print '\n Fitting Median p-Value: ' + str(np.nanmedian(p_value_sky))
    '''

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [extractionlimits_sky[0] + axesdict['wavstart'], (extractionlimits_sky[1] * 0.02) + axesdict['imgstart'],
                     extractionlimits_sky[0] + axesdict['wavstart'], (extractionlimits_sky[2] * 0.02) + axesdict['imgstart']]
        common.show_img(framedict['data'], axesdict, headparams, drawlines,
                        '2D Spectrum Overplotted with Target/Sky Boundaries')

    # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
    # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.

    # Remember extractionlimits_sky != skyextractionlims. The latter are final interpolated limits.
    skyextractionlims = common.interpolate_extraction_lims(
        extractionlimits_sky, axesdict['dispaxislen'], params['-INTERP_KIND'])

    # Convert back to standard resolution for plotting and common.subtract_sky().
    skyextractionlims[0] *= 0.02
    skyextractionlims[1] *= 0.02

    # Use a single sky bin instead using the a median value of the determined sky limits.
    if params['-SUBTRACT_SKY_SINGLE_BIN']:
        med_sky_0_array = np.repeat(
            np.nanmedian(
                skyextractionlims[0]), len(
                skyextractionlims[0]))
        med_sky_1_array = np.repeat(
            np.nanmedian(
                skyextractionlims[1]), len(
                skyextractionlims[1]))
        skyextractionlims[0] = med_sky_0_array
        skyextractionlims[1] = med_sky_1_array
    sys.stdout.write('DONE.\n')

    # DIAGNOSTICS - Plot the final extraction limits including the
    # extrapolated sections at the ends of the wavelength axis.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[0]) + axesdict['imgstart'],
                     np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[1]) + axesdict['imgstart']]
        common.show_img(framedict['data'], axesdict, headparams, drawlines,
                        '2D Spectrum Overplotted with Full Target/Sky Boundaries')

    # Actual sky subtraction takes place here.
    sys.stdout.write(' >>> Subtracting sky. ')
    sys.stdout.flush()
    framedict = common.subtract_sky(
        skyextractionlims[0],
        skyextractionlims[1],
        framedict,
        axesdict)
    sys.stdout.write('DONE.\n')
    sys.stdout.flush()

    # DIAGNOSTICS - Plot the final extraction limits including the
    # extrapolated sections at the ends of the wavelength axis.
    if params['-DIAG_PLOT_EXTRACTION_LIMITS']:
        drawlines = [np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[0]) + axesdict['imgstart'],
                     np.array(range(axesdict['dispaxislen'])) + axesdict['wavstart'], (skyextractionlims[1]) + axesdict['imgstart']]
        common.show_img(
            framedict['data'],
            axesdict,
            headparams,
            drawlines,
            'Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries')

    return framedict, skybin, skyextractionlims


###############################################################################
# MAIN ///////////////////////////////////////////////////////////////////////#
###############################################################################

if __name__ == '__main__':
    motes()
