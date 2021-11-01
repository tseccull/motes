###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import copy
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import scipy.interpolate as interp
import sys
import warnings
#from scipy.optimize import least_squares
from scipy.optimize import curve_fit
#from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr
#from scipy.stats import power_divergence, ks_2samp

###############################################################################
# FUNCTIONS //////////////////////////////////////////////////////////////////#
###############################################################################

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Calculate the FWHM and hence the extraction limits from a Moffat profile based on the distance from the central peak
# as a multiple of FWHM. A high resolution profile can also be produced if
# necessary by setting highres=True.

def extraction_limits(moffparams, width_multiplier, axesdict, highres=False):

    # If requested create a Moffat profile based on the input parameters that has 50 times the resolution of the
    # original data. Else, use the native resolution of the data
    if highres:
        r = np.linspace(0, axesdict['spataxislen'] - 0.02, num=axesdict['spataxislen'] * 50)
        moffprof = moffat(
            moffparams[0],
            moffparams[1],
            moffparams[2],
            moffparams[3],
            moffparams[4],
            moffparams[5],
            r)
        fwhm = (2 * moffparams[2] * (((2**(1 / moffparams[3])) - 1)**0.5)) * 50

    else:
        r = np.linspace(0, axesdict['spataxislen'] - 1, num=axesdict['spataxislen'])
        moffprof = moffat(moffparams[0], moffparams[1], moffparams[2], moffparams[3], 0, 0, r)
        fwhm = 2 * moffparams[2] * (((2**(1 / moffparams[3])) - 1)**0.5)

    # Locate the peak of the Moffat profile
    indexmax = np.argmax(moffprof)

    # Respectively define the upper and lower extraction limits at a distance above and below the peak of the Moffat
    # profile that equals the FWHM of the Moffat profile times a multiplying factor.
    lower_extraction_limit = indexmax - (width_multiplier * fwhm)
    upper_extraction_limit = indexmax + (width_multiplier * fwhm)

    return lower_extraction_limit, upper_extraction_limit, fwhm, moffparams[1]


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Linearly extrapolate the extraction limits at the ends of the 2D spectrum.
def extrap_extraction_lims(extlims, dispaxislen, shortend, longend):
    short_extrap_grad1 = extrap_grad(extlims[0], [0, 150, 300])
    short_extrap_grad2 = extrap_grad(extlims[1], [0, 150, 300])
    long_extrap_grad1 = extrap_grad(extlims[0], [-300, -150, -1])
    long_extrap_grad2 = extrap_grad(extlims[1], [-300, -150, -1])

    short_extrap_lim1 = extlims[0][0] - (short_extrap_grad1 * shortend)
    short_extrap_lim2 = extlims[1][0] - (short_extrap_grad2 * shortend)
    long_extrap_lim1 = extlims[0][-1] + (long_extrap_grad1 * (dispaxislen - longend))
    long_extrap_lim2 = extlims[1][-1] + (long_extrap_grad2 * (dispaxislen - longend))

    return short_extrap_lim1, short_extrap_lim2, long_extrap_lim1, long_extrap_lim2


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Given a range of data and limits to define a region of that data, calculate the region's gradient.
def extrap_grad(intextlims, median_lims):
    median_of_y_points_x_to_y = np.median(intextlims[median_lims[0]:median_lims[1]])
    median_of_y_points_y_to_z = np.median(intextlims[median_lims[1]:median_lims[2]])
    median_of_x_points_x_to_y = np.median([median_lims[0], median_lims[1]])
    median_of_x_points_y_to_z = np.median([median_lims[1], median_lims[2]])

    return (median_of_y_points_x_to_y - median_of_y_points_y_to_z) / \
        (median_of_x_points_x_to_y - median_of_x_points_y_to_z)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Define the bins of data over which Moffat profiles will be fitted. Each bin is defined such that when summed it will
# have a given signal to noise. So lower signal to noise regions will have larger bins. This function works from the
# centre of the input 2D spectrum outward to the ends in order to ensure a good start for the binning process.
# Flagged bad pixels will be replaced with appropriate values using a method described by the ESO X-Shooter
# pipeline. Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211
def get_bins(fdict, slow, shigh, dispaxislen, params, sky=False):

    # Take S/N threshold (minSNR) and minimum number of columns per dispersion bin (mincols)
    if params['-SUBTRACT_SKY'] and sky == True:
        minSNR = params['-SKY_SNR_BIN_LIM']
    elif params['-SUBTRACT_SKY'] and sky == False:
        minSNR = params['-SNR_BIN_LIM']
    else:
        minSNR = params['-SNR_BIN_LIM']

    # Minimum columns for a bin
    mincols = params['-COL_BIN_LIM']

    # Start x at the central pixel column of the dispersion axis
    x = int(dispaxislen / 2)
    width = 0
    binlocations = []

    sys.stdout.write(' >>> Determining spectrum localisation bins on dispersion axis.\n')
    sys.stdout.write('     User-defined S/N threshold = ' + str(minSNR) + '\n')

    # For this function, filter out any NaN/inf data before determining bins
    # and ensure the errors are positive.
    fdict['data'], fdict['errs'] = filter_data(fdict['data'], np.abs(fdict['errs']))

    ### Shorter Wavelengths ###

    # Start at the centre of the dispersion axis and start binning towards the
    # short wavelength end of the spectrum.
    while x - width > 0:
        snrestimate = 0.

        # If the S/N of the current bin has not yet reached the user-defined threshold (minSNR), add one more pixel
        # column to the bin.
        while snrestimate <= minSNR:
            width += 1

            # Stop the loop once the short wavelength end of the spectrum os reached.
            if x - width < 0:
                width = int(0 + x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin, continue to the next
            # iteration and add another pixel column to the bin.
            shortrows = len(list(filter(lambda x: x <= mincols, np.sum(
                fdict['qual'][:, int(x - width):int(x)], axis=1))))
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.

            # For simplicity, slice the data first.
            data_slice = fdict['data'][:, int(x - width):int(x)]
            errs_slice = fdict['errs'][:, int(x - width):int(x)]

            # Now get the collapsed profile
            bindatacol = np.nansum(data_slice, axis=1)
            binerrscol = np.sqrt(np.nansum(np.power(errs_slice, 2), axis=1))

            # If signal = 0 when slow <0 and shigh > bindatacol.shape[0]. Rectify this.
            if slow < 0:
                slow = 0
            if shigh > bindatacol.shape[0]:
                shigh = bindatacol.shape[0]

            # Now get the signal and noise
            signal = np.nansum(bindatacol[int(slow):int(shigh)])
            rssnoise = np.sqrt(np.nansum(np.power(binerrscol[int(slow):int(shigh)], 2)))

            # Should something fail eitherway, skip this column and move on.
            if ((rssnoise == 0) or (np.isnan(rssnoise)) or (np.isnan(signal))):
                snrestimate = 1
                continue

            # Estimate the S/N
            snrestimate = signal / rssnoise

        if params['-REPLACE_CRBP']:

            # Replace bad pixels, and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211

            # Here a median spatial profile for the current bin is determined by bootstrapping the good pixels
            # in each spatial pixel row with repeats to estimate a distribution of median values for each value in the
            # median spatial distribution for the bin. The mean of each of these distributions is then taken to be the
            # value of that median spatial pixel. The standard error of each of these distributions becomes the error
            # of the flux.
            meddatacol = []
            errmeddatacol = []
            for j in range(np.shape(fdict['data'])[0]):
                bootsamp = np.random.choice(fdict['data'][j, int(
                    x - width):int(x)], size=(int(width), 100), replace=True)
                meddatacol.append(np.nanmedian(bootsamp))
                meanstddist = np.nanstd(bootsamp, axis=0)
                meanstd = np.mean(meanstddist)
                errmeddatacol.append(meanstd / np.sqrt(width - 1))

            meddatacol = np.array(meddatacol)
            errmeddatacol = np.array(errmeddatacol)

            # Each pixel in the bin that was flagged to have be bad pixels or be contaminated by CRs are replaced by
            # the value at the corresponding spatial coordinate in the median spatial profile.
            for i in range(int(width)):
                if 0. in fdict['data'][:, int(x - i)]:
                    nocr = np.where(fdict['data'][:, int(x - i)] != 0.)

                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    datanocr = fdict['data'][:, int(x - i)][nocr]
                    medcolnocr = meddatacol[nocr]

                    # Add an exception to handle diving by zero / invalid value.
                    if ((np.sum(medcolnocr[slow:shigh]) == 0) or (
                            np.max(medcolnocr[slow:shigh] == 0))):
                        scale = 1
                    else:
                        scale = np.sum(datanocr[slow:shigh]) / np.sum(medcolnocr[slow:shigh])
                        if scale < 0:
                            scale = np.max(datanocr[slow:shigh]) / np.max(medcolnocr[slow:shigh])
                    meddatacol *= scale
                    errmeddatacol *= scale

                    # Substittue scaled median spatial profile values into the pixels where CRs or bad pixels have been
                    # identified.
                    multiplier = np.ones(len(meddatacol))
                    multiplier[nocr] *= 0
                    crsubstitute = multiplier * meddatacol
                    fdict['data'][:, int(x - i)] += crsubstitute
                    crerrsubstitute = multiplier * errmeddatacol
                    fdict['errs'][:, int(x - i)] += crerrsubstitute

        binlocations.append([int(x - width), int(x), snrestimate])

        # Shift to the shorter wavelength and reset the bin width. Continue the while loop.
        x -= width
        width = 0

    ### Longer Wavelengths ###

    # Repeat the same process as above, starting at the centre of the dispersion axis, but moving outward toward the
    # longer wavelength end of the 2D spectrum.
    x = int(dispaxislen / 2) + 1
    while x + width < dispaxislen:
        snrestimate = 0.

        # If the S/N of the current bin has not yet reached the user-defined threshold (minSNR), add one more pixel
        # column to the bin.
        while snrestimate <= minSNR:
            width += 1

            # Stop the loop once the long wavelength end of the spectrum is reached.
            if x + width > dispaxislen:
                width = int(dispaxislen - x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin, continue to the next
            # iteration and add another pixel column to the bin.
            shortrows = len(list(filter(lambda x: x <= mincols, np.sum(
                fdict['qual'][:, int(x):int(x + width)], axis=1))))
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.

            # For simplicity, slice the data first.
            data_slice = fdict['data'][:, int(x):int(x + width)]
            errs_slice = fdict['errs'][:, int(x):int(x + width)]

            # Now get the collapsed profile
            bindatacol = np.nansum(data_slice, axis=1)
            binerrscol = np.sqrt(np.nansum(np.power(errs_slice, 2), axis=1))

            # signal = 0 when slow <0 and shigh > bindatacol.shape[0]. Rectify this.
            if slow < 0:
                slow = 0
            if shigh > bindatacol.shape[0]:
                shigh = bindatacol.shape[0]

            # Now get the signal and noise
            signal = np.nansum(bindatacol[int(slow):int(shigh)])
            rssnoise = np.sqrt(np.nansum(np.power(binerrscol[int(slow):int(shigh)], 2)))

            # Should something fail eitherway, skip this column and move on.
            if ((rssnoise == 0) or (np.isnan(rssnoise)) or (np.isnan(signal))):
                snrestimate = 1
                continue

            # Estimate the S/N
            snrestimate = signal / rssnoise

        if params['-REPLACE_CRBP']:
            # Replace bad pixels, and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211
            meddatacol = []
            errmeddatacol = []

            # Here a median spatial profile for the current bin is determined by bootstrapping the good pixels
            # in each spatial pixel row with repeats to estimate a distribution of median values for each value in the
            # median spatial distribution for the bin. The mean of each of these distributions is then taken to be the
            # value of that median spatial pixel. The standard error of each of these distributions becomes the error
            # of the flux.
            for j in range(np.shape(fdict['data'])[0]):
                bootsamp = np.random.choice(fdict['data'][j, int(x):int(
                    x + width)], size=(int(width), 100), replace=True)
                meddatacol.append(np.nanmedian(bootsamp))
                meanstddist = np.nanstd(bootsamp, axis=0)
                meanstd = np.mean(meanstddist)
                errmeddatacol.append(meanstd / np.sqrt(width - 1))

            meddatacol = np.array(meddatacol)
            errmeddatacol = np.array(errmeddatacol)

            # Each pixel in the bin that was flagged to have be bad pixels or be contaminated by CRs are replaced by
            # the value at the corresponding spatial coordinate in the median spatial profile.
            for i in range(int(width)):
                if 0. in fdict['data'][:, int(x + i)]:
                    nocr = np.where(fdict['data'][:, int(x + i)] != 0.)

                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    datanocr = fdict['data'][:, int(x + i)][nocr]
                    medcolnocr = meddatacol[nocr]

                    # Add an exception to handle diving by zero / invalid value
                    if ((np.sum(medcolnocr[slow:shigh]) == 0) or (
                            np.max(medcolnocr[slow:shigh] == 0))):
                        scale = 1
                    else:
                        scale = np.sum(datanocr[slow:shigh]) / np.sum(medcolnocr[slow:shigh])
                        if scale < 0:
                            scale = np.max(datanocr[slow:shigh]) / np.max(medcolnocr[slow:shigh])
                    meddatacol *= scale
                    errmeddatacol *= scale

                    # Substitue scaled median spatial profile values into the pixels where CRs or bad pixels have been
                    # identified.
                    multiplier = np.ones(len(meddatacol))
                    multiplier[nocr] *= 0
                    crsubstitute = multiplier * meddatacol
                    fdict['data'][:, int(x + i)] += crsubstitute
                    crerrsubstitute = multiplier * errmeddatacol
                    fdict['errs'][:, int(x + i)] += crerrsubstitute

        binlocations.append([int(x), int(x + width), snrestimate])

        # Shift to the shorter wavelength and reset the bin width. Continue the while loop.
        x += width
        width = 0

    # Sort out the bins into the correct order on the dispersion axis.
    binlocations.sort(key=lambda x: x[0])

    return binlocations, fdict


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Print and plot the output of get_bins
def get_bins_output(binparams, params, lowext, highext, data2D, headparams, axdict):
    sys.stdout.write(' >>> ' + str(len(binparams)) +
                     ' spectrum localisation bins determined on dispersion axis.\n')

    # DIAGNOSTICS - Plot boundary locations of localisation bins on the dispersion axis.
    if params['-DIAG_PLOT_BIN_LOC']:
        drawlines = []
        for bin in binparams:
            binlineloc = np.where(np.logical_and(
                axdict['saxis'] > lowext - ((highext - lowext) * 0.2), axdict['saxis'] < highext + ((highext - lowext) * 0.2)))
            drawlines.append(np.ones(len(axdict['saxis'][binlineloc]))
                             * bin[0] + axdict['wavstart'])
            drawlines.append(axdict['saxis'][binlineloc] + axdict['imgstart'])
            drawlines.append(np.ones(len(axdict['saxis'][binlineloc]))
                             * bin[1] + axdict['wavstart'])
            drawlines.append(axdict['saxis'][binlineloc] + axdict['imgstart'])

        show_img(data2D, axdict, headparams, drawlines,
                 '2D Spectrum with Boundaries of Localisation Bins')

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Takes an input of extraction limits from the fitting of the binned data and interpolates the limits over the unbinned
# data. Limits are also linearly extrapolated towards the ends of the spectral range.
def interpolate_extraction_lims(extractionlims, dispaxislen, interpkind):
    # If the 2D spectrum was so faint that only 1 dispersion bin could be determined, set the extraction limits
    # across the unbinned spectrum to be a simple linear aperture that has the same extraction limits as was
    # determined for that one bin.
    if len(extractionlims[0]) == 1:
        finalextlims = [
            np.repeat(
                extractionlims[1][0], dispaxislen), np.repeat(
                extractionlims[2][0], dispaxislen)]

    # Otherwise, interpolate the extraction limits form the bins across the unbinned wavelength axis. Also
    # extrapolate the extraction limits in a linear fashion at the ends of the wavelength axis of the 2D spectrum
    # so that the full wavelength axis is covered.
    else:
        interpextract_1 = interp.interp1d(extractionlims[0], extractionlims[1], kind=interpkind)
        interpextract_2 = interp.interp1d(extractionlims[0], extractionlims[2], kind=interpkind)

        intermextlims = [interpextract_1(np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))).astype(int),
                         interpextract_2(np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))).astype(int)]

        shortextraplim1, shortextraplim2, longextraplim1, longextraplim2 = extrap_extraction_lims(
            intermextlims, dispaxislen, extractionlims[0][0], extractionlims[0][-1])

        extlim1 = np.insert(intermextlims[0], 0, shortextraplim1)
        extlim1 = np.append(extlim1, longextraplim1)

        extlim2 = np.insert(intermextlims[1], 0, shortextraplim2)
        extlim2 = np.append(extlim2, longextraplim2)

        nextextlims = [extlim1, extlim2]
        nextxaxis = np.insert(np.arange(extractionlims[0][0], extractionlims[0][-1]), 0, 0)
        nextxaxis = np.append(nextxaxis, dispaxislen)

        interpnextlims1 = interp.interp1d(
            nextxaxis,
            nextextlims[0],
            kind='linear',
            fill_value='extrapolate')
        interpnextlims2 = interp.interp1d(
            nextxaxis,
            nextextlims[1],
            kind='linear',
            fill_value='extrapolate')

        finalextlims = [
            interpnextlims1(
                np.array(
                    range(dispaxislen))), interpnextlims2(
                np.array(
                    range(dispaxislen)))]

    return finalextlims


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Returns a wavelength axis array using a start wavelength, the wavelength increment and the number of values along the
# axis required.

def make_wav_axis(start, increment, length):
    return np.arange(start=start, step=increment, stop=(length * increment) + start)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function detects cosmic rays in the 2D image and masks the pixels that contain them, by separating the spectrum
# region from the two background regions. Each region is sigma-clipped separately due to the different distribution of
# fluxes present in the background and the vicinity of the 2D spectrum.
def mask_cosmic_rays(data, lext, hext, multiplier=4.0):

    # Divide 2D spectrum into spectrum and background sections. This is slicing
    # data [below, centre, above] the apeture
    sections = [data[:int(lext), :], data[int(lext):int(hext), :], data[int(hext):, :]]

    # Requires similar correction to get_bins() to account for some limits
    # being outside the image (eg. individual ABBA frames reduced in STARE
    # mode of the X-Shooter ESO pipeline).
    if (lext < 0):
        lext = 0
        sections = [data[int(lext):int(hext), :], data[int(hext):, :]]
    if (hext > data.shape[0]):
        hext = data.shape[0]
        sections = [data[:int(lext), :], data[int(lext):int(hext), :]]

    # Sequentially sigma-clip each section, and create a mask marking the
    # locations of pixels to be kept.
    masks = []
    for sect in sections:

        datamed = np.nanmedian(sect)
        datastd = np.nanstd(sect)

        upperlim = datamed + (multiplier * datastd)
        lowerlim = datamed - (multiplier * datastd)

        sectmask2D = (sect < upperlim) & (sect > lowerlim)
        masks.append(sectmask2D)

    # Assemble masks for each section into a full mask for the input 2D spectrum.
    if len(masks) == 3:
        mask2D = np.vstack((masks[0], masks[1], masks[2]))
    else:
        mask2D = np.vstack((masks[0], masks[1]))

    return mask2D


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

# Creates moffat profile added a linear sloped background based on input parameters
def moffat(amp, c, alpha, beta, bglevel, bggrad, datarange):
    moff = amp * ((1 + (((datarange - c) * (datarange - c)) / (alpha * alpha))) ** -beta)

    return moff + bglevel + (datarange * bggrad)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
'''
(Deprecated)
# Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a
# Levenberg-Marquardt least squares method. Returns the best fit parameters of the Moffat function
def moffat_least_squares(r, col, seeing, pixres):
    # Set up initial conditions for the least squares fit.
    # amplitude, centre, alpha, beta, background gradient, background level
    # alpha value calculated by substituting seeing (FWHM) converted to pixels and beta = 4.765 into the equation:
    #               alpha = FWHM / (2 * ((2 ^ (1 / beta)) - 1) ^ 0.5)
    # Initial beta value comes from optimal value from atmospheric turbulence theory as described in
    # Trujillo I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [np.nanmedian(np.sort(col)[-3:]), np.argmax(col), seeing / 0.7914, 4.765, 0.,
          np.median(np.concatenate((col[:5], col[-5:])))]

    # Run the least squares fit.
    res_lsq = least_squares(moffat_resid,
                            x0,
                            bounds=([0., np.argmax(col) - 1., 0., 0., 0., -np.inf], [np.inf,
                                                                                     np.argmax(col) + 1, (5 * seeing / 0.7914), 5., np.inf, np.inf]),
                            args=(r, col),
                            method='trf',
                            ftol=1E-12
                            )

    return [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2], res_lsq.x[3], res_lsq.x[4], res_lsq.x[5]]
'''
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
'''
# Calculates residuals of fitted moffat profile and the data for the Levenberg Marquardt least squares method.
# A = x[0]
# c = x[1]
# alpha = x[2]
# beta = x[3]
# B = x[4]
# m = x[5]
# @jit

(Deprecated)
def moffat_resid(x, datarange, data):
    moff = x[0] * ((1 + ((datarange - x[1]) * (datarange - x[1])) / (x[2] * x[2])) ** -x[3])

    return moff + x[4] + (datarange * x[5]) - data
'''
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Plot the spatial profile of a collapsed spectrum or a collapsed bin therein, and plot the Moffat function fitted to
# the data on top.


def plot_fitted_spatial_profile(spataxis, bindata, hiresspataxis,
                                binmoffparams, imgstart, headparams):
    plt.figure(figsize=(7, 4.5))
    plt.plot(hiresspataxis + imgstart,
             moffat(
                 binmoffparams[0],
                 binmoffparams[1],
                 binmoffparams[2],
                 binmoffparams[3],
                 binmoffparams[4],
                 binmoffparams[5],
                 hiresspataxis),
             color='r',
             linewidth=3,
             label='Fitted Moffat Profile'
             )

    plt.plot(spataxis + imgstart, bindata, color='k', label='Spatial Profile')
    plt.grid(linestyle='dashed', color='gray')
    plt.legend()
    plt.title('Spectrum Spatial Profile and Fitted Moffat Profile')
    plt.xlabel('Spatial Axis, Pixels')
    plt.ylabel('Median Flux, ' + headparams['fluxunit'])
    plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function prints, to the terminal, a set parameters describing a fitted Moffat function.
def printmoffparams(moffparams, imgstart, datascale):
    sys.stdout.write(' >>> Fitted Moffat function parameters:\n')
    sys.stdout.write('         A = ' + str(moffparams[0]) + '\n')
    sys.stdout.write('         c = ' + str(moffparams[1] + imgstart) + '\n')
    sys.stdout.write('     alpha = ' + str(moffparams[2]) + '\n')
    sys.stdout.write('      beta = ' + str(moffparams[3]) + '\n')
    sys.stdout.write('         B = ' + str(moffparams[4]) + '\n')
    sys.stdout.write('         m = ' + str(moffparams[5]) + '\n\n')
    sys.stdout.write(' >>> Profile scaling factor used for fitting: ' + str(datascale) + '\n')
    sys.stdout.write(' >>> Plot of median spatial profile presents the orginal unscaled profile.\n')

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Takes an input image and line data to be drawn on that image and creates
# a figure to be shown on screen.
def show_img(data2D, axdict, headparams, drawlines, title):

    # Catch and suppress the following UserWarning:
    # UserWarning: Warning: 'partition' will ignore the 'mask' of the MaskedArray.
    #     a.partition(kth, axis=axis, kind=kind, order=order)
    # Numpy masked array used only when displaying the spectrograms.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        warnings.warn("partition", UserWarning)

        power = int(np.floor(np.log10(np.nanmean(data2D))))
        data2D = copy.deepcopy(data2D) / 10**power

        figwidth = 10.
        fig = plt.figure(figsize=(figwidth, figwidth / 1.9))
        gs = gridspec.GridSpec(18, 33)
        ax = plt.subplot(gs[:, :32])
        colax = plt.subplot(gs[1:, 32])
        masked_data2D = np.ma.masked_where(data2D == 0, data2D)
        cmap = matplotlib.cm.afmhot
        cmap.set_bad(color='red')
        s = ax.imshow(masked_data2D,
                      aspect='auto',
                      vmin=0,
                      vmax=np.nanmedian(masked_data2D) + (0.5 * np.nanstd(masked_data2D)),
                      origin='lower',
                      cmap=cmap,
                      extent=[axdict['wavstart'],
                              axdict['wavstart'] + len(axdict['waxis']),
                              axdict['saxis'][0] + axdict['imgstart'],
                              axdict['saxis'][-1] + axdict['imgstart']]
                      )

        for i in range(int(len(drawlines) / 2)):
            ax.plot(drawlines[i * 2], drawlines[(i * 2) + 1], color='white')
        cbar = fig.colorbar(s, cax=colax)
        cbar.ax.yaxis.set_offset_position('left')
        cbar.ax.set_ylabel('Pixel Flux, x10^' + str(power) + ' ' + headparams['fluxunit'])
        ax2 = ax.twiny()
        ax2.plot(axdict['waxis'], data2D[0, :], alpha=0)
        ax2.set_xlim(axdict['waxis'][0], axdict['waxis'][-1])
        ax2.set_xlabel('Wavelength, ' + headparams['wavunit'])
        ax.set_ylim(axdict['saxis'][0] + axdict['imgstart'],
                    axdict['saxis'][-1] + axdict['imgstart'])
        ax.set_ylabel('Spatial Axis, Pixels')
        ax.set_xlim(axdict['wavstart'], axdict['wavstart'] + len(axdict['waxis']))

        ax.set_xlabel('Dispersion Axis, Pixels')

        plt.title(title, y=1.095)

        # Add interactive scaling bar to figures if the 2D spectrum isn't flux calibrated.
        # Again for some reason teeny tiny numbers cause things to break.
        fig.subplots_adjust(bottom=0.2)
        axvmin = plt.axes([0.1, 0.05, 0.8, 0.03])
        axvmax = plt.axes([0.1, 0.01, 0.8, 0.03])
        smin = Slider(axvmin, 'LowCut', 0, np.nanmax(masked_data2D) - 1,
                      valinit=0., valstep=0.001 * np.nanmax(masked_data2D))
        smax = Slider(axvmax, 'HighCut', 1., np.nanmax(masked_data2D), valinit=np.nanmedian(
            masked_data2D) + (0.5 * np.nanstd(masked_data2D)), valstep=0.001 * np.nanmax(masked_data2D))

        def update(val):
            vmax = smax.val
            smin.valmax = vmax - 1
            vmin = smin.val
            smax.valmin = vmin + 1

            s.set_clim(vmin, vmax)
            cbar.ax.set_ylabel('Pixel Flux, x10^' + str(power) + ' ' + headparams['fluxunit'])
            fig.canvas.draw_idle()

        smin.on_changed(update)
        smax.on_changed(update)

        plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Subtracts the sky background from the 2D image by defining bg regions using limits input to the function-
#
# (and then fitting a profile to the background column by column while masking cosmic rays. The fitted linear for
# each column is subtracted from the full column to produce a background
# subtracted 2D image) - not sure if valid anymore.

def subtract_sky(bglowext, bghighext, fdict, axdict):

    fdict['data'] = fdict['data'].T
    fdict['errs'] = fdict['errs'].T
    medsky = []

    # Makes sure the limits are within the image.
    for ii in range(len(bglowext)):
        if bglowext[ii] < 0:
            bglowext[ii] = 0
        if bghighext[ii] > axdict['saxis'][-1]:
            bghighext[ii] = axdict['saxis'][-1]

        datacol = fdict['data'][ii]
        colrange = range(len(datacol))
        skypix = datacol[np.where(np.logical_or(colrange < bglowext[ii], colrange > bghighext[ii]))]

        bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
        skysamp = np.nanmedian(bootsky, axis=0)
        skylevel = np.nanmean(skysamp)
        medsky.append(skylevel)

        skyerr = np.std(skysamp) / (99**0.5)
        fdict['data'][ii] -= skylevel
        fdict['errs'][ii] = ((fdict['errs'][ii]**2) + (skyerr**2))**0.5

    medsky = np.array(medsky)
    fdict['skymod'] = np.tile(medsky, (np.shape(fdict['data'])[1], 1))
    fdict['data'] = fdict['data'].T
    fdict['errs'] = fdict['errs'].T
    chipgaps = np.where(np.median(fdict['data'], axis=0) == 0)
    fdict['data'][:, chipgaps[0]] += 1

    return fdict


###############################################################################
# DEV FUNCTIONS //////////////////////////////////////////////////////////////#
###############################################################################

def qual_to_molecfit_compat(qual2D):

    # This is a method which ensures that MOTES doesn’t ignore the pixels that
    # have a QUAL extension which is deemed to be interpolated or have a bad
    # B-Spline fit as a result of ESO pipeline.

    # Those pixels have a special number in QUAL2D above 0. However, since ESO
    # standards (and their instruments such as ESO Molecfit) determines any pixel
    # above 0 as 'bad', a large number of pixels will be not read by other ESO
    # instruments, especially when using B-spline data - preventing from
    # interpolation. In other words, we want to keep only the ’truly’ bad
    # pixels inside.

    # This is folollowed by converting qual2D to binary GME standards, to
    # ensure compatibility with the cosmic mask/removal code.

    # Input = Original qual2D
    # Output = Corrected qual2D

    # Replace interpolated values so that its NOT a 'bad' pixel.
    # Interpolated flux during standard extraction - not really bad pixels
    qual2D[qual2D == 4194304] = 1
    # Both below prevent NIR molecfit operation on BSPLINE stare subtracted frames, even on SAs
    # Pixel where BSPLINE model sky fit is inaccurate
    qual2D[qual2D == 8388608] = 1
    # Outliers of BSPLINE sky model fit - again, not really bad pixels
    qual2D[qual2D == 16777216] = 1

    # Convert X-Shooter qual frame to GME standard boolean.
    # ie. Good pixels = 1; Bad pixels = -1 and 0
    qual2D[qual2D > 0] *= - 1
    qual2D[qual2D == 0] = 1
    qual2D[qual2D < 0] = 0

    return qual2D


def filter_data(data2D, errs2D):

    # This method takes in the data2D and errs2D and outputs frames where any
    # NaN or Inf values are 0.0 (this can be later filtered to a median of the column)

    # This is used before fitting profiles, determining bins before sky
    # subtraction and extraction procedures(which ensures that the S / N in
    # the bin is numerical), as well as before optimal extraction which rely
    # on formal errors to be numerical. Ultimately, this method makes sure the
    # final spectrum is fully numerical.

    # Input = Original data2D, errs2D
    # Output = Corrected data2D, errs2D

    # Zero out any NaNs and Infs from the error frame and corresponding data
    errs2D[np.logical_or(np.isnan(errs2D), np.isinf(errs2D))] = 0.0
    data2D[np.logical_or(np.isnan(errs2D), np.isinf(errs2D))] = 0.0
    # Likewise, ensure that data is numerical and  corresponding errors.
    data2D[np.logical_or(np.isnan(data2D), np.isinf(data2D))] = 0.0
    errs2D[np.logical_or(np.isnan(data2D), np.isinf(data2D))] = 0.0

    return data2D, errs2D


def optimal_extraction(data2D, errs2D, extractionlimits,
                       binparameters, datascale=10**18, seeing=1.0):

    # Perform optimal extraction using a modified version of Horne (1986) where S=0, G=0 and errors are not 'revised' since we already have the 'variance frame' (and we are not working with counts).
    # Ideally, this extraction reduces errors with the weighting and conserves the flux (when compared to the standard extraction)
    # This subroutine uses an analytic Moffat profile with a residual sky constant (post subtraction), instead of polynominals as in Horne (1986)
    # Furthermore, profiling takes place bin by bin, accounting for spatial
    # profile variations across the dispersion axis; extraction takes place
    # column by column within the bin limits and within the extraction limits previously calculated

    # Input = data2D, errs2D
    # Input = binparameters (contains the bin limits across the dispersion axis, to enable slicing the data across dispersion axis)
    # Input = finalextractionlimits (contains limits at each dispersion pixel)
    # Input = datascale and seeing
    # Outputs = optdata1D, opterrs1D, stddata1D, stderrs1D

    # Filter any NaNs and Inf for data/errs AND ensure the errors are positive for this extraction.
    data2D, errs2D = filter_data(data2D, np.abs(errs2D))

    # initiate the arrays to store the extracted spectra
    optdata1D = np.empty((0))
    opterrs1D = np.empty((0))
    stddata1D = np.empty((0))
    stderrs1D = np.empty((0))

    # See below
    bin_index1_list = np.zeros((1))
    bin_index0_list = np.zeros((1))
    centre_corrected = False

    # First loop (bin by bin)
    for i, bin in zip(range(len(binparameters)), binparameters):

        # In GME, the dispersion bin limits are determined in two parts.
        # From centre towards lower and then through higher wavelengths.
        # However, this results in 'a gap' between the bins near the centre, ommited by by this
        # extraction.
        # eg. [117, 125, 114.85173647422386] and [126, 136, 121.81477181579963].
        # Pixel at 125 will NOT the extracted in the first bin. Usually that is done in the next bin
        # Thus the extracted spectra is always one shorter of full wavaxis
        # I correct this by creating arrays to store bin[1] bin[0]. When the
        # bin[1] of the previous bin does not equal bin[0] of present, do
        # bin[0] -1 to include the missing pixel column from the previous bin.
        # The centre corrected bool

        # 22/11/20 - Added a correction to NOT include the very first bin by
        # mistake (second if statement)
        if (centre_corrected == False):
            if ((bin_index1_list[i - 1] != bin_index0_list[i]) and not (bin[0] <= 0)):
                bin[0] -= 1
                centre_corrected = True
            else:
                pass
        else:
            pass

        bin_index0_list = np.append(bin_index0_list, bin[0])
        bin_index1_list = np.append(bin_index1_list, bin[1])

        # Assume that previously derived binparameters and the finalextractionlimits are correct.
        # Define bin limits in the dispersion direction, in pixels
        # Take all values of extractionlimits within the bin and get their median (should be the same value since these limits are taken in as constants for each bin)
        # Ensure they are within the image if the spatial profile is at the borders.
        lower_limit = int(np.floor(np.median(extractionlimits[0][bin[0]:bin[1]])))
        if (lower_limit < 0):
            lower_limit = 0
        upper_limit = int(np.ceil(np.median(extractionlimits[1][bin[0]:bin[1]])))
        if (upper_limit > len(data2D[0])):
            upper_limit = len(data2D[0])

        # Slice the bin on the spatial AND dispersion axis. That part will be extracted.
        spectra_slice_data = data2D[lower_limit:upper_limit, bin[0]:bin[1]]
        spectra_slice_err = errs2D[lower_limit:upper_limit, bin[0]:bin[1]]

        # Collapse the data in the bin in dispersion direction by taking a sum
        # (and errors in quadrature, assuming they are 1st standard dev as in the
        # harvester functions)
        spectra_slice_data_sum = np.nansum(spectra_slice_data, axis=1)
        spectra_slice_errs_sum = np.sqrt(np.nansum(np.power(spectra_slice_err, 2), axis=1))

        # The profile needs to be refittted to satisfy both normality and positivity.
        # Firstly, the spatial axis for the profile
        spatial_axis = np.arange(0, len(spectra_slice_data_sum))

        # For simplicity, assume a Moffat profile with some small residual background level and without background gradient at this stage or reduction. (ie. abba=True)
        # Tests have shown that in NIR spectra of XShooter this works better for
        # now and its more robust against peaks that may occur beside the main
        # outliers in faint targets that may show as secondary peaks.
        profile_parameters, _ = moffat_weighted_curvefit(spatial_axis,
                                                         spectra_slice_data_sum * datascale,
                                                         seeing,
                                                         spectra_slice_errs_sum * datascale,
                                                         abba=True)

        # Descale the amplitude parameter and background
        profile_parameters[0] /= datascale
        profile_parameters[4] /= datascale
        profile_parameters[5] /= datascale
        # Define the profile values on the spatial axis
        profile = GME_Moffat1D_curvefit(spatial_axis, *profile_parameters)
        # Enforce positivity as in the Hornes 1986 paper - Zero out any negative/inf/none values on the profile.
        # This would be biased under normal circumstances, however, the profile is there to scale.
        profile[np.logical_or(profile < 0, np.isinf(profile))] = 0.0
        profile[profile is None] = 0.0
        # Normalise the bin optimal extraction profile using actual sum of flux in the bin
        normalised_profile = profile / np.nansum(profile)

        '''
        # Plot the profiles to check
        plt.title('Optimal Extraction Profile')
        plt.plot(
            spatial_axis,
            profile,
            color='green')
        plt.plot(spatial_axis, spectra_slice_data_sum, 'k.')
        plt.axhline(y=0, color='orange')
        plt.tight_layout()
        plt.show()
        '''

        # Second loop (column by column; extraction with the formula using the profile)
        for k in range(np.shape(spectra_slice_data)[1]):

            y = spectra_slice_data[:, k]
            y_err = spectra_slice_err[:, k]
            x = spatial_axis

            # Ensure that any zeros in errors are the median of the y_err array
            # to get some non-zero estimate for the optimal extraction.
            y_err[np.where(y_err == 0.0)[0]] = np.nanmedian(y_err)
            y_err[np.where(y_err is None)[0]] = np.nanmedian(y_err)

            # Convert the input error (one gaussian sigma) to variance. The variance is sigma**2.
            y_err = np.power(y_err, 2)

            # Perform standard extraction in the single column. The result error is one standard
            # deviation (1 sigma).
            standard_extraction_flux = np.sum(y)
            standard_extraction_err = np.sqrt(np.sum(y_err))

            # Append to 1D arrays.
            stddata1D = np.append(stddata1D, standard_extraction_flux)
            stderrs1D = np.append(stderrs1D, standard_extraction_err)

            # Step 5 of Horne 1986 - Cosmic ray masking with a conservative 5 sigmas.
            value_term = np.power(y - np.multiply(normalised_profile, standard_extraction_flux), 2)
            sigma_term = np.multiply(y_err, 25)
            cosmic_mask = np.invert([value_term > sigma_term])[0]
            # With invert, it gives False (ie. 0) if the above condition is satisfied which
            # means to mask the rays

            # Add an exepction to handle scenarious where the mask is all False If that happens, ignore the
            # cosmic ray correction.
            if not np.any(cosmic_mask):
                cosmic_mask = np.full(shape=cosmic_mask.shape, fill_value=True)

            # Extract the optimal spectrum (Step 8 of the algorithm)
            # Page 4 of Horne 1986: These formulae are equivalent to determining the
            # OPT spectrum by scaling a known spatial profile P, to fit the sky
            # subtracted data, whose variances are V.

            # Numerator and demoninator of the optimal extraction.
            optimal_numerator = np.sum(
                np.divide(
                    np.multiply(
                        normalised_profile,
                        y),
                    y_err)[cosmic_mask])
            optimal_denominator = np.sum(
                np.divide(
                    np.power(
                        normalised_profile,
                        2),
                    y_err)[cosmic_mask])

            # Finally, get the optimal values for flux and error. The resultant error is
            # one standard deviation (1 Gaussian sigma, again). If it fails
            # nontheless, use the std_extraction flux and errors.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                try:
                    optimal_extraction_flux = optimal_numerator / optimal_denominator
                except RuntimeWarning:
                    optdata1D = np.append(optdata1D, standard_extraction_flux)
                    opterrs1D = np.append(opterrs1D, standard_extraction_err)
                    continue

            # If the flux is obtained (denominator is not zero, get the corresponding error.
            optimal_extraction_err = np.sqrt(np.sum(normalised_profile[cosmic_mask]) /
                                             np.sum(np.divide(np.power(normalised_profile, 2), y_err)[cosmic_mask]))

            # Append the optimal extraction values to arrays
            optdata1D = np.append(optdata1D, optimal_extraction_flux)
            opterrs1D = np.append(opterrs1D, optimal_extraction_err)

    return optdata1D, opterrs1D, stddata1D, stderrs1D


def GME_Moffat1D_curvefit(x, amp, c, alpha, beta, bglevel, bggrad):
    # Standard MOTES fitting profile.
    return (amp * ((1 + (((x - c) * (x - c)) / (alpha * alpha))) ** -beta) + bglevel + (x * bggrad))


def GME_Moffat1D_profile(x, amp, c, alpha, beta, bglevel):
    # For getting the optimal extraction profile and NODDING mode. Assume some
    # residual constant background.
    return (amp * ((1 + (((x - c) * (x - c)) / (alpha * alpha))) ** -beta) + bglevel)


def moffat_weighted_curvefit(r, col, seeing, col_err, abba, maxfev=10000):

    # Non-linear weighted fitting using scipy's curve_fit. The data and errors
    # should be already scaled. Important: the input errors are scaled one sigma errors as output from the harvester.

    # Input = r (spatial axis), col (data column), seeing and col_err (err column), abba flag for double pass sky subtraction profile
    # Output = A list of moffat parameters ie. [A, c, alpha, beta, B, m]

    # Filter any NaNs and Inf in data/errs AND ensure the errors are positive for the fitting
    col, col_err = filter_data(col, np.abs(col_err))
    # Zero/None errors give an infinite weight, make them np.median of the column
    col_err[np.where(col_err == 0.0)[0]] = np.nanmedian(col_err)
    col_err[np.where(col_err is None)[0]] = np.nanmedian(col_err)

    # Initialise the initial values.
    # Changing to 1.05 x A as the upper limit fits the moffat more efficiently.
    if (not abba):
        x0 = [np.nanmax(col), r[np.argmax(col)], seeing, 4.675, 0.1,
              np.nanmedian(np.concatenate((col[:5], col[-5:])))]
        bounds_curvefit = ([0, x0[1] - 1, 0, 0, -np.inf, -np.inf], [(x0[0]) * 1.05,
                                                                    x0[1] + 1, 5 * (x0[2] / 0.7914), 10, np.inf, np.inf])
        moff_function = GME_Moffat1D_curvefit
    else:
        x0 = [np.nanmax(col), r[np.argmax(col)], seeing, 4.675, 0.1]
        bounds_curvefit = ([0, x0[1] - 1, 0, 0, 0], [(x0[0]) * 1.05,
                                                     x0[1] + 1, 5 * (x0[2] / 0.7914), 10, np.inf])
        moff_function = GME_Moffat1D_profile

    # Fit the profile with a Moffat, remember to scale the data.
    # Should something go wrong (eg. ValueError), try fitting without errors.
    try:
        popt, _ = curve_fit(
            f=moff_function,
            xdata=r,
            ydata=col,
            sigma=col_err,
            absolute_sigma=True,
            p0=x0,
            bounds=bounds_curvefit,
            maxfev=maxfev)
    except BaseException:
        popt, _ = curve_fit(
            f=moff_function,
            xdata=r,
            ydata=col,
            p0=x0,
            bounds=bounds_curvefit,
            maxfev=maxfev)

    # y is the expected spatial flux distribution as given by the fitted profile (say across the bin)
    # x is 'col' the actual data used to fit.
    y = moff_function(r, *popt)
    fitting_stats = pearsonr(x=col, y=y)
    moffatparams = popt.tolist()
    if abba:
        # Append the zero gradient to ensure compatibility with the rest of the software.
        moffatparams.append(0.0)

    return moffatparams, float(fitting_stats[1])
