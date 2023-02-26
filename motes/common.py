"""
common.py - Common functions for the MOTES pipeline.
"""

import sys
import warnings
import copy
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from matplotlib.widgets import Slider
from scipy.optimize import least_squares


def extraction_limits(moffparams, axesdict, width_multiplier=3.0):
    """Calculate the FWHM (i.e., the extraction limits) from a Moffat profile based on the distance from the central peak as a multiple of FWHM.

    Args:
        moffparams (list): list containing the parameters of the Moffat profile to be created and measured. The parameters are: [amplitude, location, scale, power].
        axesdict (dict): dictionary containing the axes of the data frame 2D spectrum. Here only the length of the spatial axis is retrieved.
        width_multiplier (float, optional): defines the distance from the center of the spatial profile at which to set the extraction limits, in multiples of the FWHM. Defaults to 3.0.

    Returns:
        lower_extraction_limit (numpy.float64) - the lower bound of the region to be extracted.
        upper_extraction_limit (numpy.float64) - the upper bound of the region to be extracted.
        fwhm (numpy.float64)                   - the Full Width at Half Maximum of the Moffat profile.
        moffparams[1] (numpy.float64)          - location of the center of the Moffat profile.
    """

    # Create a Moffat profile based on the input parameters.
    r = np.linspace(0, axesdict["spataxislen"] - 1, num=axesdict["spataxislen"])
    fwhm = 2 * moffparams[2] * (((2 ** (1 / moffparams[3])) - 1) ** 0.5)

    # Respectively define the upper and lower extraction limits at a distance above and below the peak of the Moffat
    # profile that equals the FWHM of the Moffat profile times a multiplying factor.
    lower_extraction_limit = moffparams[1] - (width_multiplier * fwhm)
    upper_extraction_limit = moffparams[1] + (width_multiplier * fwhm)

    return lower_extraction_limit, upper_extraction_limit, fwhm, moffparams[1]


def extrap_extraction_lims(extlims, dispaxislen, shortend, longend):
    """Linearly extrapolate the extraction limits at the ends of the 2D spectrum.

    Description:
    inearly extrapolate the extraction limits at the ends of the 2D spectrum.

    Args:
        extlims (list): A list containing the lower and upper extraction limits for each spatial pixel.
        dispaxislen (int): The length of the dispersion axis.
        shortend (int): The number of pixels at the short end of the dispersion axis to be excluded from the extraction.
        longend (int): The number of pixels at the long end of the dispersion axis to be excluded from the extraction.

    Returns:
    pars (dict): A dictionary containing the parameters read in from motesparams.txt.

    """
    short_extrap_grad1 = extrap_grad(extlims[0], [0, 150, 300])
    short_extrap_grad2 = extrap_grad(extlims[1], [0, 150, 300])
    long_extrap_grad1 = extrap_grad(extlims[0], [-300, -150, -1])
    long_extrap_grad2 = extrap_grad(extlims[1], [-300, -150, -1])

    short_extrap_lim1 = extlims[0][0] - (short_extrap_grad1 * shortend)
    short_extrap_lim2 = extlims[1][0] - (short_extrap_grad2 * shortend)
    long_extrap_lim1 = extlims[0][-1] + (
        long_extrap_grad1 * (dispaxislen - longend)
    )
    long_extrap_lim2 = extlims[1][-1] + (
        long_extrap_grad2 * (dispaxislen - longend)
    )

    return (
        short_extrap_lim1,
        short_extrap_lim2,
        long_extrap_lim1,
        long_extrap_lim2,
    )


def extrap_grad(intextlims, median_lims):
    """Calculate the gradient of a region of data.

    Description:
    Given a range of data and limits to define a region of that data, calculate the region's gradient.

    Args:
        intextlims (numpy.ndarray): A range of data.
        median_lims (list): A list containing the limits of the region of data to be used to calculate the gradient.

    Returns:
       gradient (int) = The gradient of the region of data.

    """
    median_of_y_points_x_to_y = np.median(
        intextlims[median_lims[0] : median_lims[1]]
    )
    median_of_y_points_y_to_z = np.median(
        intextlims[median_lims[1] : median_lims[2]]
    )
    median_of_x_points_x_to_y = np.median([median_lims[0], median_lims[1]])
    median_of_x_points_y_to_z = np.median([median_lims[1], median_lims[2]])

    gradient = (median_of_y_points_x_to_y - median_of_y_points_y_to_z) / (
        median_of_x_points_x_to_y - median_of_x_points_y_to_z
    )

    return gradient


def filter_data(data2D, errs2D):
    """Filters the data2D and errs2D to remove NaNs and Infs.

    Description:
    This function takes in the data2D and errs2D and outputs frames where any
    NaN or Inf values are 0.0 (this can be later filtered to a median of the column)

    This is used before the extraction procedures to ensure the S/N in the bin is numerical (required for optimal extraction).

    Args:
        data2D (numpy.ndarray): Original data2D
        errs2D (numpy.ndarray): Original errs2D

    Returns:
        data2D (numpy.ndarray): Filtered data2D
        errs2D (numpy.ndarray): Filtered errs2D
    """

    errs2D[np.isfinite(errs2D) == False] = 0.0
    data2D[np.isfinite(errs2D) == False] = 0.0
    data2D[np.isfinite(data2D) == False] = 0.0
    errs2D[np.isfinite(data2D) == False] = 0.0

    return data2D, errs2D


def get_bins(
    fdict, slow, shigh, dispaxislen, params, sky=False, replace_crbp=False
):
    """Define the bins of data over which Moffat profiles will be fitted.

    Description:

    Define the bins of data over which Moffat profiles will be fitted. Each bin is defined such that when summed it will have a given signal to noise (S/N). So lower S/N regions will have larger bins. This function works from the centre of the input 2D spectrum outward to the ends in order to ensure a good start for the binning process.

    Flagged bad pixels will be replaced with appropriate values using a method described by the ESO X-Shooter pipeline. Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211

    Args:
        fdict (dict): Dictionary containing all the data, error, and quality frames.
        slow (int): Lower spatial limit of the region where the S/N will be measured when defining the extent of a bin.
        shigh (int): Upper spatial limit of the region where the S/N will be measured when defining the extent of a bin.
        dispaxislen (int): Length of the dispersion axis.
        params (dict): Dictionary of parameters ready in from the motesparams.txt configuration file.
        sky (bool): Optional kwarg used to tell get_bins whether the sky has been subtracted yet or not, and to set the minSNR threshold accordingly.
        replace_crbp (bool): Optional kwarg; when True, get_bins will try to replace bad pixels with values estimated using the median spatial profile of the current bin.

    Returns:
        binlocations (list): A list containing the details for each bin determined by get_bins. The boundaries and S/N of each bin are recorded here.
        fdict (dict): Returns fdict with bad pixels replaced if replace_crbp=True. If replace_crbp=False, fdict is returned unchanged.
    """

    # Take S/N threshold (minSNR) and minimum number of columns per dispersion bin (mincols)
    if params["-SUBTRACT_SKY"] and sky == True:
        minSNR = params["-SKY_SNR_BIN_LIM"]
    elif params["-SUBTRACT_SKY"] and sky == False:
        minSNR = params["-SNR_BIN_LIM"]
    else:
        minSNR = params["-SNR_BIN_LIM"]

    # Minimum columns for a bin
    mincols = params["-COL_BIN_LIM"]

    # Start x at the central pixel column of the dispersion axis
    x = int(dispaxislen / 2)
    width = 0
    binlocations = []

    sys.stdout.write(
        " >>> Determining spectrum localisation bins on dispersion axis.\n"
    )
    sys.stdout.write("     User-defined S/N threshold = " + str(minSNR) + "\n")

    # Start at the centre of the dispersion axis and start binning towards the short wavelength end of the spectrum.
    while x - width > 0:
        snrestimate = 0.0

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
            shortrows = len(
                list(
                    filter(
                        lambda x: x <= mincols,
                        np.sum(
                            fdict["qual"][:, int(x - width) : int(x)], axis=1
                        ),
                    )
                )
            )
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.
            datacol = fdict["data"][:, int(x - width) : int(x)]
            errscol = fdict["errs"][:, int(x - width) : int(x)]
            bindatacol = np.nansum(datacol, axis=1)
            binerrscol = np.sqrt(np.nansum(np.power(errscol, 2), axis=1))
            signal = np.nansum(bindatacol[int(slow) : int(shigh)])
            rssnoise = np.sqrt(
                np.nansum(np.power(binerrscol[int(slow) : int(shigh)], 2))
            )

            # Estimate the S/N
            snrestimate = signal / rssnoise

        if replace_crbp:
            # Replace bad pixels and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211

            # Here a median spatial profile for the current bin is determined by bootstrapping the good pixels
            # in each spatial pixel row with repeats to estimate a distribution of median values for each value in the
            # median spatial distribution for the bin. The mean of each of these distributions is then taken to be the
            # value of that median spatial pixel. The standard error of each of these distributions becomes the error
            # of the flux.

            meddatacol = np.nanmedian(
                fdict["data"][:, int(x - width) : int(x)], axis=1
            )
            medtile = np.tile(meddatacol, (width, 1)).T
            madtile = np.abs(
                fdict["data"][:, int(x - width) : int(x)] - medtile
            )
            errmeddatacol = np.nanmedian(madtile, axis=1)

            nmeddatacol = meddatacol / np.sum(meddatacol)

            for i in range(int(width)):
                # If there are both good and bad pixels in the current column (i.e. not all pixels are bad), repair the bad ones.
                # Leave columns of all bad pixels as they are.
                if (
                    0.0 in fdict["qual"][:, int(x - i)]
                    and all(x == 0 for x in fdict["qual"][:, int(x - i)])
                    == False
                ):
                    # print(fdict[fdict['qual'][:, int(x-i)])
                    # print(np.isfinite(fdict['data'][:, int(x - i)]))
                    cr = np.where(fdict["qual"][:, int(x - i)] == 0)
                    nocr = np.where(fdict["qual"][:, int(x - i)] == 1)
                    proportion_nocr = np.sum(nmeddatacol[nocr])
                    total_nocr = np.sum(meddatacol[nocr])

                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    fdict["data"][:, int(x - i)][cr] = total_nocr * (
                        nmeddatacol[cr] / proportion_nocr
                    )
                    fdict["errs"][:, int(x - i)][cr] = np.median(
                        fdict["errs"][:, int(x - i)][nocr]
                    )

        binlocations.append([int(x - width), int(x), snrestimate])

        x -= width
        width = 0

    x = int(dispaxislen / 2)

    # Repeat the same process as above, starting at the centre of the dispersion axis, but moving outward toward the
    # longer wavelength end of the 2D spectrum.
    while x + width < dispaxislen:
        snrestimate = 0.0

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
            shortrows = len(
                list(
                    filter(
                        lambda x: x <= mincols,
                        np.sum(
                            fdict["qual"][:, int(x) : int(x + width)], axis=1
                        ),
                    )
                )
            )
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.
            bindatacol = np.nansum(
                fdict["data"][:, int(x) : int(x + width)], axis=1
            )
            binerrscol = np.sqrt(
                np.nansum(
                    np.power(fdict["errs"][:, int(x) : int(x + width)], 2),
                    axis=1,
                )
            )
            signal = np.nansum(bindatacol[int(slow) : int(shigh)])
            rssnoise = np.sqrt(
                np.nansum(np.power(binerrscol[int(slow) : int(shigh)], 2))
            )

            # Estimate the S/N
            snrestimate = signal / rssnoise

        if replace_crbp:
            # Replace bad pixels, and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211

            meddatacol = np.nanmedian(
                fdict["data"][:, int(x - width) : int(x)], axis=1
            )
            medtile = np.tile(meddatacol, (width, 1)).T
            madtile = np.abs(
                fdict["data"][:, int(x - width) : int(x)] - medtile
            )
            errmeddatacol = np.nanmedian(madtile, axis=1)

            nmeddatacol = meddatacol / np.sum(meddatacol)

            for i in range(int(width)):
                if 0.0 in fdict["qual"][:, int(x + i)]:
                    cr = np.where(fdict["qual"][:, int(x + i)] == 0)
                    proportion_nocr = np.sum(nmeddatacol[nocr])
                    total_nocr = np.sum(meddatacol[nocr])

                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    fdict["data"][:, int(x - i)][cr] = total_nocr * (
                        nmeddatacol[cr] / proportion_nocr
                    )
                    fdict["errs"][:, int(x - i)][cr] = np.median(
                        fdict["errs"][:, int(x - i)][nocr]
                    )

        binlocations.append([int(x), int(x + width), snrestimate])

        x += width
        width = 0

    # Sort out the bins into the correct order on the dispersion axis.
    binlocations.sort(key=lambda x: x[0])

    return binlocations, fdict


def get_bins_output(
    binparams, params, lowext, highext, data2D, headparams, axdict
):
    """Print and plot the output of get_bins.

    Description: Print and plot the output of get_bins.

    Args:
        binparams (list): list containing the locations and S/N of each bin
        params (dict): dictionary of parameters read in from the motesparams.txt confifuration file.
        lowext (int): lower limit of the spatial region measured for S/N in get_bins()
        highext (int): upper limit of the spatial region measured for S/N in get_bins()
        data2D (numpy.ndarray): 2D spectroscopic data
        headparams (dict): dictionary of parameters full from the datafile header.
        axdict (dict): dictionary containing axes and axis metadata for the current extraction

    Returns:
        None
    """

    sys.stdout.write(
        " >>> "
        + str(len(binparams))
        + " spectrum localisation bins determined on dispersion axis.\n"
    )

    # DIAGNOSTICS - Plot boundary locations of localisation bins on the dispersion axis.
    if params["-DIAG_PLOT_BIN_LOC"]:
        drawlines = []
        for b in binparams:
            binlineloc = np.where(
                np.logical_and(
                    axdict["saxis"] > lowext - ((highext - lowext) * 0.2),
                    axdict["saxis"] < highext + ((highext - lowext) * 0.2),
                )
            )
            drawlines.append(
                np.ones(len(axdict["saxis"][binlineloc])) * b[0]
                + axdict["wavstart"]
            )
            drawlines.append(axdict["saxis"][binlineloc] + axdict["imgstart"])
            drawlines.append(
                np.ones(len(axdict["saxis"][binlineloc])) * b[1]
                + axdict["wavstart"]
            )
            drawlines.append(axdict["saxis"][binlineloc] + axdict["imgstart"])

        show_img(
            data2D,
            axdict,
            headparams,
            drawlines,
            "2D Spectrum with Boundaries of Localisation Bins",
        )

    return None


def interpolate_extraction_lims(extractionlims, dispaxislen):
    """Interpolate extraction limits over unbinned 2D spectrum to create the final extraction limits for the unbinned spectrum.

    Description:
    Interpolate extraction limits over unbinned 2D spectrum to create the final extraction limits for the unbinned spectrum.
    Takes an input of extraction limits from the fitting of the binned data and interpolates the limits over the unbinned data.
    Limits are also linearly extrapolated towards the ends of the spectral range.

    Args:
        extractionlims (numpy.ndarray): Extraction limits for each bin.
        dispaxislen (float): Length of the dispersion axis of the unbinned 2D spectrum.

    Returns:
        finalextlims (list): List containing the extraction limits for the unbinned 2D spectrum.
    """

    # If the 2D spectrum was so faint that only 1 dispersion bin could be determined, set the extraction limits
    # across the unbinned spectrum to be a simple linear aperture that has the same extraction limits as was
    # determined for that one bin.
    if len(extractionlims[0]) == 1:
        finalextlims = [
            np.repeat(extractionlims[1][0], dispaxislen),
            np.repeat(extractionlims[2][0], dispaxislen),
        ]

    # Otherwise, interpolate the extraction limits form the bins across the unbinned wavelength axis. Also
    # extrapolate the extraction limits in a linear fashion at the ends of the wavelength axis of the 2D spectrum
    # so that the full wavelength axis is covered.
    else:
        interpextract_1 = interp.interp1d(
            extractionlims[0], extractionlims[1], kind="linear"
        )
        interpextract_2 = interp.interp1d(
            extractionlims[0], extractionlims[2], kind="linear"
        )

        intermextlims = [
            interpextract_1(
                np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))
            ),
            interpextract_2(
                np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))
            ),
        ]

        (
            shortextraplim1,
            shortextraplim2,
            longextraplim1,
            longextraplim2,
        ) = extrap_extraction_lims(
            intermextlims,
            dispaxislen,
            extractionlims[0][0],
            extractionlims[0][-1],
        )

        extlim1 = np.insert(intermextlims[0], 0, shortextraplim1)
        extlim1 = np.append(extlim1, longextraplim1)

        extlim2 = np.insert(intermextlims[1], 0, shortextraplim2)
        extlim2 = np.append(extlim2, longextraplim2)

        nextextlims = [extlim1, extlim2]
        nextxaxis = np.insert(
            np.arange(extractionlims[0][0], extractionlims[0][-1]), 0, 0
        )
        nextxaxis = np.append(nextxaxis, dispaxislen)

        interpnextlims1 = interp.interp1d(
            nextxaxis, nextextlims[0], kind="linear", fill_value="extrapolate"
        )
        interpnextlims2 = interp.interp1d(
            nextxaxis, nextextlims[1], kind="linear", fill_value="extrapolate"
        )

        finalextlims = [
            interpnextlims1(np.array(range(dispaxislen))),
            interpnextlims2(np.array(range(dispaxislen))),
        ]

    return finalextlims


# Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a
# Levenberg-Marquardt least squares method. Returns the best fit parameters of the Moffat function
def linear_least_squares(r, col):
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(linear_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1]]


def linear_resid(x, datarange, data):
    """Calculate residuals of fitted linear profile and the data for the Levenberg-Marquardt least squares method.
    Args:
        x (list): List containing the best fit parameters of the linear profile.
        datarange (_type_):
        data (_type_):

    Returns:
        residual (_type_): The residuals of the fitted linear profile and the data.
    """
    residual = (x[0] * datarange) + x[1] - data
    return residual


def make_wav_axis(start, increment, length):
    """Returns a wavelength axis array using a start wavelength, the wavelength increment and the number of values along the axis required.

    Args:
        start (float): The start wavelength of the axis, in nm.
        increment (float): The wavelength increment between each value along the axis, in nm.
        length (int): Number of values along the axis required.

    Returns:
        wavaxis (numpy.ndarray): A wavelength axis array.
    """
    wavaxis = np.arange(
        start=start, step=increment, stop=(length * increment) + start
    )
    return wavaxis


def moffat(amp, c, alpha, beta, bglevel, bggrad, datarange):
    """Creates moffat profile added a linear sloped background based on input parameters

    Args:
        amp (float64): amplitude of the Moffat profile
        c (float64): location of the center/peak of the Moffat profile on the spatial axis.
        alpha (float64): he main parameter that defines the width of the Moffat profile.
        beta (float64): plays a role in defining the width of the Moffat profile in the outer wings far from the profile's peak.
        bglevel (float64): height of the linearly varying background
        bggrad (float64): gradient of the linearly varying background
        datarange (numpy.ndarray): x axis of the Moffat profile.

    Returns:
        moffat_background (numpy.ndarray): a moffat profile defined on the datarange axis using the parameters input to the function, with added background flux.
    """

    moffat = amp * (
        (1 + (((datarange - c) * (datarange - c)) / (alpha * alpha))) ** -beta
    )
    moffat_background = moffat + bglevel + (datarange * bggrad)

    return moffat_background


def moffat_least_squares(r, col, seeing, pixres):
    """Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a least squares method. Returns the best fit parameters of the Moffat function

    Args:
        r (_type_): spatial axis of the data being fit
        col (_type_): data being fit
        seeing (_type_): estimated FWHM of the spatial profile
        pixres (_type_): spatial resolution of each pixel in arcsec/pixel

    Returns:
        param_list (list): list of best fit output parameters returned by the least squares routine.
    """

    # Set up initial conditions for the least squares fit.
    # x0 = [amplitude, centre, alpha, beta, background gradient, background level]
    # Initial beta estimate comes from optimal value from atmospheric turbulence theory as described in
    # Trujillo, I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [
        np.nanmedian(np.sort(col)[-3:]),
        np.argmax(col),
        seeing / pixres,
        4.765,
        0.0,
        np.median(np.concatenate((col[:5], col[-5:]))),
    ]

    # Run the least squares fit.
    res_lsq = least_squares(
        moffat_resid,
        x0,
        bounds=(
            [0.0, np.argmax(col) - 1.0, 0.0, 0.0, 0.0, -np.inf],
            [
                np.inf,
                np.argmax(col) + 1,
                (5 * seeing / pixres),
                5.0,
                np.inf,
                np.inf,
            ],
        ),
        args=(r, col),
        method="trf",
        ftol=1e-12,
    )
    param_list = [
        res_lsq.x[0],
        res_lsq.x[1],
        res_lsq.x[2],
        res_lsq.x[3],
        res_lsq.x[4],
        res_lsq.x[5],
    ]
    return param_list


def moffat_resid(x, datarange, data):
    """Calculates residuals of fitted moffat profile and the data for the least squares fitting.

    Description:
        A = x[0]
        c = x[1]
        alpha = x[2]
        beta = x[3]
        B = x[4]
        m = x[5]

    Args:
        x (numpy.ndarray): an array of parameters defining the shape of the model moffat profile
        datarange (numpy.ndarray): spatial axis of the data
        data (numpy.ndarray): the data

    Returns:
        residual (numpy.ndarray): the residual array between the model moffat profile and the data
    """

    moff = x[0] * (
        (1 + ((datarange - x[1]) * (datarange - x[1])) / (x[2] * x[2])) ** -x[3]
    )
    residual = moff + x[4] + (datarange * x[5]) - data
    return residual


def optimal_extraction(data2D, errs2D, extractionlimits, binparameters, axdict):
    """Performs optimal extraction on the data2D and errs2D arrays using the extraction limits and bin parameters.

    Description:
        Perform optimal extraction using a modified version of Horne (1986) where S=0, G=0 and errors are not 'revised'
        since we already have the 'variance frame' (and we are not working with counts). Ideally, this extraction reduces
        errors with the weighting and conserves the flux (when compared to the standard extraction). This subroutine uses an
        analytic Moffat profile (post sky-subtraction), instead of polynominals as in Horne (1986). Furthermore, profiling
        takes place bin by bin, accounting for spatial profile variations across the dispersion axis; extraction takes place
        column by column within the bin limits and within the extraction limits previously calculated.

    Args:
        data2D (numpy.ndarray): Input data frame
        errs2D (numpy.ndarray): Input error frame
        extractionlimits (numpy.ndarray): contains limits at each dispersion pixel
        binparameters (list): Contains the bin limits across the dispersion axis, to enable slicing the data across dispersion axis.
        axdict (dict): dictionary containing the spatial axis array and other relevant information about the size and shape of the data frame

    Returns:
        optidata1D (numpy.ndarray): 1D array of the optimally extracted spectrum
        optierrs1D (numpy.ndarray): 1D array of the uncertainties of the optimally extracted spectrum
        aperdata1D (numpy.ndarray): 1D array of the aperture extracted spectrum
        apererrs1D (numpy.ndarray): 1D array of the uncertainties of the aperture extracted spectrum
    """

    print("Performing optimal extraction...")

    # Filter any NaNs and Inf for data/errs AND ensure the errors are positive for this extraction.
    data2D, errs2D = filter_data(data2D, np.abs(errs2D))

    # Set up output arrays for the optimally and aperture extracted spectra and their respective uncertainties
    optidata1D = np.zeros(np.shape(data2D)[0])
    optierrs1D = np.zeros(np.shape(data2D)[0])
    aperdata1D = np.zeros(np.shape(data2D)[0])
    apererrs1D = np.zeros(np.shape(data2D)[0])

    # Set up bin identification parameters for the following loop.
    bin_number = -1
    b = np.zeros(8)

    # Loop through each dispersion element of the spectrum.
    for i, col in enumerate(data2D):
        # Identify the location of the current element in the original 2D spectrum
        dpix = i + axdict["wavstart"]

        # If the current element belongs in the next bin as defined by getbins, use the new bin's parameters and increment the bin number.
        if (
            bin_number < len(binparameters) - 1
            and dpix == binparameters[bin_number + 1][-2]
        ):
            bin_number += 1
            b = binparameters[bin_number]

        # Get the extraction limits for the current dispersion element and define the spatial axis. Where the extraction limits include partial pixels on the edge of the aperture, those pixels are included in their entirety.
        loextlim = extractionlimits[0][i]
        hiextlim = extractionlimits[1][i]
        ax = axdict["saxis"][
            int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))
        ]

        # Use the extraction limits to define the data column for this wavelength element.
        col = col[int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))]

        # Use the extraction limits to define the errs column for this wavelength element.
        # Where errs have a value of 0, set them to the median err value for the entire column.
        err = errs2D[i]
        err[np.where(err == 0)] = np.median(err)
        err = err[int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))]

        # Use the err column to get the variance column
        var = err * err

        # Perform the standard aperture extraction on the data in the column, and add this value to the aperture 1D output array.
        # Get the root sum square of the err column and add this value to the uncertainty array for the 1D aperture extraction.
        f = np.sum(col)
        aperdata1D[i] += f
        apererrs1D[i] += np.sqrt(np.sum(var))

        # Step 5 of the Horne 1986 algorithm - Defining the spatial profile.
        # Use the average value of the extraction limits in the current column to estimate the location of the centre of the peak of the spatial profile.
        profcent = (hiextlim + loextlim) * 0.5

        # Use the Moffat profile parameters for the current bin to make a moffat profile that approximates the shape of the spectrum's spatial profile in the current column.
        # The estimated centre of the extraction limits is used to shift the profile to the correct location along the spatial axis.
        # The background is assumed to have been subtracted by now, so the background level and gradient are set to 0.
        profile = moffat(b[0], profcent, b[2], b[3], 0.0, 0.0, ax)
        # Enforce positivity of the data and normalise the spatial profile.
        nprof = profile / np.sum(profile)
        col[np.where(col < 0.0)] = 0.0

        # Step 7 of Horne 1986 - Cosmic ray masking with a conservative 5 sigmas.
        value = np.abs((col - (f * nprof)) ** 2)
        sigma = 5
        condition = var * sigma**2
        # With invert, it gives False (ie. 0) if the above condition is satisfied which
        # means to mask the rays
        cosmic_mask = np.invert([value > condition])[0]

        # Add an exception to handle scenarios where the mask is all False If that happens, ignore the
        # cosmic ray correction.
        if not np.any(cosmic_mask):
            cosmic_mask = np.full(shape=cosmic_mask.shape, fill_value=True)

        # Extract the optimal spectrum (Step 8 of the algorithm)
        # Page 4 of Horne 1986: These formulae are equivalent to determining the
        # OPT spectrum by scaling a known spatial profile P, to fit the sky
        # subtracted data, whose variances are V.

        # If column is in a GMOS chip gap set the optimal flux and uncertainty to 0.
        if all(x == 0.0 for x in col) or all(x == 1.0 for x in col):
            fopt = 0.0
            vopt = 0.0
        else:
            fopt = np.sum(((nprof * col) / var)[cosmic_mask]) / np.sum(
                ((nprof * nprof) / var)[cosmic_mask]
            )
            vopt = np.sum(nprof[cosmic_mask]) / np.sum(
                ((nprof * nprof) / var)[cosmic_mask]
            )
        optidata1D[i] += fopt
        optierrs1D[i] += np.sqrt(vopt)

    return optidata1D, optierrs1D, aperdata1D, apererrs1D


def plot_fitted_spatial_profile(
    spataxis, bindata, hiresspataxis, binmoffparams, imgstart, headparams
):
    """Plot the spatial profile of a collapsed spectrum or a collapsed bin therein, and plot the Moffat function fitted to the data on top.

    Args:
        spataxis (_type_): the spatial, or x, axis of the profile.
        bindata (_type_): the binned data that has been fitted with a Moffat profile.
        hiresspataxis (_type_): A supersampled spatial axis used only for plotting purposes.
        binmoffparams (_type_): Parameters defining the Moffat profiel that was fitted to bindata.
        imgstart (_type_): ower limit of the spatial axis after the original 2D spectrum was cut down to the region defined in reg.txt
        headparams (_type_): dictionary of parameters pulled from the header of the current datafile

    Returns:
        None
    """
    plt.figure(figsize=(7, 4.5))
    plt.plot(
        hiresspataxis + imgstart,
        moffat(
            binmoffparams[0],
            binmoffparams[1],
            binmoffparams[2],
            binmoffparams[3],
            binmoffparams[4],
            binmoffparams[5],
            hiresspataxis,
        ),
        color="r",
        linewidth=3,
        label="Fitted Moffat Profile",
    )

    plt.plot(spataxis + imgstart, bindata, color="k", label="Spatial Profile")
    plt.grid(linestyle="dashed", color="gray")
    plt.legend()
    plt.title("Spectrum Spatial Profile and Fitted Moffat Profile")
    plt.xlabel("Spatial Axis, Pixels")
    plt.ylabel("Median Flux, " + headparams["fluxunit"])
    plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #


def poly2_least_squares(r, col):
    """Fits a second order polynomial to a data column using a Levenberg-Marquardt least squares method.

    Description:
    Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a Levenberg-Marquardt least squares method.

    Args:
        r (_type_): the spatial axis of the data column.
        col (_type_): the data column.

    Returns:
       params (list): the parameters of the fitted polynomial.

    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly2_resid, x0, args=(r, col), method="trf")

    params = [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2]]
    return params


def poly2_resid(x, datarange, data):
    """Calculates residuals of fitted linear profile and the data for the Levenberg Marquardt least squares method.

    Description:
    m = x[0]
    c = x[1]

    Args:
        x (_type_): the parameters of the fitted polynomial.
        datarange (_type_): the spatial axis of the data column.
        data (_type_): the data column.

    Returns:
        resid (_type_): the residuals of the fitted polynomial and the data.
    """
    resid = (x[0] * datarange * datarange) + (x[1] * datarange) + x[2] - data

    return resid


def poly3_least_squares(r, col):
    """Fits a third order polynomial to a data column using a Levenberg-Marquardt least squares method.

    Description:
    Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a Levenberg-Marquardt least squares method.

    Args:
        r (_type_): the spatial axis of the data column.
        col (_type_): the data column.

    Returns:
         params (list): the parameters of the fitted polynomial.

    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly3_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2], res_lsq.x[3]]


def poly3_resid(x, datarange, data):
    """Calculates residuals of fitted linear profile and the data for the Levenberg Marquardt least squares method.

    Description:
    m = x[0]
    c = x[1]

    Args:
        x (_type_): the parameters of the fitted polynomial.
        datarange (_type_): the spatial axis of the data column.
        data (_type_): the data column.

    Returns:
        resid (_type_): the residuals of the fitted polynomial and the data.
    """
    resid = (
        (x[0] * datarange * datarange * datarange)
        + (x[1] * datarange * datarange)
        + (x[2] * datarange)
        + x[3]
        - data
    )
    return resid


def printmoffparams(moffparams, imgstart, datascale):
    """Prints, to the terminal, a set parameters describing a fitted Moffat function.

    Description:
    Takes a list of Moffat profile parameters, the lower limit of the spatial axis after the original 2D spectrum was cut down to the region defined in reg.txt and the multiplier used to scale the spatial profile so it could be fit with a Moffat profile using scipy least_squares and prints them to the terminal.

    Args:
        moffparams (list): the parameters of the fitted Moffat profile.
        imgstart (int): the lower limit of the spatial axis after the original 2D spectrum was cut down to the region defined in reg.txt.
        datascale (float): the multiplier used to scale the spatial profile so it could be fit with a Moffat profile using scipy least_squares.

    Returns:
        None
    """
    sys.stdout.write(" >>> Fitted Moffat function parameters:\n")
    sys.stdout.write("         A = " + str(moffparams[0]) + "\n")
    sys.stdout.write("         c = " + str(moffparams[1] + imgstart) + "\n")
    sys.stdout.write("     alpha = " + str(moffparams[2]) + "\n")
    sys.stdout.write("      beta = " + str(moffparams[3]) + "\n")
    sys.stdout.write("         B = " + str(moffparams[4]) + "\n")
    sys.stdout.write("         m = " + str(moffparams[5]) + "\n\n")
    sys.stdout.write(
        " >>> Profile scaling factor used for fitting: " + str(datascale) + "\n"
    )
    sys.stdout.write(
        " >>> Plot of median spatial profile presents the orginal unscaled profile.\n"
    )

    return None


# Takes an input image and line data to be drawn on that image and creates a figure to be shown on screen.
def show_img(data2D, axdict, headparams, drawlines, title):
    """Creates a figure to be shown on screen.

    Description:
    Takes an input image and line data to be drawn on that image and creates a figure to be shown on screen.

    Args:
        data2D (numpy.ndarray): the input image.
        axdict (dict): a dictionary containing the spatial and spectral axes of the input image.
        headparams (dict): a dictionary containing the header parameters of the input image.
        drawlines (list): a list of line data to be drawn on the input image.
        title (str): the title of the figure.

    Returns:
        None
    """

    # Catch and suppress the following UserWarning:
    # UserWarning: Warning: 'partition' will ignore the 'mask' of the MaskedArray.
    #     a.partition(kth, axis=axis, kind=kind, order=order)
    # Numpy masked array used only when displaying the spectrograms.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        warnings.warn("partition", UserWarning)

        power = int(np.floor(np.log10(np.abs(np.nanmean(data2D)))))
        data2D = copy.deepcopy(data2D) / 10**power

        figwidth = 10.0
        fig = plt.figure(figsize=(figwidth, figwidth / 1.9))
        gs = gridspec.GridSpec(18, 33)
        ax = plt.subplot(gs[:, :32])
        colax = plt.subplot(gs[1:, 32])
        masked_data2D = np.ma.masked_where(data2D == 0, data2D)
        cmap = matplotlib.cm.inferno
        cmap.set_bad(color="red")
        s = ax.imshow(
            masked_data2D,
            aspect="auto",
            vmin=0,
            vmax=np.nanmedian(masked_data2D) + (0.5 * np.nanstd(masked_data2D)),
            origin="lower",
            cmap=cmap,
            extent=[
                axdict["wavstart"],
                axdict["wavstart"] + len(axdict["waxis"]),
                axdict["saxis"][0] + axdict["imgstart"],
                axdict["saxis"][-1] + axdict["imgstart"],
            ],
        )

        for i in range(int(len(drawlines) / 2)):
            ax.plot(drawlines[i * 2], drawlines[(i * 2) + 1], color="white")
        cbar = fig.colorbar(s, cax=colax)
        cbar.ax.yaxis.set_offset_position("left")
        cbar.ax.set_ylabel(
            "Pixel Flux, x10^" + str(power) + " " + headparams["fluxunit"]
        )
        ax2 = ax.twiny()
        ax2.plot(axdict["waxis"], data2D[0, :], alpha=0)
        ax2.set_xlim(axdict["waxis"][0], axdict["waxis"][-1])
        ax2.set_xlabel("Wavelength, " + headparams["wavunit"])
        ax.set_ylim(
            axdict["saxis"][0] + axdict["imgstart"],
            axdict["saxis"][-1] + axdict["imgstart"],
        )
        ax.set_ylabel("Spatial Axis, Pixels")
        ax.set_xlim(
            axdict["wavstart"], axdict["wavstart"] + len(axdict["waxis"])
        )

        ax.set_xlabel("Dispersion Axis, Pixels")

        plt.title(title, y=1.095)

        # Add interactive scaling bar to figures if the 2D spectrum isn't flux calibrated.
        # Again for some reason teeny tiny numbers cause things to break.
        fig.subplots_adjust(bottom=0.2)
        axvmin = plt.axes([0.1, 0.05, 0.8, 0.03])
        axvmax = plt.axes([0.1, 0.01, 0.8, 0.03])
        smin = Slider(
            axvmin,
            "LowCut",
            0,
            np.nanmax(masked_data2D) - 1,
            valinit=0.0,
            valstep=0.001 * np.nanmax(masked_data2D),
        )
        smax = Slider(
            axvmax,
            "HighCut",
            1.0,
            np.nanmax(masked_data2D),
            valinit=np.nanmedian(masked_data2D)
            + (3.0 * np.nanstd(masked_data2D)),
            valstep=0.001 * np.nanmax(masked_data2D),
        )

        def update(val):
            vmax = smax.val
            smin.valmax = vmax - 1
            vmin = smin.val
            smax.valmin = vmin + 1

            s.set_clim(vmin, vmax)
            cbar.ax.set_ylabel(
                "Pixel Flux, x10^" + str(power) + " " + headparams["fluxunit"]
            )
            fig.canvas.draw_idle()

        smin.on_changed(update)
        smax.on_changed(update)

        plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
def subtract_sky(bglowext, bghighext, fdict, axdict, pars, hpars):
    """Subtracts the sky background from the 2D image.

    Description:
        Subtracts the sky background from the 2D image by defining bg regions using limits input to the function and then
        fitting a profile to the background column by column while masking cosmic rays. The fitted linear for
        each column is subtracted from the full column to produce a background subtracted 2D image.

    Args:
        bglowext (numpy.ndarray): Lower limits of the background regions.
        bghighext (numpy.ndarray): Upper limits of the background regions.
        fdict (dict): Dictionary containing the 2D image.
        axdict (dict): Dictionary containing the axis information.
        pars (dict): Dictionary containing MOTES parameters.
        hpars (dict): Dictionary containing header information.

    Returns:
        fdict (dict): Dictionary containing the background subtracted 2D image.
    """

    print("Subtracting sky background...")

    fdict["data"] = fdict["data"].T
    fdict["errs"] = fdict["errs"].T

    medsky = []
    colnum = len(bglowext)

    # Makes sure the limits are within the image.
    for ii in range(colnum):
        if bglowext[ii] < 0:
            bglowext[ii] = 0
        if bghighext[ii] > axdict["saxis"][-1]:
            bghighext[ii] = axdict["saxis"][-1]

        datacol = fdict["data"][ii]
        colrange = np.array(range(len(datacol)))
        skypix = datacol[
            np.where(
                np.logical_or(colrange < bglowext[ii], colrange > bghighext[ii])
            )
        ]

        # Kill MOTES if there isn't enough background sky to perform sky subtraction. Should probably be made more nuanced later.
        if len(skypix) == 0:
            sys.stdout.write(" >>> No pixels contained inside sky region.\n")
            sys.stdout.write(
                "     -BG_FWHM_MULTIPLIER in motesparams.txt is probably too large.\n"
            )
            sys.stdout.write(
                "     Please reduce -BG_FWHM_MULTIPLIER and try again.\n"
            )
            sys.stdout.write(
                "     -BG_FWHM_MULTIPLIER < "
                + str(
                    round(
                        np.min(
                            [
                                np.shape(fdict["data"])[0],
                                np.shape(fdict["data"])[1],
                            ]
                        )
                        / (2 * hpars["seeing"]),
                        1,
                    )
                )
                + " recommended in this case.\n"
            )
            sys.stdout.write(
                "     Enlarging the 2D spectrum region in reg.txt is also a viable solution.\n"
            )
            sys.stdout.write("     Terminating MOTES.\n")
            sys.exit()

        skyrange = colrange[
            np.where(
                np.logical_or(colrange < bglowext[ii], colrange > bghighext[ii])
            )
        ]
        if len(set(skypix)) == 1:
            continue

        else:
            medpix = np.nanmedian(skypix)
            stdpix = np.nanstd(skypix)
            loc = np.where(
                np.logical_and(
                    skypix > medpix - (10 * stdpix),
                    skypix < medpix + (10 * stdpix),
                )
            )
            skypix = skypix[loc]
            skyrange = skyrange[loc]

        if pars["-SKYSUB_MODE"] == "MEDIAN":
            bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
            skysamp = np.nanmedian(bootsky, axis=0)
            skylevel = np.nanmean(skysamp)
            medsky.append(skylevel)
            skyerr = np.std(skysamp) / (99**0.5)

        if pars["-SKYSUB_MODE"] == "LINEAR":
            bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
            grads = []
            intercepts = []
            for jj in bootsky.T:
                linpars = linear_least_squares(skyrange, jj)
                grads.append(linpars[0])
                intercepts.append(linpars[1])
            intercepts = np.array(intercepts)
            grads = np.array(grads)
            skygrad = np.mean(grads)
            skygraderr = np.std(grads) / (99**0.5)
            skyint = np.mean(intercepts)
            skyinterr = np.std(intercepts) / (99**0.5)
            skylevel = (skygrad * colrange) + skyint
            medsky.append(skylevel)
            skyerr = (
                (skygraderr * colrange * skygraderr * colrange)
                + (skyinterr * skyinterr)
            ) ** 0.5

        if pars["-SKYSUB_MODE"] == "POLY2":
            bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
            grads = []
            intercepts = []
            quads = []
            for jj in bootsky.T:
                linpars = poly2_least_squares(skyrange, jj)
                quads.append(linpars[0])
                grads.append(linpars[1])
                intercepts.append(linpars[2])
            quads = np.array(quads)
            intercepts = np.array(intercepts)
            grads = np.array(grads)
            skyquad = np.mean(quads)
            skyquaderr = np.std(quads) / (99**0.5)
            skygrad = np.mean(grads)
            skygraderr = np.std(grads) / (99**0.5)
            skyint = np.mean(intercepts)
            skyinterr = np.std(intercepts) / (99**0.5)
            skylevel = (
                (skyquad * colrange * colrange) + (skygrad * colrange) + skyint
            )
            medsky.append(skylevel)
            skyerr = (
                (
                    skyquaderr
                    * skyquaderr
                    * colrange
                    * colrange
                    * colrange
                    * colrange
                )
                + (skygraderr * colrange * skygraderr * colrange)
                + (skyinterr * skyinterr)
            ) ** 0.5

        if pars["-SKYSUB_MODE"] == "POLY3":
            bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
            grads = []
            intercepts = []
            quads = []
            trips = []
            for jj in bootsky.T:
                linpars = poly3_least_squares(skyrange, jj)
                trips.append(linpars[0])
                quads.append(linpars[1])
                grads.append(linpars[2])
                intercepts.append(linpars[3])
            trips = np.array(trips)
            quads = np.array(quads)
            intercepts = np.array(intercepts)
            grads = np.array(grads)
            skytrip = np.mean(trips)
            skytriperr = np.std(trips) / (99**0.5)
            skyquad = np.mean(quads)
            skyquaderr = np.std(quads) / (99**0.5)
            skygrad = np.mean(grads)
            skygraderr = np.std(grads) / (99**0.5)
            skyint = np.mean(intercepts)
            skyinterr = np.std(intercepts) / (99**0.5)
            skylevel = (
                (skytrip * colrange * colrange * colrange)
                + (skyquad * colrange * colrange)
                + (skygrad * colrange)
                + skyint
            )
            medsky.append(skylevel)
            skyerr = (
                (
                    skytriperr
                    * skytriperr
                    * colrange
                    * colrange
                    * colrange
                    * colrange
                    * colrange
                    * colrange
                )
                + (
                    skyquaderr
                    * skyquaderr
                    * colrange
                    * colrange
                    * colrange
                    * colrange
                )
                + (skygraderr * colrange * skygraderr * colrange)
                + (skyinterr * skyinterr)
            ) ** 0.5

        fdict["data"][ii] -= skylevel
        fdict["errs"][ii] = ((fdict["errs"][ii] ** 2) + (skyerr**2)) ** 0.5

        sys.stdout.write(
            "     " + str(ii + 1) + "/" + str(colnum) + " columns completed.\r"
        )

    medsky = np.array(medsky)
    if pars["-SKYSUB_MODE"] == "MEDIAN":
        fdict["skymod"] = np.tile(medsky, (np.shape(fdict["data"])[1], 1))
    else:
        fdict["skymod"] = medsky

    fdict["data"] = fdict["data"].T
    fdict["errs"] = fdict["errs"].T

    chipgaps = np.where(np.nanmedian(fdict["data"], axis=0) == 0)
    fdict["data"][:, chipgaps[0]] += 1

    return fdict
