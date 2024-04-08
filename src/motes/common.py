"""
common.py - Common functions for the MOTES pipeline.
"""

import copy
import sys
import warnings

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from matplotlib.widgets import Slider
from scipy.optimize import least_squares


def extraction_limits(moffat_parameters, width_multiplier=3.0):
    """
    Calculate the extraction limits from a Moffat profile based on the distance from the central
    peak as a multiple of FWHM.

    Args:
        moffat_parameters (list)                  : list containing the parameters of the Moffat profile
                                             to be created and measured. The parameters are:
                                             [amplitude, location, scale, power].
        width_multiplier (float, optional) : defines the distance from the center of the spatial
                                             profile at which to set the extraction limits, in
                                             multiples of the FWHM. Defaults to 3.0.

    Returns:
        lower_extraction_limit (numpy.float64) : the lower bound of the region to be extracted.
        upper_extraction_limit (numpy.float64) : the upper bound of the region to be extracted.
        fwhm (numpy.float64)                   : the Full Width at Half Maximum of the Moffat
                                                 profile.
        moffat_parameters[1] (numpy.float64)          : location of the center of the Moffat profile.
    """

    # Create a Moffat profile based on the input parameters.
    fwhm = 2 * moffat_parameters[2] * (((2 ** (1 / moffat_parameters[3])) - 1) ** 0.5)

    # Respectively define the upper and lower extraction limits at a distance above and below the
    # peak of the Moffat profile that equals the FWHM of the Moffat profile times a multiplying
    # factor.
    lower_extraction_limit = moffat_parameters[1] - (width_multiplier * fwhm)
    upper_extraction_limit = moffat_parameters[1] + (width_multiplier * fwhm)

    return lower_extraction_limit, upper_extraction_limit, fwhm, moffat_parameters[1]


def extrapolate_extraction_limits(extraction_limits, dispersion_axis_length, short_end, long_end):
    """
    Linearly extrapolate the extraction limits at the ends of the 2D spectrum.

    Args:
        extraction_limits (list)    : A list containing the lower and upper extraction limits for each
                            spatial pixel.
        dispersion_axis_length (int) : The length of the dispersion axis.
        short_end (int)    : The number of pixels at the short end of the dispersion axis to be
                            excluded from the extraction.
        long_end (int)     : The number of pixels at the long end of the dispersion axis to be
                            excluded from the extraction.

    Returns:
        parameters (dict):

    """
    extrapolation_gradient_short_1 = extrapolate_gradient(extraction_limits[0], [0, 150, 300])
    extrapolation_gradient_short_2 = extrapolate_gradient(extraction_limits[1], [0, 150, 300])
    extrapolation_gradient_long_1 = extrapolate_gradient(extraction_limits[0], [-300, -150, -1])
    extrapolation_gradient_long_2 = extrapolate_gradient(extraction_limits[1], [-300, -150, -1])

    extrapolation_limit_short_1 = extraction_limits[0][0] - (extrapolation_gradient_short_1 * short_end)
    extrapolation_limit_short_2 = extraction_limits[1][0] - (extrapolation_gradient_short_2 * short_end)
    extrapolation_limit_long_1 = extraction_limits[0][-1] + (extrapolation_gradient_long_1 * (dispersion_axis_length - long_end))
    extrapolation_limit_long_2 = extraction_limits[1][-1] + (extrapolation_gradient_long_2 * (dispersion_axis_length - long_end))

    return (
        extrapolation_limit_short_1,
        extrapolation_limit_short_2,
        extrapolation_limit_long_1,
        extrapolation_limit_long_2,
    )


def extrapolate_gradient(data_array, data_limits):
    """
    Given a range of data and limits to define a region of that data, calculate the region's
    gradient.

    Args:
        data_array (numpy.ndarray) : A 1D array of data.
        data_limits (list)         : A list containing the limits of the region of data to be used
                                     to calculate the gradient.

    Returns:
       gradient (int): The gradient of the region of data.

    """
    median_of_y_points_x_to_y = np.median(data_array[data_limits[0] : data_limits[1]])
    median_of_y_points_y_to_z = np.median(data_array[data_limits[1] : data_limits[2]])
    median_of_x_points_x_to_y = np.median([data_limits[0], data_limits[1]])
    median_of_x_points_y_to_z = np.median([data_limits[1], data_limits[2]])

    gradient = (median_of_y_points_x_to_y - median_of_y_points_y_to_z) / (
        median_of_x_points_x_to_y - median_of_x_points_y_to_z
    )

    return gradient


def filter_data(data_2D, errs_2D):
    """
    This function takes in the data_2D and errs_2D and outputs frames where any NaN or Inf values are
    0.0 (this can be later filtered to a median of the column). This is used before the extraction
    procedures to ensure the S/N in the bin is numerical (required for optimal extraction).

    Args:
        data_2D (numpy.ndarray) : Original data_2D
        errs_2D (numpy.ndarray) : Original errs_2D

    Returns:
        data_2D (numpy.ndarray) : Filtered data_2D
        errs_2D (numpy.ndarray) : Filtered errs_2D
    """

    errs_2D[~np.isfinite(errs_2D)] = 0.0
    data_2D[~np.isfinite(errs_2D)] = 0.0
    data_2D[~np.isfinite(data_2D)] = 0.0
    errs_2D[~np.isfinite(data_2D)] = 0.0

    return data_2D, errs_2D


def get_bins(frame_dict, spatial_lo_limit, spatial_hi_limit, dispersion_axis_length, parameter_dict, has_sky=False):
    """
    Define the bins of data over which Moffat profiles will be fitted. Each bin is defined such
    that when summed it will have a given signal to noise (S/N). So lower S/N regions will have
    larger bins. This function works from the centre of the input 2D spectrum outward to the ends
    in order to ensure a good start for the binning process.

    Args:
        frame_dict (dict)                  : Dictionary containing all the data, error, and quality
                                        frames.
        spatial_lo_limit (int)                    : Lower spatial limit of the region where the S/N will be
                                        measured when defining the extent of a bin.
        spatial_hi_limit (int)                   : Upper spatial limit of the region where the S/N will be
                                        measured when defining the extent of a bin.
        dispersion_axis_length (int)             : Length of the dispersion axis.
        parameter_dict (dict)                 : Dictionary of parameters ready in from the motesparams.txt
                                        configuration file.
        has_sky (bool, optional)      : Used to tell get_bins whether the sky has been subtracted
                                        yet or not, and to set the minimum_SNR threshold accordingly.
                                        False by default.

    Returns:
        bin_locations (list) : A list containing the details for each bin determined by get_bins.
                              The boundaries and S/N of each bin are recorded here.
        frame_dict (dict)        : Returns frame_dict.
    """

    # Take S/N threshold (minimum_SNR) and minimum number of columns per dispersion bin (minimum_columns)
    if parameter_dict["-SUBTRACT_SKY"] and has_sky:
        minimum_SNR = parameter_dict["-SKY_SNR_BIN_LIM"]
    else:
        minimum_SNR = parameter_dict["-SNR_BIN_LIM"]

    # Minimum columns for a bin
    minimum_columns = parameter_dict["-COL_BIN_LIM"]

    # Start x at the central pixel column of the dispersion axis
    x = int(dispersion_axis_length / 2)
    width = 0
    bin_locations = []

    sys.stdout.write(
        " >>> Determining spectrum localisation bins on dispersion axis.\n"
    )
    sys.stdout.write("     User-defined S/N threshold = " + str(minimum_SNR) + "\n")

    # Start at the centre of the dispersion axis and start binning towards the short wavelength end
    # of the spectrum.
    while x - width > 0:
        snr_estimate = 0.0

        # If the S/N of the current bin has not yet reached the user-defined threshold (minimum_SNR),
        # add one more pixel column to the bin.
        while snr_estimate <= minimum_SNR:
            width += 1

            # Stop the loop once the short wavelength end of the spectrum os reached.
            if x - width < 0:
                width = int(0 + x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin,
            # continue to the next iteration and add another pixel column to the bin.
            short_rows = len(
                list(
                    filter(
                        lambda x: x <= minimum_columns,
                        np.sum(frame_dict["qual"][:, int(x - width) : int(x)], axis=1),
                    )
                )
            )
            if short_rows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the
            # sum of the flux in the bin and the noise is the root sum square of the errors on the
            # flux in the bin. Errors are taken from the spectrum's error frame.
            data_columns = frame_dict["data"][:, int(x - width) : int(x)]
            errs_columns = frame_dict["errs"][:, int(x - width) : int(x)]
            binned_data_columns = np.nansum(data_columns, axis=1)
            binned_errs_columns = np.sqrt(np.nansum(np.power(errs_columns, 2), axis=1))
            signal = np.nansum(binned_data_columns[int(spatial_lo_limit) : int(spatial_hi_limit)])
            root_sum_square_noise = np.sqrt(
                np.nansum(np.power(binned_errs_columns[int(spatial_lo_limit) : int(spatial_hi_limit)], 2))
            )

            # Estimate the S/N
            snr_estimate = signal / root_sum_square_noise

        bin_locations.append([int(x - width), int(x), snr_estimate])

        x -= width
        width = 0

    x = int(dispersion_axis_length / 2)

    # Repeat the same process as above, starting at the centre of the dispersion axis, but moving
    # outward toward the longer wavelength end of the 2D spectrum.
    while x + width < dispersion_axis_length:
        snr_estimate = 0.0

        # If the S/N of the current bin has not yet reached the user-defined threshold (minimum_SNR),
        # add one more pixel column to the bin.
        while snr_estimate <= minimum_SNR:
            width += 1

            # Stop the loop once the long wavelength end of the spectrum is reached.
            if x + width > dispersion_axis_length:
                width = int(dispersion_axis_length - x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin,
            # continue to the next iteration and add another pixel column to the bin.
            short_rows = len(
                list(
                    filter(
                        lambda x: x <= minimum_columns,
                        np.sum(frame_dict["qual"][:, int(x) : int(x + width)], axis=1),
                    )
                )
            )
            if short_rows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the
            # sum of the flux in the bin and the noise is the root sum square of the errors on the
            # flux in the bin. Errors are taken from the spectrum's error frame.
            binned_data_columns = np.nansum(frame_dict["data"][:, int(x) : int(x + width)], axis=1)
            binned_errs_columns = np.sqrt(
                np.nansum(
                    np.power(frame_dict["errs"][:, int(x) : int(x + width)], 2),
                    axis=1,
                )
            )
            signal = np.nansum(binned_data_columns[int(spatial_lo_limit) : int(spatial_hi_limit)])
            root_sum_square_noise = np.sqrt(
                np.nansum(np.power(binned_errs_columns[int(spatial_lo_limit) : int(spatial_hi_limit)], 2))
            )

            # Estimate the S/N
            snr_estimate = signal / root_sum_square_noise

        bin_locations.append([int(x), int(x + width), snr_estimate])

        x += width
        width = 0

    # Sort out the bins into the correct order on the dispersion axis.
    bin_locations.sort(key=lambda x: x[0])

    return bin_locations, frame_dict


def get_bins_output(binparams, params, lowext, highext, data2D, headparams, axdict):
    """
    Print and plot the output of get_bins.

    Args:
        binparams (list)       : A list containing the locations and S/N of each bin
        params (dict)          : A dictionary of parameters read in from the motesparams.txt
                                 configuration file.
        lowext (int)           : lower limit of the spatial region measured for S/N in get_bins()
        highext (int)          : upper limit of the spatial region measured for S/N in get_bins()
        data2D (numpy.ndarray) : 2D spectroscopic data
        headparams (dict)      : A dictionary of parameters full from the datafile header.
        axdict (dict)          : A dictionary containing axes and axis metadata for the current
                                 extraction

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
                    axdict["spatial_axis"] > lowext - ((highext - lowext) * 0.2),
                    axdict["spatial_axis"] < highext + ((highext - lowext) * 0.2),
                )
            )
            drawlines.append(
                np.ones(len(axdict["spatial_axis"][binlineloc])) * b[0] + axdict["wavelength_start"]
            )
            drawlines.append(axdict["spatial_axis"][binlineloc] + axdict["data_spatial_floor"])
            drawlines.append(
                np.ones(len(axdict["spatial_axis"][binlineloc])) * b[1] + axdict["wavelength_start"]
            )
            drawlines.append(axdict["spatial_axis"][binlineloc] + axdict["data_spatial_floor"])

        show_img(
            data2D,
            axdict,
            headparams,
            drawlines,
            "2D Spectrum with Boundaries of Localisation Bins",
        )

    return None


def interpolate_extraction_lims(extractionlims, dispaxislen):
    """
    Interpolate extraction limits over unbinned 2D spectrum to create the final extraction limits
    for the unbinned spectrum. Takes an input of extraction limits from the fitting of the binned
    data and interpolates the limits over the unbinned data. Limits are also linearly extrapolated
    towards the ends of the spectral range.

    Args:
        extractionlims (numpy.ndarray) : Extraction limits for each bin.
        dispaxislen (float)            : Length of the dispersion axis of the unbinned 2D spectrum.

    Returns:
        finalextlims (list)            : List containing the extraction limits for the unbinned 2D
                                         spectrum.
    """

    # If the 2D spectrum was so faint that only 1 dispersion bin could be determined, set the
    # extraction limits across the unbinned spectrum to be a simple linear aperture that has the
    # same extraction limits as was determined for that one bin.
    if len(extractionlims[0]) == 1:
        finalextlims = [
            np.repeat(extractionlims[1][0], dispaxislen),
            np.repeat(extractionlims[2][0], dispaxislen),
        ]

    # Otherwise, interpolate the extraction limits form the bins across the unbinned wavelength
    # axis. Also extrapolate the extraction limits in a linear fashion at the ends of the
    # wavelength axis of the 2D spectrum so that the full wavelength axis is covered.
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
        ) = extrapolate_extraction_limits(
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


# Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the
# column using a Levenberg-Marquardt least squares method. Returns the best fit parameters of the
# Moffat function
def linear_least_squares(r, col):
    """
    Fit a linear profile to a data column using a Levenberg-Marquardt least squares method.

    Args:
        r (numpy.ndarray)   : Spatial axis of the data column being fitted.
        col (numpy.ndarray) : Data column being fitted.

    Returns:
        [res_lsq.x[0], res_lsq.x[1]] (list) : list containing the best fit parameters of the linear
                                              profile.
    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(linear_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1]]


def linear_resid(x, datarange, data):
    """
    Calculate residuals of fitted linear profile and the data for the Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the linear
                                    profile.
        datarange (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray)  : The residuals of the fitted line and the data.
    """
    residual = (x[0] * datarange) + x[1] - data
    return residual


def poly2_least_squares(r, col):
    """
    Fits a second-order polynomial to a data column using a Levenberg-Marquardt least squares
    method.

    Args:
        r (numpy.ndarray)   : Spatial axis of the data column being fitted.
        col (numpy.ndarray) : Data column being fitted.

    Returns:
       params (list)        : the parameters of the fitted polynomial.
    """

    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly2_resid, x0, args=(r, col), method="trf")

    params = [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2]]
    return params


def poly2_resid(x, datarange, data):
    """
    Calculates residuals of a fitted second-order polynomial and the data for the
    Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the polynomial.
        datarange (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray) : The residuals of the fitted polynominal and the data.
    """
    residual = (x[0] * datarange * datarange) + (x[1] * datarange) + x[2] - data

    return residual


def poly3_least_squares(r, col):
    """
    Fits a third-order polynomial to a data column using a Levenberg-Marquardt least squares
    method.

    Args:
        r (numpy.ndarray)   : Spatial axis of the data column being fitted.
        col (numpy.ndarray) : Data column being fitted.

    Returns:
        params (list) : the best-fit parameters of the fitted polynomial.
    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly3_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2], res_lsq.x[3]]


def poly3_resid(x, datarange, data):
    """
    Calculates residuals of a fitted third-order polynomial and the data for the
    Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the polynomial.
        datarange (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray) : The residuals of the fitted polynominal and the data.
    """
    residual = (
        (x[0] * datarange * datarange * datarange)
        + (x[1] * datarange * datarange)
        + (x[2] * datarange)
        + x[3]
        - data
    )
    return residual


def make_wav_axis(start, increment, length):
    """
    Returns a wavelength axis array using a start wavelength, the wavelength increment and the
    number of values required along the axis.

    Args:
        start (float)     : The start wavelength of the axis, in nm.
        increment (float) : The wavelength increment between each value along the axis.
        length (int)      : Number of values along the axis required.

    Returns:
        wavaxis (numpy.ndarray) : A wavelength axis array.
    """
    wavaxis = np.arange(start=start, step=increment, stop=(length * increment) + start)
    return wavaxis


def moffat(amp, c, alpha, beta, bglevel, bggrad, datarange):
    """
    Creates moffat profile added a linear sloped background based on input parameters

    Args:
        amp (float64)             : amplitude of the Moffat profile
        c (float64)               : location of the center/peak of the Moffat profile on the
                                    spatial axis.
        alpha (float64)           : he main parameter that defines the width of the Moffat profile.
        beta (float64)            : plays a role in defining the width of the Moffat profile in the
                                    outer wings far from the profile's peak.
        bglevel (float64)         : height of the linearly varying background
        bggrad (float64)          : gradient of the linearly varying background
        datarange (numpy.ndarray) : x axis of the Moffat profile.

    Returns:
        moffat_background (numpy.ndarray) : a moffat profile defined on the datarange axis using
                                            the parameters input to the function, with added
                                            background flux.
    """

    moffat = amp * (
        (1 + (((datarange - c) * (datarange - c)) / (alpha * alpha))) ** -beta
    )
    moffat_background = moffat + bglevel + (datarange * bggrad)

    return moffat_background


def moffat_least_squares(r, col, seeing, pixres):
    """
    Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to
    the column using a least squares method. Returns the best fit parameters of the Moffat
    function.

    Args:
        r (numpy.ndarray)   : spatial axis of the data being fit
        col (numpy.ndarray) : data being fitted
        seeing (float)      : estimated FWHM of the spatial profile
        pixres (float)      : spatial resolution of each pixel in arcsec/pixel

    Returns:
        param_list (list) : list of best fit output parameters returned by the least squares
                            routine.
    """

    # Set up initial conditions for the least squares fit.
    # x0 = [amplitude, centre, alpha, beta, background gradient, background level]
    # Initial beta estimate comes from optimal value from atmospheric turbulence theory as
    # described in Trujillo, I. et al. (2001), MNRAS, 328, 977-985
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
        x (numpy.ndarray)         : an array of parameters defining the shape of the model moffat
                                    profile
        datarange (numpy.ndarray) : spatial axis of the data
        data (numpy.ndarray)      : the data

    Returns:
        residual (numpy.ndarray) : the residual array between the model moffat profile and the data
    """

    moff = x[0] * (
        (1 + ((datarange - x[1]) * (datarange - x[1])) / (x[2] * x[2])) ** -x[3]
    )
    residual = moff + x[4] + (datarange * x[5]) - data
    return residual


def optimal_extraction(data2D, errs2D, extractionlimits, binparameters, axdict):
    """
    Perform optimal extraction using a modified version of Horne (1986) where S=0, G=0 and errors
    are not 'revised' since we already have the 'variance frame'. Ideally, this extraction reduces
    errors with the weighting and conserves the flux (when compared to the standard extraction).
    This subroutine uses an analytic Moffat profile (post sky-subtraction), instead of polynominals
    as in Horne (1986). Furthermore, profiling takes place bin by bin, accounting for spatial
    profile variations across the dispersion axis; extraction takes place column by column within
    the bin limits and within the extraction limits previously calculated.

    Args:
        data2D (numpy.ndarray)           : Input data frame
        errs2D (numpy.ndarray)           : Input error frame
        extractionlimits (numpy.ndarray) : An array containing limits at each dispersion pixel
        binparameters (list)             : A list containing the bin limits across the dispersion
                                           axis, to enable slicing the data across dispersion axis.
        axdict (dict)                    : A dictionary containing the spatial axis array and other
                                           relevant information about the size and shape of the
                                           data frame

    Returns:
        optidata1D (numpy.ndarray) : 1D array of the optimally extracted spectrum
        optierrs1D (numpy.ndarray) : 1D array of the uncertainties of the optimally extracted
                                     spectrum
        aperdata1D (numpy.ndarray) : 1D array of the aperture extracted spectrum
        apererrs1D (numpy.ndarray) : 1D array of the uncertainties of the aperture extracted
                                     spectrum
    """

    sys.stdout.write("Performing optimal extraction...")

    # Filter any NaNs and Inf for data/errs AND ensure the errors are positive for this extraction.
    data2D, errs2D = filter_data(data2D, np.abs(errs2D))

    # Set up output arrays for the optimally and aperture extracted spectra and their respective
    # uncertainties
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
        dpix = i + axdict["wavelength_start"]

        # If the current element belongs in the next bin as defined by getbins, use the new bin's
        # parameters and increment the bin number.
        if (
            bin_number < len(binparameters) - 1
            and dpix == binparameters[bin_number + 1][-2]
        ):
            bin_number += 1
            b = binparameters[bin_number]

        # Get the extraction limits for the current dispersion element and define the spatial axis.
        # Where the extraction limits include partial pixels on the edge of the aperture, those
        # pixels are included in their entirety.
        loextlim = extractionlimits[0][i]
        hiextlim = extractionlimits[1][i]
        ax = axdict["spatial_axis"][int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))]

        # Use the extraction limits to define the data column for this wavelength element.
        col = col[int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))]

        # Use the extraction limits to define the errs column for this wavelength element. Where
        # errs have a value of 0, set them to the median err value for the entire column.
        err = errs2D[i]
        err[np.where(err == 0)] = np.median(err)
        err = err[int(np.floor(loextlim)) : int(np.ceil(hiextlim + 1))]

        # Use the err column to get the variance column
        var = err * err

        # Perform the standard aperture extraction on the data in the column, and add this value to
        # the aperture 1D output array. Get the root sum square of the err column and add this
        # value to the uncertainty array for the 1D aperture extraction.
        f = np.sum(col)
        aperdata1D[i] += f
        apererrs1D[i] += np.sqrt(np.sum(var))

        # Step 5 of the Horne 1986 algorithm - Defining the spatial profile. Use the average value
        # of the extraction limits in the current column to estimate the location of the centre of
        # the peak of the spatial profile.
        profcent = (hiextlim + loextlim) * 0.5

        # Use the Moffat profile parameters for the current bin to make a moffat profile that
        # approximates the shape of the spectrum's spatial profile in the current column. The
        # estimated centre of the extraction limits is used to shift the profile to the correct
        # location along the spatial axis. The background is assumed to have been subtracted by
        # now, so the background level and gradient are set to 0. Because the profile is a model
        # PSF, none of its values are negative and positivity of the profile need not be enforced.
        profile = moffat(b[0], profcent, b[2], b[3], 0.0, 0.0, ax)
        # Normalize the profile such that its sum equals unity.
        nprof = profile / np.sum(profile)

        # Step 7 of Horne 1986 - Cosmic ray masking with a conservative 5 sigmas.
        value = np.abs((col - (f * nprof)) ** 2)
        sigma = 5
        condition = var * sigma**2
        # With invert, it gives False (ie. 0) if the above condition is satisfied which means to
        # mask the rays
        cosmic_mask = np.invert([value > condition])[0]

        # Add an exception to handle scenarios where the mask is all False If that happens, ignore
        # the cosmic ray correction.
        if not np.any(cosmic_mask):
            cosmic_mask = np.full(shape=cosmic_mask.shape, fill_value=True)

        # Extract the optimal spectrum (Step 8 of the algorithm) Page 4 of Horne 1986: These
        # formulae are equivalent to determining the OPT spectrum by scaling a known spatial
        # profile P, to fit the sky subtracted data, whose variances are V.

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
    """
    Plot the spatial profile of a collapsed spectrum or a collapsed bin therein, and plot the
    Moffat function fitted to the data on top.

    Args:
        spataxis (numpy.ndarray)      : The spatial, or x, axis of the profile.
        bindata (numpy.ndarray)       : The binned data that has been fitted with a Moffat profile.
        hiresspataxis (numpy.ndarray) : The supersampled spatial axis used only for plotting
                                        purposes.
        binmoffparams (list)          : Parameters defining the Moffat profiel that was fitted to
                                        bindata.
        imgstart (int)                : The limit of the spatial axis after the original 2D
                                        spectrum was cut down to the region defined in reg.txt
        headparams (dict)             : A dictionary of parameters pulled from the header of the
                                        current datafile

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
    plt.ylabel("Median Flux, " + headparams["flux_unit"])
    plt.show()

    return None


def print_moffat_parameters(moffparams, imgstart, datascale):
    """
    Takes a list of Moffat profile parameters, the lower limit of the spatial axis after the
    original 2D spectrum was cut down to the region defined in reg.txt and the multiplier used to
    scale the spatial profile so it could be fit with a Moffat profile using scipy least_squares
    and prints them to the terminal.

    Args:
        moffparams (list) : The parameters of the fitted Moffat profile.
        imgstart (int)    : The lower limit of the spatial axis after the original 2D spectrum was
                            cut down to the region defined in reg.txt.
        datascale (float) : The multiplier used to scale the spatial profile so it could be fit
                            with a Moffat profile using scipy least_squares.

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


def show_img(data2D, axdict, headparams, drawlines, title):
    """
    Takes an input image and line data to be drawn on that image and creates a figure to be shown
    on screen.

    Args:
        data2D (numpy.ndarray) : The input image.
        axdict (dict)          : A dictionary containing the spatial and spectral axes of the input
                                 image.
        headparams (dict)      : A dictionary containing the header parameters of the input image.
        drawlines (list)       : A list of line data to be drawn on the input image.
        title (str)            : The title of the figure.

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
                axdict["wavelength_start"],
                axdict["wavelength_start"] + len(axdict["wavelength_axis"]),
                axdict["spatial_axis"][0] + axdict["data_spatial_floor"],
                axdict["spatial_axis"][-1] + axdict["data_spatial_floor"],
            ],
        )

        for i in range(int(len(drawlines) / 2)):
            ax.plot(drawlines[i * 2], drawlines[(i * 2) + 1], color="white")
        cbar = fig.colorbar(s, cax=colax)
        cbar.ax.yaxis.set_offset_position("left")
        cbar.ax.set_ylabel(
            "Pixel Flux, x10^" + str(power) + " " + headparams["flux_unit"]
        )
        ax2 = ax.twiny()
        ax2.plot(axdict["wavelength_axis"], data2D[0, :], alpha=0)
        ax2.set_xlim(axdict["wavelength_axis"][0], axdict["wavelength_axis"][-1])
        ax2.set_xlabel("Wavelength, " + headparams["wavelength_unit"])
        ax.set_ylim(
            axdict["spatial_axis"][0] + axdict["data_spatial_floor"],
            axdict["spatial_axis"][-1] + axdict["data_spatial_floor"],
        )
        ax.set_ylabel("Spatial Axis, Pixels")
        ax.set_xlim(axdict["wavelength_start"], axdict["wavelength_start"] + len(axdict["wavelength_axis"]))

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
            valinit=np.nanmedian(masked_data2D) + (3.0 * np.nanstd(masked_data2D)),
            valstep=0.001 * np.nanmax(masked_data2D),
        )

        def update(val):
            vmax = smax.val
            smin.valmax = vmax - 1
            vmin = smin.val
            smax.valmin = vmin + 1

            s.set_clim(vmin, vmax)
            cbar.ax.set_ylabel(
                "Pixel Flux, x10^" + str(power) + " " + headparams["flux_unit"]
            )
            fig.canvas.draw_idle()

        smin.on_changed(update)
        smax.on_changed(update)

        plt.show()

    return None


def subtract_sky(bglowext, bghighext, fdict, axdict, pars, hpars):
    """
    Subtracts the sky background from the 2D image by defining bg regions using limits input to the
    function and then fitting a profile to the background column by column while masking cosmic
    rays. The background level or profile for each column is subtracted from the full column to
    produce a background subtracted 2D image.

    Args:
        bglowext (numpy.ndarray)  : Lower limits of the background regions.
        bghighext (numpy.ndarray) : Upper limits of the background regions.
        fdict (dict)              : A dictionary containing the 2D image.
        axdict (dict)             : A dictionary containing the axis information.
        pars (dict)               : A dictionary containing MOTES parameters.
        hpars (dict)              : A dictionary containing header information.

    Returns:
        fdict (dict) : Dictionary containing the background subtracted 2D image.
    """

    sys.stdout.write("Subtracting sky background...\n")

    fdict["data"] = fdict["data"].T
    fdict["errs"] = fdict["errs"].T

    medsky = []
    colnum = len(bglowext)

    # Makes sure the limits are within the image.
    for ii in range(colnum):
        if bglowext[ii] < 0:
            bglowext[ii] = 0
        if bghighext[ii] > axdict["spatial_axis"][-1]:
            bghighext[ii] = axdict["spatial_axis"][-1]

        datacol = fdict["data"][ii]
        colrange = np.array(range(len(datacol)))
        skypix = datacol[
            np.where(np.logical_or(colrange < bglowext[ii], colrange > bghighext[ii]))
        ]

        # Kill MOTES if there isn't enough background sky to perform sky subtraction. Should
        # probably be made more nuanced later.
        if len(skypix) == 0:
            sys.stdout.write(" >>> No pixels contained inside sky region.\n")
            sys.stdout.write(
                "     -BG_FWHM_MULTIPLIER in motesparams.txt is probably too large.\n"
            )
            sys.stdout.write("     Please reduce -BG_FWHM_MULTIPLIER and try again.\n")
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
            np.where(np.logical_or(colrange < bglowext[ii], colrange > bghighext[ii]))
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
            skylevel = (skyquad * colrange * colrange) + (skygrad * colrange) + skyint
            medsky.append(skylevel)
            skyerr = (
                (skyquaderr * skyquaderr * colrange * colrange * colrange * colrange)
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
                + (skyquaderr * skyquaderr * colrange * colrange * colrange * colrange)
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
        fdict["sky_model"] = np.tile(medsky, (np.shape(fdict["data"])[1], 1))
    else:
        fdict["sky_model"] = medsky

    fdict["data"] = fdict["data"].T
    fdict["errs"] = fdict["errs"].T

    chipgaps = np.where(np.nanmedian(fdict["data"], axis=0) == 0)
    fdict["data"][:, chipgaps[0]] += 1

    return fdict
