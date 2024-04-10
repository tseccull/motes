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


def get_bins_output(bin_parameters, parameters, spatial_lo_limit, spatial_hi_limit, data_2D, header_parameters, axes_dict):
    """
    Print and plot the output of get_bins.

    Args:
        bin_parameters (list)       : A list containing the locations and S/N of each bin
        parameters (dict)          : A dictionary of parameters read in from the motesparams.txt
                                 configuration file.
        spatial_lo_limit (int)           : lower limit of the spatial region measured for S/N in get_bins()
        spatial_hi_limit (int)          : upper limit of the spatial region measured for S/N in get_bins()
        data_2D (numpy.ndarray) : 2D spectroscopic data
        header_parameters (dict)      : A dictionary of parameters full from the datafile header.
        axes_dict (dict)          : A dictionary containing axes and axis metadata for the current
                                 extraction

    Returns:
        None
    """

    sys.stdout.write(
        " >>> "
        + str(len(bin_parameters))
        + " spectrum localisation bins determined on dispersion axis.\n"
    )

    # DIAGNOSTICS - Plot boundary locations of localisation bins on the dispersion axis.
    if parameters["-DIAG_PLOT_BIN_LOC"]:
        draw_lines = []
        for b in bin_parameters:
            bin_line_location = np.where(
                np.logical_and(
                    axes_dict["spatial_axis"] > spatial_lo_limit - ((spatial_hi_limit - spatial_lo_limit) * 0.2),
                    axes_dict["spatial_axis"] < spatial_hi_limit + ((spatial_hi_limit - spatial_lo_limit) * 0.2),
                )
            )
            draw_lines.append(
                np.ones(len(axes_dict["spatial_axis"][bin_line_location])) * b[0] + axes_dict["wavelength_start"]
            )
            draw_lines.append(axes_dict["spatial_axis"][bin_line_location] + axes_dict["data_spatial_floor"])
            draw_lines.append(
                np.ones(len(axes_dict["spatial_axis"][bin_line_location])) * b[1] + axes_dict["wavelength_start"]
            )
            draw_lines.append(axes_dict["spatial_axis"][bin_line_location] + axes_dict["data_spatial_floor"])

        show_img(
            data_2D,
            axes_dict,
            header_parameters,
            draw_lines,
            "2D Spectrum with Boundaries of Localisation Bins",
        )

    return None


def interpolate_extraction_lims(extraction_limits, dispersion_axis_length):
    """
    Interpolate extraction limits over unbinned 2D spectrum to create the final extraction limits
    for the unbinned spectrum. Takes an input of extraction limits from the fitting of the binned
    data and interpolates the limits over the unbinned data. Limits are also linearly extrapolated
    towards the ends of the spectral range.

    Args:
        extraction_limits (numpy.ndarray) : Extraction limits for each bin.
        dispersion_axis_length (float)            : Length of the dispersion axis of the unbinned 2D spectrum.

    Returns:
        final_extraction_limits (list)            : List containing the extraction limits for the unbinned 2D
                                         spectrum.
    """

    # If the 2D spectrum was so faint that only 1 dispersion bin could be determined, set the
    # extraction limits across the unbinned spectrum to be a simple linear aperture that has the
    # same extraction limits as was determined for that one bin.
    if len(extraction_limits[0]) == 1:
        final_extraction_limits = [
            np.repeat(extraction_limits[1][0], dispersion_axis_length),
            np.repeat(extraction_limits[2][0], dispersion_axis_length),
        ]

    # Otherwise, interpolate the extraction limits form the bins across the unbinned wavelength
    # axis. Also extrapolate the extraction limits in a linear fashion at the ends of the
    # wavelength axis of the 2D spectrum so that the full wavelength axis is covered.
    else:
        interpolate_extraction_limits_1 = interp.interp1d(
            extraction_limits[0], extraction_limits[1], kind="linear"
        )
        interpolate_extraction_limits_2 = interp.interp1d(
            extraction_limits[0], extraction_limits[2], kind="linear"
        )

        intermediate_extraction_limits = [
            interpolate_extraction_limits_1(
                np.array(np.arange(extraction_limits[0][0], extraction_limits[0][-1]))
            ),
            interpolate_extraction_limits_2(
                np.array(np.arange(extraction_limits[0][0], extraction_limits[0][-1]))
            ),
        ]

        (
            short_extrapolated_limit_1,
            short_extrapolated_limit_2,
            long_extrapolated_limit_1,
            long_extrapolated_limit_2,
        ) = extrapolate_extraction_limits(
            intermediate_extraction_limits,
            dispersion_axis_length,
            extraction_limits[0][0],
            extraction_limits[0][-1],
        )

        extraction_limit_1 = np.insert(intermediate_extraction_limits[0], 0, short_extrapolated_limit_1)
        extraction_limit_1 = np.append(extraction_limit_1, long_extrapolated_limit_1)

        extraction_limit_2 = np.insert(intermediate_extraction_limits[1], 0, short_extrapolated_limit_2)
        extraction_limit_2 = np.append(extraction_limit_2, long_extrapolated_limit_2)

        new_extraction_limits = [extraction_limit_1, extraction_limit_2]
        interpolation_x_axis = np.insert(
            np.arange(extraction_limits[0][0], extraction_limits[0][-1]), 0, 0
        )
        interpolation_x_axis = np.append(interpolation_x_axis, dispersion_axis_length)

        interpolated_new_extraction_limit_1 = interp.interp1d(
            interpolation_x_axis, new_extraction_limits[0], kind="linear", fill_value="extrapolate"
        )
        interpolated_new_extraction_limit_2 = interp.interp1d(
            interpolation_x_axis, new_extraction_limits[1], kind="linear", fill_value="extrapolate"
        )

        final_extraction_limits = [
            interpolated_new_extraction_limit_1(np.array(range(dispersion_axis_length))),
            interpolated_new_extraction_limit_2(np.array(range(dispersion_axis_length))),
        ]

    return final_extraction_limits


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


def linear_resid(x, data_range, data):
    """
    Calculate residuals of fitted linear profile and the data for the Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the linear
                                    profile.
        data_range (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray)  : The residuals of the fitted line and the data.
    """
    residual = (x[0] * data_range) + x[1] - data
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


def poly2_resid(x, data_range, data):
    """
    Calculates residuals of a fitted second-order polynomial and the data for the
    Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the polynomial.
        data_range (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray) : The residuals of the fitted polynominal and the data.
    """
    residual = (x[0] * data_range * data_range) + (x[1] * data_range) + x[2] - data

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


def poly3_resid(x, data_range, data):
    """
    Calculates residuals of a fitted third-order polynomial and the data for the
    Levenberg-Marquardt least squares method.

    Args:
        x (list)                  : A list containing the best fit parameters of the polynomial.
        data_range (numpy.ndarray) : the spatial axis of the data column.
        data (numpy.ndarray)      : the data column.

    Returns:
        residual (numpy.ndarray) : The residuals of the fitted polynominal and the data.
    """
    residual = (
        (x[0] * data_range * data_range * data_range)
        + (x[1] * data_range * data_range)
        + (x[2] * data_range)
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


def moffat(amplitude, center, alpha, beta, background_level, background_gradient, data_range):
    """
    Creates moffat profile added a linear sloped background based on input parameters

    Args:
        amplitude (float64)             : amplitude of the Moffat profile
        center (float64)               : location of the center/peak of the Moffat profile on the
                                    spatial axis.
        alpha (float64)           : he main parameter that defines the width of the Moffat profile.
        beta (float64)            : plays a role in defining the width of the Moffat profile in the
                                    outer wings far from the profile's peak.
        background_level (float64)         : height of the linearly varying background
        background_gradient (float64)          : gradient of the linearly varying background
        data_range (numpy.ndarray) : x axis of the Moffat profile.

    Returns:
        moffat_background (numpy.ndarray) : a moffat profile defined on the data_range axis using
                                            the parameters input to the function, with added
                                            background flux.
    """

    moffat = amplitude * (
        (1 + (((data_range - center) * (data_range - center)) / (alpha * alpha))) ** -beta
    )
    moffat_background = moffat + background_level + (data_range * background_gradient)

    return moffat_background


def moffat_least_squares(r, col, seeing, pixel_resolution):
    """
    Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to
    the column using a least squares method. Returns the best fit parameters of the Moffat
    function.

    Args:
        r (numpy.ndarray)   : spatial axis of the data being fit
        col (numpy.ndarray) : data being fitted
        seeing (float)      : estimated FWHM of the spatial profile
        pixel_resolution (float)      : spatial resolution of each pixel in arcsec/pixel

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
        seeing / pixel_resolution,
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
                (5 * seeing / pixel_resolution),
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


def moffat_resid(x, data_range, data):
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
        data_range (numpy.ndarray) : spatial axis of the data
        data (numpy.ndarray)      : the data

    Returns:
        residual (numpy.ndarray) : the residual array between the model moffat profile and the data
    """

    moff = x[0] * (
        (1 + ((data_range - x[1]) * (data_range - x[1])) / (x[2] * x[2])) ** -x[3]
    )
    residual = moff + x[4] + (data_range * x[5]) - data
    return residual


def optimal_extraction(data_2D, errs_2D, extraction_limits, bin_parameters, axes_dict):
    """
    Perform optimal extraction using a modified version of Horne (1986) where S=0, G=0 and errors
    are not 'revised' since we already have the 'variance frame'. Ideally, this extraction reduces
    errors with the weighting and conserves the flux (when compared to the standard extraction).
    This subroutine uses an analytic Moffat profile (post sky-subtraction), instead of polynominals
    as in Horne (1986). Furthermore, profiling takes place bin by bin, accounting for spatial
    profile variations across the dispersion axis; extraction takes place column by column within
    the bin limits and within the extraction limits previously calculated.

    Args:
        data_2D (numpy.ndarray)           : Input data frame
        errs_2D (numpy.ndarray)           : Input error frame
        extraction_limits (numpy.ndarray) : An array containing limits at each dispersion pixel
        bin_parameters (list)             : A list containing the bin limits across the dispersion
                                           axis, to enable slicing the data across dispersion axis.
        axes_dict (dict)                    : A dictionary containing the spatial axis array and other
                                           relevant information about the size and shape of the
                                           data frame

    Returns:
        optimal_1D_data (numpy.ndarray) : 1D array of the optimally extracted spectrum
        optimal_1D_errs (numpy.ndarray) : 1D array of the uncertainties of the optimally extracted
                                     spectrum
        aperture_1D_data (numpy.ndarray) : 1D array of the aperture extracted spectrum
        aperture_1D_errs (numpy.ndarray) : 1D array of the uncertainties of the aperture extracted
                                     spectrum
    """

    sys.stdout.write("Performing optimal extraction...")

    # Filter any NaNs and Inf for data/errs AND ensure the errors are positive for this extraction.
    data_2D, errs_2D = filter_data(data_2D, np.abs(errs_2D))

    # Set up output arrays for the optimally and aperture extracted spectra and their respective
    # uncertainties
    optimal_1D_data = np.zeros(np.shape(data_2D)[0])
    optimal_1D_errs = np.zeros(np.shape(data_2D)[0])
    aperture_1D_data = np.zeros(np.shape(data_2D)[0])
    aperture_1D_errs = np.zeros(np.shape(data_2D)[0])

    # Set up bin identification parameters for the following loop.
    bin_number = -1
    b = np.zeros(8)

    # Loop through each dispersion element of the spectrum.
    for i, col in enumerate(data_2D):
        # Identify the location of the current element in the original 2D spectrum
        og_column_location = i + axes_dict["wavelength_start"]

        # If the current element belongs in the next bin as defined by getbins, use the new bin's
        # parameters and increment the bin number.
        if (
            bin_number < len(bin_parameters) - 1
            and og_column_location == bin_parameters[bin_number + 1][-2]
        ):
            bin_number += 1
            b = bin_parameters[bin_number]

        # Get the extraction limits for the current dispersion element and define the spatial axis.
        # Where the extraction limits include partial pixels on the edge of the aperture, those
        # pixels are included in their entirety.
        lo_extraction_limit = extraction_limits[0][i]
        hi_extraction_limit = extraction_limits[1][i]
        ax = axes_dict["spatial_axis"][int(np.floor(lo_extraction_limit)) : int(np.ceil(hi_extraction_limit + 1))]

        # Use the extraction limits to define the data column for this wavelength element.
        col = col[int(np.floor(lo_extraction_limit)) : int(np.ceil(hi_extraction_limit + 1))]

        # Use the extraction limits to define the errs column for this wavelength element. Where
        # errs have a value of 0, set them to the median err value for the entire column.
        err = errs_2D[i]
        err[np.where(err == 0)] = np.median(err)
        err = err[int(np.floor(lo_extraction_limit)) : int(np.ceil(hi_extraction_limit + 1))]

        # Use the err column to get the variance column
        var = err * err

        # Perform the standard aperture extraction on the data in the column, and add this value to
        # the aperture 1D output array. Get the root sum square of the err column and add this
        # value to the uncertainty array for the 1D aperture extraction.
        f = np.sum(col)
        aperture_1D_data[i] += f
        aperture_1D_errs[i] += np.sqrt(np.sum(var))

        # Step 5 of the Horne 1986 algorithm - Defining the spatial profile. Use the average value
        # of the extraction limits in the current column to estimate the location of the centre of
        # the peak of the spatial profile.
        profile_center = (hi_extraction_limit + lo_extraction_limit) * 0.5

        # Use the Moffat profile parameters for the current bin to make a moffat profile that
        # approximates the shape of the spectrum's spatial profile in the current column. The
        # estimated centre of the extraction limits is used to shift the profile to the correct
        # location along the spatial axis. The background is assumed to have been subtracted by
        # now, so the background level and gradient are set to 0. Because the profile is a model
        # PSF, none of its values are negative and positivity of the profile need not be enforced.
        profile = moffat(b[0], profile_center, b[2], b[3], 0.0, 0.0, ax)
        # Normalize the profile such that its sum equals unity.
        normalized_profile = profile / np.sum(profile)

        # Step 7 of Horne 1986 - Cosmic ray masking with a conservative 5 sigmas.
        value = np.abs((col - (f * normalized_profile)) ** 2)
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
            fopt = np.sum(((normalized_profile * col) / var)[cosmic_mask]) / np.sum(
                ((normalized_profile * normalized_profile) / var)[cosmic_mask]
            )
            vopt = np.sum(normalized_profile[cosmic_mask]) / np.sum(
                ((normalized_profile * normalized_profile) / var)[cosmic_mask]
            )
        optimal_1D_data[i] += fopt
        optimal_1D_errs[i] += np.sqrt(vopt)

    return optimal_1D_data, optimal_1D_errs, aperture_1D_data, aperture_1D_errs


def plot_fitted_spatial_profile(
    spatial_axis, bin_data, hi_resolution_spatial_axis, bin_moffat_parameters, image_start, header_parameters
):
    """
    Plot the spatial profile of a collapsed spectrum or a collapsed bin therein, and plot the
    Moffat function fitted to the data on top.

    Args:
        spatial_axis (numpy.ndarray)      : The spatial, or x, axis of the profile.
        bin_data (numpy.ndarray)       : The binned data that has been fitted with a Moffat profile.
        hi_resolution_spatial_axis (numpy.ndarray) : The supersampled spatial axis used only for plotting
                                        purposes.
        bin_moffat_parameters (list)          : Parameters defining the Moffat profiel that was fitted to
                                        bin_data.
        image_start (int)                : The limit of the spatial axis after the original 2D
                                        spectrum was cut down to the region defined in reg.txt
        header_parameters (dict)             : A dictionary of parameters pulled from the header of the
                                        current datafile

    Returns:
        None
    """

    plt.figure(figsize=(7, 4.5))
    plt.plot(
        hi_resolution_spatial_axis + image_start,
        moffat(
            bin_moffat_parameters[0],
            bin_moffat_parameters[1],
            bin_moffat_parameters[2],
            bin_moffat_parameters[3],
            bin_moffat_parameters[4],
            bin_moffat_parameters[5],
            hi_resolution_spatial_axis,
        ),
        color="r",
        linewidth=3,
        label="Fitted Moffat Profile",
    )

    plt.plot(spatial_axis + image_start, bin_data, color="k", label="Spatial Profile")
    plt.grid(linestyle="dashed", color="gray")
    plt.legend()
    plt.title("Spectrum Spatial Profile and Fitted Moffat Profile")
    plt.xlabel("Spatial Axis, Pixels")
    plt.ylabel("Median Flux, " + header_parameters["flux_unit"])
    plt.show()

    return None


def print_moffat_parameters(moffat_parameters, image_start, data_scale):
    """
    Takes a list of Moffat profile parameters, the lower limit of the spatial axis after the
    original 2D spectrum was cut down to the region defined in reg.txt and the multiplier used to
    scale the spatial profile so it could be fit with a Moffat profile using scipy least_squares
    and prints them to the terminal.

    Args:
        moffat_parameters (list) : The parameters of the fitted Moffat profile.
        image_start (int)    : The lower limit of the spatial axis after the original 2D spectrum was
                            cut down to the region defined in reg.txt.
        data_scale (float) : The multiplier used to scale the spatial profile so it could be fit
                            with a Moffat profile using scipy least_squares.

    Returns:
        None
    """
    sys.stdout.write(" >>> Fitted Moffat function parameters:\n")
    sys.stdout.write("         A = " + str(moffat_parameters[0]) + "\n")
    sys.stdout.write("         c = " + str(moffat_parameters[1] + image_start) + "\n")
    sys.stdout.write("     alpha = " + str(moffat_parameters[2]) + "\n")
    sys.stdout.write("      beta = " + str(moffat_parameters[3]) + "\n")
    sys.stdout.write("         B = " + str(moffat_parameters[4]) + "\n")
    sys.stdout.write("         m = " + str(moffat_parameters[5]) + "\n\n")
    sys.stdout.write(
        " >>> Profile scaling factor used for fitting: " + str(data_scale) + "\n"
    )
    sys.stdout.write(
        " >>> Plot of median spatial profile presents the orginal unscaled profile.\n"
    )

    return None


def show_img(data_2D, axes_dict, header_parameters, draw_lines, title):
    """
    Takes an input image and line data to be drawn on that image and creates a figure to be shown
    on screen.

    Args:
        data_2D (numpy.ndarray) : The input image.
        axes_dict (dict)          : A dictionary containing the spatial and spectral axes of the input
                                 image.
        header_parameters (dict)      : A dictionary containing the header parameters of the input image.
        draw_lines (list)       : A list of line data to be drawn on the input image.
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

        power = int(np.floor(np.log10(np.abs(np.nanmean(data_2D)))))
        data_2D = copy.deepcopy(data_2D) / 10**power

        figwidth = 10.0
        fig = plt.figure(figsize=(figwidth, figwidth / 1.9))
        gs = gridspec.GridSpec(18, 33)
        ax = plt.subplot(gs[:, :32])
        color_axis = plt.subplot(gs[1:, 32])
        masked_data2D = np.ma.masked_where(data_2D == 0, data_2D)
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
                axes_dict["wavelength_start"],
                axes_dict["wavelength_start"] + len(axes_dict["wavelength_axis"]),
                axes_dict["spatial_axis"][0] + axes_dict["data_spatial_floor"],
                axes_dict["spatial_axis"][-1] + axes_dict["data_spatial_floor"],
            ],
        )

        for i in range(int(len(draw_lines) / 2)):
            ax.plot(draw_lines[i * 2], draw_lines[(i * 2) + 1], color="white")
        cbar = fig.colorbar(s, cax=color_axis)
        cbar.ax.yaxis.set_offset_position("left")
        cbar.ax.set_ylabel(
            "Pixel Flux, x10^" + str(power) + " " + header_parameters["flux_unit"]
        )
        ax2 = ax.twiny()
        ax2.plot(axes_dict["wavelength_axis"], data_2D[0, :], alpha=0)
        ax2.set_xlim(axes_dict["wavelength_axis"][0], axes_dict["wavelength_axis"][-1])
        ax2.set_xlabel("Wavelength, " + header_parameters["wavelength_unit"])
        ax.set_ylim(
            axes_dict["spatial_axis"][0] + axes_dict["data_spatial_floor"],
            axes_dict["spatial_axis"][-1] + axes_dict["data_spatial_floor"],
        )
        ax.set_ylabel("Spatial Axis, Pixels")
        ax.set_xlim(axes_dict["wavelength_start"], axes_dict["wavelength_start"] + len(axes_dict["wavelength_axis"]))

        ax.set_xlabel("Dispersion Axis, Pixels")

        plt.title(title, y=1.095)

        # Add interactive scaling bar to figures if the 2D spectrum isn't flux calibrated.
        # Again for some reason teeny tiny numbers cause things to break.
        fig.subplots_adjust(bottom=0.2)
        ax_vmin = plt.axes([0.1, 0.05, 0.8, 0.03])
        ax_vmax = plt.axes([0.1, 0.01, 0.8, 0.03])
        slider_min = Slider(
            ax_vmin,
            "LowCut",
            0,
            np.nanmax(masked_data2D) - 1,
            valinit=0.0,
            valstep=0.001 * np.nanmax(masked_data2D),
        )
        slider_max = Slider(
            ax_vmax,
            "HighCut",
            1.0,
            np.nanmax(masked_data2D),
            valinit=np.nanmedian(masked_data2D) + (3.0 * np.nanstd(masked_data2D)),
            valstep=0.001 * np.nanmax(masked_data2D),
        )

        def update(val):
            vmax = slider_max.val
            slider_min.valmax = vmax - 1
            vmin = slider_min.val
            slider_max.valmin = vmin + 1

            s.set_clim(vmin, vmax)
            cbar.ax.set_ylabel(
                "Pixel Flux, x10^" + str(power) + " " + header_parameters["flux_unit"]
            )
            fig.canvas.draw_idle()

        slider_min.on_changed(update)
        slider_max.on_changed(update)

        plt.show()

    return None


def subtract_sky(background_spatial_lo_limit, background_spatial_hi_limit, frame_dict, axes_dict, parameters, header_parameters):
    """
    Subtracts the sky background from the 2D image by defining bg regions using limits input to the
    function and then fitting a profile to the background column by column while masking cosmic
    rays. The background level or profile for each column is subtracted from the full column to
    produce a background subtracted 2D image.

    Args:
        background_spatial_lo_limit (numpy.ndarray)  : Lower limits of the background regions.
        background_spatial_hi_limit (numpy.ndarray) : Upper limits of the background regions.
        frame_dict (dict)              : A dictionary containing the 2D image.
        axes_dict (dict)             : A dictionary containing the axis information.
        parameters (dict)               : A dictionary containing MOTES parameters.
        header_parameters (dict)              : A dictionary containing header information.

    Returns:
        frame_dict (dict) : Dictionary containing the background subtracted 2D image.
    """

    sys.stdout.write("Subtracting sky background...\n")

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    sky_model = []
    number_of_columns = len(background_spatial_lo_limit)

    # Makes sure the limits are within the image.
    for ii in range(number_of_columns):
        if background_spatial_lo_limit[ii] < 0:
            background_spatial_lo_limit[ii] = 0
        if background_spatial_hi_limit[ii] > axes_dict["spatial_axis"][-1]:
            background_spatial_hi_limit[ii] = axes_dict["spatial_axis"][-1]

        data_column = frame_dict["data"][ii]
        column_axis = np.array(range(len(data_column)))
        sky_pixels = data_column[
            np.where(np.logical_or(column_axis < background_spatial_lo_limit[ii], column_axis > background_spatial_hi_limit[ii]))
        ]

        # Kill MOTES if there isn't enough background sky to perform sky subtraction. Should
        # probably be made more nuanced later.
        if len(sky_pixels) == 0:
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
                                np.shape(frame_dict["data"])[0],
                                np.shape(frame_dict["data"])[1],
                            ]
                        )
                        / (2 * header_parameters["seeing"]),
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

        sky_axis = column_axis[
            np.where(np.logical_or(column_axis < background_spatial_lo_limit[ii], column_axis > background_spatial_hi_limit[ii]))
        ]
        if len(set(sky_pixels)) == 1:
            continue

        else:
            sky_pixels_median = np.nanmedian(sky_pixels)
            sky_pixels_sigma = np.nanstd(sky_pixels)
            good_sky_pixels_location = np.where(
                np.logical_and(
                    sky_pixels > sky_pixels_median - (10 * sky_pixels_sigma),
                    sky_pixels < sky_pixels_median + (10 * sky_pixels_sigma),
                )
            )
            good_sky_pixels = good_sky_pixels[good_sky_pixels_location]
            sky_axis = sky_axis[good_sky_pixels_location]

        if parameters["-SKYSUB_MODE"] == "MEDIAN":
            bootsky = np.random.choice(good_sky_pixels, (len(good_sky_pixels), 100), replace=True)
            skysamp = np.nanmedian(bootsky, axis=0)
            skylevel = np.nanmean(skysamp)
            sky_model.append(skylevel)
            skyerr = np.std(skysamp) / (99**0.5)

        if parameters["-SKYSUB_MODE"] == "LINEAR":
            bootsky = np.random.choice(good_sky_pixels, (len(good_sky_pixels), 100), replace=True)
            grads = []
            intercepts = []
            for jj in bootsky.T:
                linpars = linear_least_squares(sky_axis, jj)
                grads.append(linpars[0])
                intercepts.append(linpars[1])
            intercepts = np.array(intercepts)
            grads = np.array(grads)
            skygrad = np.mean(grads)
            skygraderr = np.std(grads) / (99**0.5)
            skyint = np.mean(intercepts)
            skyinterr = np.std(intercepts) / (99**0.5)
            skylevel = (skygrad * column_axis) + skyint
            sky_model.append(skylevel)
            skyerr = (
                (skygraderr * column_axis * skygraderr * column_axis)
                + (skyinterr * skyinterr)
            ) ** 0.5

        if parameters["-SKYSUB_MODE"] == "POLY2":
            bootsky = np.random.choice(good_sky_pixels, (len(good_sky_pixels), 100), replace=True)
            grads = []
            intercepts = []
            quads = []
            for jj in bootsky.T:
                linpars = poly2_least_squares(sky_axis, jj)
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
            skylevel = (skyquad * column_axis * column_axis) + (skygrad * column_axis) + skyint
            sky_model.append(skylevel)
            skyerr = (
                (skyquaderr * skyquaderr * column_axis * column_axis * column_axis * column_axis)
                + (skygraderr * column_axis * skygraderr * column_axis)
                + (skyinterr * skyinterr)
            ) ** 0.5

        if parameters["-SKYSUB_MODE"] == "POLY3":
            bootsky = np.random.choice(good_sky_pixels, (len(good_sky_pixels), 100), replace=True)
            grads = []
            intercepts = []
            quads = []
            trips = []
            for jj in bootsky.T:
                linpars = poly3_least_squares(sky_axis, jj)
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
                (skytrip * column_axis * column_axis * column_axis)
                + (skyquad * column_axis * column_axis)
                + (skygrad * column_axis)
                + skyint
            )
            sky_model.append(skylevel)
            skyerr = (
                (
                    skytriperr
                    * skytriperr
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                )
                + (skyquaderr * skyquaderr * column_axis * column_axis * column_axis * column_axis)
                + (skygraderr * column_axis * skygraderr * column_axis)
                + (skyinterr * skyinterr)
            ) ** 0.5

        frame_dict["data"][ii] -= skylevel
        frame_dict["errs"][ii] = ((frame_dict["errs"][ii] ** 2) + (skyerr**2)) ** 0.5

        sys.stdout.write(
            "     " + str(ii + 1) + "/" + str(number_of_columns) + " columns completed.\r"
        )

    sky_model = np.array(sky_model)
    if parameters["-SKYSUB_MODE"] == "MEDIAN":
        frame_dict["sky_model"] = np.tile(sky_model, (np.shape(frame_dict["data"])[1], 1))
    else:
        frame_dict["sky_model"] = sky_model

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    chipgaps = np.where(np.nanmedian(frame_dict["data"], axis=0) == 0)
    frame_dict["data"][:, chipgaps[0]] += 1

    return frame_dict
