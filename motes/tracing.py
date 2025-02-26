#!/usr/bin/env python3

import logging
import motes.common as common
import numpy as np
import scipy.interpolate as interp

from scipy.optimize import least_squares

logger = logging.getLogger("motes")


def end_caps(start, middle, end):
    '''
    Combines three array-likes into an array.
    
    Args:
    -- start (numpy array-like: numpy.ndarray, int)
         Int or array to prepend to the start of the middle array.
    -- middle (numpy.ndarray)
         Array in the middle of the output array.
    -- end (numpy array-like: numpy.ndarray, int)
         Int or array to append to the end of the middle array.
    
    Returns:
     -- capped_array (numpy.ndarray)
          Array produced by combining the three provided as args.
    '''
	
    capped_array = np.insert(middle, 0, start)
    capped_array = np.append(capped_array, end)
    
    return capped_array


def extrapolate_extraction_limits(
    extraction_limits,
    dispersion_axis_length,
    short_end_distance,
    long_end_distance
):
    """
    Linearly extrapolate the extraction limits at the ends of the 2D
    spectrum.

    Args:
     -- extraction_limits (list)
          A list containing the lower and upper extraction limits for
          each spatial pixel.
     -- dispersion_axis_length (int)
          The length of the dispersion axis.
     -- short_end_distance (int)
          The number of pixels at the short end of the dispersion axis
          to be excluded from the extraction.
     -- long_end_distance (int)
          The number of pixels at the long end of the dispersion axis to
          be excluded from the extraction.

    Returns:
     -- extrapolated_limits (dict)
          List of extrapolated sections of extraction limits.

    """
    extrap_grad_short_1 = extrapolate_gradient(extraction_limits[0], [0, 150, 300])
    extrap_grad_short_2 = extrapolate_gradient(extraction_limits[1], [0, 150, 300])
    extrap_grad_long_1 = extrapolate_gradient(extraction_limits[0], [-300, -150, -1])
    extrap_grad_long_2 = extrapolate_gradient(extraction_limits[1], [-300, -150, -1])

    long_end = dispersion_axis_length - long_end_distance
    extrap_limit_short_1 = extraction_limits[0][0] - (extrap_grad_short_1 * short_end_distance)
    extrap_limit_short_2 = extraction_limits[1][0] - (extrap_grad_short_2 * short_end_distance)
    extrap_limit_long_1 = extraction_limits[0][-1] + (extrap_grad_long_1 * long_end)
    extrap_limit_long_2 = extraction_limits[1][-1] + (extrap_grad_long_2 * long_end)
    
    extrapolated_limits = [
        extrap_limit_short_1,
        extrap_limit_short_2,
        extrap_limit_long_1,
        extrap_limit_long_2
    ]
    
    return extrapolated_limits


def extrapolate_gradient(data_array, data_limits):
    """
    Given a range of data and limits to define a region of that data,
    calculate the region's gradient.

    Args:
     -- data_array (numpy.ndarray)
          A 1D array of data.
     -- data_limits (list)
          A list containing the limits of the region of data to be used
          to calculate the gradient.

    Returns:
     -- gradient (int)
          The gradient of the region of data.
    """
    median_of_y_points_x_to_y = np.median(data_array[data_limits[0] : data_limits[1]])
    median_of_y_points_y_to_z = np.median(data_array[data_limits[1] : data_limits[2]])
    median_of_x_points_x_to_y = np.median([data_limits[0], data_limits[1]])
    median_of_x_points_y_to_z = np.median([data_limits[1], data_limits[2]])

    grad_numerator = median_of_y_points_x_to_y - median_of_y_points_y_to_z
    grad_denominator = median_of_x_points_x_to_y - median_of_x_points_y_to_z

    gradient = grad_numerator / grad_denominator

    return gradient


def get_bins(
    frame_dict,
    spatial_lo_limit,
    spatial_hi_limit,
    dispersion_axis_length,
    parameter_dict,
    has_sky=False
):
    """
    Define the bins of data over which Moffat profiles will be fitted.
    Each bin is defined such that when summed it will have a given
    signal to noise (S/N). So lower S/N regions will have larger bins.
    This function works from the centre of the input 2D spectrum outward
    to the ends in order to ensure a good start for the binning process.

    Args:
     -- frame_dict (dict)
          Dictionary containing all the data, error, and quality frames.
     -- spatial_lo_limit (int)
          Lower spatial limit of the region where the S/N will be
          measured when defining the extent of a bin.
     -- spatial_hi_limit (int)
          Upper spatial limit of the region where the S/N will be
          measured when defining the extent of a bin.
     -- dispersion_axis_length (int)
          Length of the dispersion axis.
     -- parameter_dict (dict)
          Dictionary of parameters ready in from the motesparams.txt
          configuration file.
     -- has_sky (bool, optional)
          Used to tell get_bins whether the sky has been subtracted yet
          or not, and to set the minimum_snr threshold accordingly.
          False by default.

    Returns:
     -- bin_locations (list)
          A list containing the details for each bin determined by
          get_bins. The boundaries and S/N of each bin are recorded
          here.
     -- frame_dict (dict)
          Returns frame_dict.
    """

    # Take S/N threshold (minimum_snr) and minimum number of columns per
    # dispersion bin (minimum_columns)
    if parameter_dict.subtract_sky and has_sky:
        minimum_snr = parameter_dict.sky_snr_limit
    else:
        minimum_snr = parameter_dict.extraction_snr_limit

    # Minimum columns for a bin
    minimum_columns = parameter_dict.minimum_column_limit

    # Start x at the central pixel column of the dispersion axis
    x = int(dispersion_axis_length / 2)
    width = 0
    bin_locations = []
    
    logger.info("Determining spectrum localisation bins on dispersion axis.")
    logger.info("User defined S/N threshold = %s.", minimum_snr)

    # Start at the centre of the dispersion axis and start binning
    # towards the short wavelength end of the spectrum.
    while x - width > 0:
        snr_estimate = 0.0

        # If the S/N of the current bin has not yet reached the
        # user-defined threshold (minimum_snr), add one more pixel
        # column to the bin.
        while snr_estimate <= minimum_snr:
            width += 1

            # Stop the loop once the short wavelength end of the
            # spectrum os reached.
            if x - width < 0:
                width = int(0 + x)
                break

            # If there aren't enough good pixels in each spatial column
            # of the current bin, continue to the next iteration and add
            # another pixel column to the bin.
            short_rows = get_short_rows(
                frame_dict["qual"], int(x-width), int(x), minimum_columns
            )
            
            if len(short_rows) > 0:
                continue

            # Determine the S/N of the bin.
            bin_data = frame_dict["data"][spatial_lo_limit:spatial_hi_limit, x-width:x]
            bin_errs = frame_dict["errs"][spatial_lo_limit:spatial_hi_limit, x-width:x]
            snr_estimate = get_bin_snr(bin_data, bin_errs)

        bin_locations.append([int(x - width), int(x), snr_estimate])

        x -= width
        width = 0

    x = int(dispersion_axis_length / 2)

    # Repeat the same process as above, starting at the centre of the
    # dispersion axis, but moving outward toward the longer wavelength
    # end of the 2D spectrum.
    while x + width < dispersion_axis_length:
        snr_estimate = 0.0

        # If the S/N of the current bin has not yet reached the
        # user-defined threshold (minimum_snr), add one more pixel
        # column to the bin.
        while snr_estimate <= minimum_snr:
            width += 1

            # Stop the loop once the long wavelength end of the spectrum
            # is reached.
            if x + width > dispersion_axis_length:
                width = int(dispersion_axis_length - x)
                break

            # If there aren't enough good pixels in each column of the 
            # current bin, continue to the next iteration and add
            # another pixel column to the bin.
            short_rows = get_short_rows(
                frame_dict["qual"], int(x), int(x + width), minimum_columns
            )
            if len(short_rows) > 0:
                continue

            # Determine the S/N of the bin.
            bin_data = frame_dict["data"][spatial_lo_limit:spatial_hi_limit, x:x+width]
            bin_errs = frame_dict["errs"][spatial_lo_limit:spatial_hi_limit, x:x+width]
            snr_estimate = get_bin_snr(bin_data, bin_errs)

        bin_locations.append([int(x), int(x + width), snr_estimate])

        x += width
        width = 0

    # Sort out the bins into the correct order on the dispersion axis.
    bin_locations.sort(key=lambda x: x[0])

    return bin_locations, frame_dict


def get_bin_snr(data, errs):
    """
    Calculates the signal-to-noise ratio of the supplied piece of data.
    
    Args:
     -- data (numpy.ndarray)
          A 2D data array that's summed to calculate the signal.
     -- errs (numpy.ndarray)
          A 2D array containing the uncertainies of the data array.
          The root-sum-square of the errs array is taken as the noise
          in the signal-to-noise calculation.
    
    Returns:
     -- snr_estimate (float)
          The estimated signal-to-noise ratio.
    """
    
    signal = np.nansum(data)
    root_sum_square_noise = np.sqrt(np.nansum(errs*errs))
    snr_estimate = signal / root_sum_square_noise
    
    return snr_estimate


def get_short_rows(qual, bin_bottom, bin_top, min_cols):
    """
    Determines whether any spatial rows in the current bin in get_bins
    have a number of good data points below the minimum number set by
    the user (min_cols).
    
    Args:
     -- qual (numpy.ndarray)
          2D array of quality flags associated with the current data bin
          in get_bins(). 1=GOOD, 0=BAD.
     -- bin_bottom (int)
          Index on the dispersion axis marking the lower limit of the
          data bin.
     -- bin_top (int)
          Index on the dispersion axis marking the upper limit of the
          data bin.
     -- min_cols (int)
          User-defined minimum number of data columns to include in a
          data bin in get_bins() before estimating its signal-to-noise.
    
    Returns:
     -- short_rows (list)
          List containing the number of good pixels present along the
          dispersion axis for spatial pixel rows where this number is
          below the minimum limit set by the user.
    """
    
    qual_dispersion_sum = np.sum(qual[:, bin_bottom : bin_top], axis=1)
    find_short_row = lambda x: x<= min_cols
    short_rows = list(filter(find_short_row, qual_dispersion_sum))
	
    return short_rows
    

def interpolate_extraction_lims(extraction_limits, dispersion_axis_length):
    """
    Interpolate extraction limits over unbinned 2D spectrum to create
    the final extraction limits for the unbinned spectrum. Takes an
    input of extraction limits from the fitting of the binned data and
    interpolates the limits over the unbinned data. Limits are also
    linearly extrapolated towards the ends of the spectral range.

    Args:
     -- extraction_limits (numpy.ndarray)
          Extraction limits for each bin.
     -- dispersion_axis_length (float)
          Length of the dispersion axis of the unbinned 2D spectrum.

    Returns:
     -- final_extraction_limits (list)
          List containing the extraction limits for the unbinned 2D
          spectrum.
    """

    # If the 2D spectrum was so faint that only 1 dispersion bin could
    # be determined, set the extraction limits across the unbinned
    # spectrum to be a simple linear aperture that has the same
    # extraction limits as was determined for that one bin.
    if len(extraction_limits[0]) == 1:
        final_extraction_limits = [
            np.repeat(extraction_limits[1][0], dispersion_axis_length),
            np.repeat(extraction_limits[2][0], dispersion_axis_length),
        ]

    # Otherwise, interpolate the extraction limits from the bins across
    # the unbinned wavelength axis. Also extrapolate the extraction
    # limits in a linear fashion at the ends of the wavelength axis of
    # the 2D spectrum so that the full wavelength axis is covered.
    else:
        interp_dispersion_axis = np.arange(extraction_limits[0][0], extraction_limits[0][-1])
        middle_extraction_limits = interpolation_loop(interp_dispersion_axis, extraction_limits)
        
        extrapolated_x_axis = end_caps(0, interp_dispersion_axis, dispersion_axis_length)
        extrapolated_limits = extrapolate_extraction_limits(
            middle_extraction_limits,
            dispersion_axis_length,
            extraction_limits[0][0],
            extraction_limits[0][-1],
        )
        
        new_extraction_limits = [extrapolated_x_axis]
        for i in range(1,3):
            capped_limits = end_caps(
                extrapolated_limits[i-1], 
                middle_extraction_limits[i-1], 
                extrapolated_limits[i+1]
            )
            new_extraction_limits.append(capped_limits)
        
        ext_lim_pixel_axis = np.array(range(dispersion_axis_length))
        final_extraction_limits = (
            interpolation_loop(ext_lim_pixel_axis, new_extraction_limits, filler="extrapolate")
        )
        
    return final_extraction_limits


def interpolation_loop(new_x_axis, bounds, filler=np.nan):
    """
    Interpolates the supplied extraction boundaries to a new dispersion
    pixel axis.
    
    Args:
     -- new_x_axis (numpy.ndarray)
          The new dispersion axis that the extraction limits will be
          interpolated to.
     -- bounds (numpy.ndarray)
          Extraction limits to be interpolated.
    
    Return:
     -- interpolated_boundaries (list)
          List of interpolated extraction limits.
    """
	
    interpolated_boundaries = [np.zeros(len(new_x_axis)), np.zeros(len(new_x_axis))]
    for i in range(2):
       interp_bound = interp.interp1d(bounds[0], bounds[i+1], kind="linear", fill_value=filler)
       interpolated_boundaries[i] += interp_bound(new_x_axis)
    
    return interpolated_boundaries


def median_moffat_fitting(median_profile, header_parameters, spatial_axis):
    """
    Fits a Moffat function to the median-collapsed spatial profile of
    of the 2D spectroscopic data. The parameters of the fitted Moffat
    profile and the data scaling factor used in the fit are returned.
    
    Args:
     -- median_profile (numpy.ndarray)
          Numpy array containing the median spatial profile of the 2D
          spectrum.
     -- header_parameters (dict)
          Dict of parameters harvested from the header of the input
          data file.
     -- spatial_axis (numpy.ndarray)
          Numpy array containing the spatial axis of the 2D
          spectroscopic data.
    
    Returns:
     -- moffat_profile_parameters (list)
          Parameters defining the shape of the Moffat profile fitted
          to the data.
     -- data_scaling_factor (numpy.float32)
          Factor used to scale median_profile so that its Moffat
          parameters can be fitted accurately. 
    """
    # Scipy least squares doesn't like really tiny numbers like fluxes 
    # in erg/s/cm^2/Angstrom, so it's necessary to scale the data to a
    # size that least squares can handle. The shape of the profile
    # fitted to the scaled spatial profile is the same as the unscaled,
    # but to to get a model profile that matches the original profile,
    # the profile amplitude (A), background level (B), and background 
    # gradient (m) all need to be scaled down again after the fitting.
    median_of_median_profile = np.nanmedian(median_profile)
    scaling_exponent = np.log10(np.abs(median_of_median_profile))
    data_scaling_factor = 10 ** np.abs(np.floor(scaling_exponent))
    # Fit the median spatial profile with a Moffat function.
    moffat_profile_parameters = moffat_least_squares(
        spatial_axis,
        median_profile * data_scaling_factor,
        header_parameters["seeing"],
        header_parameters["pixel_resolution"]
    )
    median_alpha = moffat_profile_parameters[2]
    median_beta = moffat_profile_parameters[3]
    # Get an improved estimate of the FWHM of the spectrum from the best 
    # fit Moffat profile.
    header_parameters["seeing"] = common.moffat_fwhm(median_alpha, median_beta)
    # Scale the amplitude, background gradient, and background level of
    # the model Moffat profile down.
    moffat_profile_parameters[0] /= data_scaling_factor
    moffat_profile_parameters[4] /= data_scaling_factor
    moffat_profile_parameters[5] /= data_scaling_factor
    
    return moffat_profile_parameters, data_scaling_factor


def moffat_fitting(
    each_bin,
    data,
    axes_dict,
    data_scaling_factor,
    header_parameters,
    fwhm_multiplier,
    extraction_limits
):
    """
    Fits a Moffat function to a median-collapsed bin of the 2D
    spectroscopic data defined by get_bins(). The parameters of the
    fitted Moffat profile and the extraction limits for the bin
    determined from those Moffat parameters are returned.
    
    Args:
     -- each_bin (list)
          List containing the dimensions of the current data bin being
          fit.
     -- data (numpy.ndarray)
          2D data frame that is the subsection of the full 2D spectrum
          frame demarked by the bin dimensions in each_bin.
     -- axes_dict (dict)
          A dictionary containing axes and axis metadata.
     -- data_scaling_factor (numpy.float32)
          Factor used to scale the data so that its Moffat parameters
          can be fitted accurately. 
     -- header_parameters (dict)
          Dict of parameters harvested from the header of the input
          data file.
     -- fwhm_multiplier (float)
          Factor used to scale median_profile so that its Moffat
          parameters can be fitted accurately. 
     -- extraction_limits (list)
          Extraction limits for each bin.
    
    Returns:
     -- bin_moffat_parameters (list)
          Parameters defining the shape of the Moffat profile fitted
          to the data.
     -- extraction_limits (list)
          Extraction limits for each bin.
    """

    # Take the median spatial profile of the dispersion bin, and leave
    # out pixel columns in the chip gaps if this is a GMOS spectrum.
    raw_bin_data = data[:, each_bin[0] : each_bin[1]]
    chip_gap_location = np.where(np.median(raw_bin_data, axis=0) != 1)
    bin_data = np.nanmedian(raw_bin_data[:, chip_gap_location[0]], axis=1)
    
    # Use a Levenberg-Marquardt Least Squares method to fit a Moffat
    # function to the median spatial profile and return its parameters.
    bin_moffat_parameters = moffat_least_squares(
        axes_dict["spatial_axis"],
        bin_data * data_scaling_factor,
        header_parameters["seeing"],
        header_parameters["pixel_resolution"],
    )
    bin_moffat_parameters[0] /= data_scaling_factor
    bin_moffat_parameters[4] /= data_scaling_factor
    bin_moffat_parameters[5] /= data_scaling_factor
    
    # Define the extraction limits of the current dispersion bin based
    # on the parameters of the Moffat profile previously fitted to it.
    limits_and_center = set_extraction_limits(
        bin_moffat_parameters,
        width_multiplier=fwhm_multiplier
    )
    bin_lower_extraction_limit = limits_and_center[0]
    bin_upper_extraction_limit = limits_and_center[1]
    moffat_center = limits_and_center[2]
    
    bin_dispersion_center = (each_bin[0] + each_bin[1]) * 0.5
    extraction_limits.append(
        [
            bin_dispersion_center,
            bin_lower_extraction_limit,
            bin_upper_extraction_limit,
            moffat_center
        ]
    )
    
    # Record the Moffat function parameters for each dispersion bin and
    # add the wavstart offset to the bin locations so they can be saved
    # as metadata along with the extracted spectrum.
    bin_moffat_parameters.append(each_bin[0] + axes_dict["wavelength_start"])
    bin_moffat_parameters.append(each_bin[1] + axes_dict["wavelength_start"])
    
    return bin_moffat_parameters, extraction_limits, bin_data


def moffat_least_squares(r, col, seeing, pixel_resolution):
    """
    Takes a data column, spatial axis and seeing of the observation and
    fits a Moffat function to the column using a least squares method. 
    Returns the best fit parameters of the Moffat function.

    Args:
     -- r (numpy.ndarray)
          Spatial axis of the data being fit.
     -- col (numpy.ndarray)
          The data being fitted.
     -- seeing (float)
          Estimated FWHM of the spatial profile.
     -- pixel_resolution (float)
          Spatial resolution of each pixel in arcsec/pixel.

    Returns:
     -- param_list (list)
          List of best fit output parameters returned by the least
          squares routine.
    """

    # Set up initial conditions for the least squares fit.
    # x0 = [amplitude, centre, alpha, beta, background level, background gradient]
    # Initial beta estimate comes from optimal value from atmospheric 
    # turbulence theory as described in Trujillo, I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    
    a0 = np.nanmedian(np.sort(col)[-3:])
    c0 = np.argmax(col)
    alpha0 = seeing / pixel_resolution
    beta0 = 4.765
    bg_level0 = np.median(np.concatenate((col[:5], col[-5:])))
    bg_grad0 = 0.

    x0 = [a0, c0, alpha0, beta0, bg_level0, bg_grad0]
    
    lower_bound = [0.0, np.argmax(col) - 1.0, 0.0, 0.0, -np.inf, -np.inf]
    upper_bound = [np.inf, c0+1, 5*alpha0, 5.0, np.inf, np.inf]

    # Run the least squares fit.
    res_lsq = least_squares(
        moffat_resid,
        x0,
        bounds=(lower_bound, upper_bound),
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
    """
    Calculates residuals of fitted moffat profile and the data for the
    least squares fitting.

    Description:
        A = x[0]
        c = x[1]
        alpha = x[2]
        beta = x[3]
        B = x[4]
        m = x[5]

    Args:
     -- x (numpy.ndarray)
          An array of parameters defining the shape of the model moffat
          profile.
     -- data_range (numpy.ndarray)
          The spatial axis of the data
     -- data (numpy.ndarray)
          The data

    Returns:
     -- residual (numpy.ndarray)
          The residual array between the model moffat profile and the
          data.
    """

    residual = common.moffat(x, data_range) - data
    return residual


def set_extraction_limits(moffat_parameters, width_multiplier=3.0):
    """
    Calculate the extraction limits from a Moffat profile based on the
    distance from the central peak as a multiple of FWHM.

    Args:
     -- moffat_parameters (list)
          List containing the parameters of the fitted Moffat profile. 
     -- width_multiplier (float, default=3.0)
          Defines the distance from the center of the spatial profile at
          which to set the extraction limits, in multiples of the FWHM.

    Returns:
     -- lower_extraction_limit (numpy.float64)
          The lower bound of the region to be extracted.
     -- upper_extraction_limit (numpy.float64)
          The upper bound of the region to be extracted.
     -- fwhm (numpy.float64)
          The Full Width at Half Maximum of the Moffat profile.
     -- moffat_parameters[1] (numpy.float64)
          Location of the center of the Moffat profile.
    """

    # Create a Moffat profile based on the input parameters.
    fwhm = common.moffat_fwhm(moffat_parameters[2], moffat_parameters[3])

    # Respectively define the upper and lower extraction limits at a
    # distance above and below the peak of the Moffat profile that
    # equals the FWHM of the Moffat profile times a multiplying factor.
    lower_extraction_limit = moffat_parameters[1] - (width_multiplier * fwhm)
    upper_extraction_limit = moffat_parameters[1] + (width_multiplier * fwhm)

    return [lower_extraction_limit, upper_extraction_limit, moffat_parameters[1]]
