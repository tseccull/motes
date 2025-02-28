#!/usr/bin/env python3

import logging
import motes.common as common
import numpy as np

logger = logging.getLogger("motes")


def filter_data(data_2d, errs_2d):
    """
    This function takes in the data_2d and errs_2d and outputs frames
    where any NaN or Inf values are 0.0 (this can be later filtered to a
    median of the column). This is used before the extraction procedures
    to ensure the S/N in the bin is numerical (required for optimal
    extraction).

    Args:
     -- data_2d (numpy.ndarray)
          Original data_2d
     -- errs_2d (numpy.ndarray)
          Original errs_2d

    Returns:
     -- data_2d (numpy.ndarray)
          Filtered data_2d
     -- errs_2d (numpy.ndarray)
          Filtered errs_2d
    """

    errs_2d[~np.isfinite(errs_2d)] = 0.0
    data_2d[~np.isfinite(errs_2d)] = 0.0
    data_2d[~np.isfinite(data_2d)] = 0.0
    errs_2d[~np.isfinite(data_2d)] = 0.0

    return data_2d, errs_2d


def optimal_extraction(
    data_2d, errs_2d, extraction_limits, bin_parameters, axes_dict
):
    """
    Perform optimal extraction using a modified version of Horne (1986)
    where S=0, G=0 and errors are not 'revised' since we already have
    the 'variance frame'. Ideally, this extraction reduces errors with
    the weighting and conserves the flux (when compared to the standard
    extraction). This subroutine uses an analytic Moffat profile (post
    sky-subtraction), instead of polynominals as in Horne (1986).
    Furthermore, profiling takes place bin by bin, accounting for
    spatial profile variations across the dispersion axis; extraction
    takes place column by column within the bin limits and within the
    extraction limits previously calculated.

    Args:
     -- data_2d (numpy.ndarray)
          Input data frame
     -- errs_2d (numpy.ndarray)
          Input error frame
     -- extraction_limits (numpy.ndarray)
          An array containing limits at each dispersion pixel
     -- bin_parameters (list)
          A list containing the bin limits across the dispersion axis,
          to enable slicing the data across dispersion axis.
     -- axes_dict (dict)
          A dictionary containing the spatial axis array and other
          relevant information about the size and shape of the data
          frame.

    Returns:
     -- optimal_1d_data (numpy.ndarray)
          1D array of the optimally extracted spectrum.
     -- optimal_1d_errs (numpy.ndarray)
          1D array of the uncertainties of the optimally extracted
          spectrum.
     -- aperture_1d_data (numpy.ndarray)
          1D array of the aperture extracted spectrum.
     -- aperture_1d_errs (numpy.ndarray)
          1D array of the uncertainties of the aperture extracted
          spectrum.
    """

    logger.info("Performing optimal extraction of 1D spectrum.")

    # Filter any NaNs and Inf for data/errs AND ensure the errors are
    # positive for this extraction.
    data_2d, errs_2d = filter_data(data_2d, np.abs(errs_2d))

    # Set up output arrays for the optimally and aperture extracted
    # spectra and their respective uncertainties
    optimal_1d_data = np.zeros(np.shape(data_2d)[0])
    optimal_1d_errs = np.zeros(np.shape(data_2d)[0])
    aperture_1d_data = np.zeros(np.shape(data_2d)[0])
    aperture_1d_errs = np.zeros(np.shape(data_2d)[0])

    # Set up bin identification parameters for the following loop.
    bin_number = -1
    b = np.zeros(8)
    end_bin_number = len(bin_parameters) - 1

    # Loop through each dispersion element of the spectrum.
    for i, col in enumerate(data_2d):
        # Identify the location of the current element in the original
        # 2D spectrum.
        first_column = i + axes_dict["wavelength_start"]
        next_bin_number = bin_number + 1
        
        # If the current element belongs in the next bin as defined by
        # getbins, use the new bin's parameters and increment the bin 
        # number.
        if (bin_number < end_bin_number and first_column == bin_parameters[next_bin_number][-2]):
            bin_number += 1
            b = bin_parameters[bin_number]

        # Get the extraction limits for the current dispersion element
        # and define the spatial axis. Where the extraction limits
        # include partial pixels on the edge of the aperture, those
        # pixels are included in their entirety.
        fine_lo_extraction_limit = extraction_limits[0][i]
        fine_hi_extraction_limit = extraction_limits[1][i]
        lo_extraction_limit = int(np.floor(extraction_limits[0][i]))
        hi_extraction_limit = int(np.ceil(extraction_limits[1][i])) + 1
        ax = axes_dict["spatial_axis"][lo_extraction_limit : hi_extraction_limit]

        # Use the extraction limits to define the data column for this
        # wavelength element.
        col = col[lo_extraction_limit : hi_extraction_limit]

        # Use the extraction limits to define the errs column for this
        # wavelength element. Where errs have a value of 0, set them to
        # the median err value for the entire column. If all errs have a
        # value of zero, the column is likely to be bad and should have
        # the optimal extraction outputs set to zero.
        err = errs_2d[i]
        if all(x == 0.0 for x in err):
            optimal_1d_data[i] += 0.0
            optimal_1d_errs[i] += 0.0
            continue

        err[np.where(err == 0)] = np.median(err[np.where(err != 0)])
        err = err[lo_extraction_limit : hi_extraction_limit]

        # Use the err column to get the variance column
        var = err * err

        # Perform the standard aperture extraction on the data in the
        # column, and add this value to the aperture 1D output array.
        # Get the root sum square of the err column and add this value
        # to the uncertainty array for the 1D aperture extraction.
        f = np.sum(col)
        aperture_1d_data[i] += f
        aperture_1d_errs[i] += np.sqrt(np.sum(var))

        # Step 5 of the Horne 1986 algorithm - Defining the spatial
        # profile. Use the average value of the extraction limits in the
        # current column to estimate the location of the centre of the 
        # peak of the spatial profile.
        profile_center = (fine_hi_extraction_limit + fine_lo_extraction_limit) * 0.5

        # Use the Moffat profile parameters for the current bin to make
        # a moffat profile that approximates the shape of the spectrum's
        # spatial profile in the current column. The estimated centre of
        # the extraction limits is used to shift the profile to the
        # correct location along the spatial axis. The background is
        # assumed to have been subtracted by now, so the background
        # level and gradient are set to 0. Because the profile is a
        # model PSF, none of its values are negative and positivity of
        # the profile need not be enforced.
        profile = common.moffat([b[0], profile_center, b[2], b[3], 0.0, 0.0], ax)
        # Normalize the profile such that its sum equals unity.
        normalized_profile = profile / np.sum(profile)

        # Step 7 of Horne 1986 - Cosmic ray masking with a conservative
        # 5 sigmas.
        value = np.abs((col - (f * normalized_profile)) ** 2)
        sigma = 5
        condition = var * sigma**2
        # With invert, it gives False (ie. 0) if the above condition is
        # satisfied which means to mask the rays
        cosmic_mask = np.invert([value > condition])[0]

        # Add an exception to handle scenarios where the mask is all
        # False If that happens, ignore the cosmic ray correction.
        if not np.any(cosmic_mask):
            cosmic_mask = np.full(shape=cosmic_mask.shape, fill_value=True)

        # Extract the optimal spectrum (Step 8 of the algorithm) Page 4
        # of Horne 1986: These formulae are equivalent to determining
        # the OPT spectrum by scaling a known spatial profile P, to fit
        # the sky subtracted data, whose variances are V.

        # If column is in a GMOS chip gap set the optimal flux and
        # uncertainty to 0.
        if all(x == 0.0 for x in col) or all(x == 1.0 for x in col):
            fopt = 0.0
            vopt = 0.0
        else:
            fopt_numerator = np.sum(((normalized_profile * col) / var)[cosmic_mask])
            fopt_denominator = np.sum(((normalized_profile ** 2) / var)[cosmic_mask])
            fopt = fopt_numerator / fopt_denominator
            
            vopt_numerator = np.sum(normalized_profile[cosmic_mask])
            vopt_denominator = np.sum(((normalized_profile ** 2) / var)[cosmic_mask])
            vopt = vopt_numerator / vopt_denominator

        optimal_1d_data[i] += fopt
        optimal_1d_errs[i] += np.sqrt(vopt)
        
    extracted_data = [optimal_1d_data, optimal_1d_errs, aperture_1d_data, aperture_1d_errs]

    return extracted_data


