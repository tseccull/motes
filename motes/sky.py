#!/home/tom/miniforge3/envs/work/bin/python

"""
sky.py - Contains all functions required for the sky subtraction
         capability of MOTES.
"""


import logging
import motes.diagnostics as diagnostics
import motes.tracing as tracing
import numpy as np

logger = logging.getLogger("motes")


def sky_median(good_sky_pixels):
    """
    Bootstrap the median sky level and the median's uncertainty from an
    array of sky pixels.
    
    Args:
     -- good_sky_pixels (numpy.ndarray)
          Numpy array/column of pixels sampling the sky background.
          
    Returns:
     -- column_sky_model (numpy.float32)
          Median sky level bootstrapped from the sky pixels.
     -- column_sky_model_err (numpy.float32)
          Uncertainty of the median sky level bootstrappted from the
          sky pixels.
    """
    
    bootstrapped_sky_pixels = np.random.choice(
        good_sky_pixels,
        (len(good_sky_pixels), 100), 
        replace=True
    )
    median_sky_sample = np.nanmedian(bootstrapped_sky_pixels, axis=0)
    column_sky_model = np.nanmean(median_sky_sample)
    column_sky_model_err = np.std(median_sky_sample) / (99**0.5)
    
    return column_sky_model, column_sky_model_err


def sky_locator(
    frame_dict,
    axes_dict,
    data_scaling_factor,
    header_parameters,
    bin_parameters,
    motes_parameters,
    input_file_name
):
    """
    Perform sky subtraction on the 2D spectrum. Localise the spectrum
    in the same way done for the extraction, and then use the regions
    outside the boundaries defined by that process to characterise and 
    subtract background sky emission.

    Args:
     -- frame_dict (dict)
          A dictionary containing the 2D spectrum and its associated
          errors and quality arrays.
     -- axes_dict (dict)
          A dictionary containing the wavelength and spatial axes of the
          2D spectrum.
     -- data_scaling_factor (float)
          A flux scale factor to convert the flux units of the 2D
          spectrum to the same units as the sky.
     -- header_parameters (dict)
          A dictionary containing the header parameters of the 2D 
          spectrum.
     -- bin_parameters (dict)
          A dictionary containing the bin parameters of the 2D spectrum.
     -- motes_parameters (dict)
          A dictionary containing the parameters of the extraction.
     -- input_file_name (str)
          Name of the current file being processed.

    Returns:
     -- frame_dict (dict)
          A dictionary containing the 2D spectrum and its associated
          errors and quality arrays.
     -- moffat_parameters_all_sky_bins (list)
          A list containing bins for the sky background.
     -- sky_extraction_limits (list)
          A list containing the extraction limits for the sky 
          background.
    """

    logger.info("Fitting Moffat Functions to each bin to localise sky.")

    moffat_parameters_all_sky_bins = []
    extraction_limits = []
    for each_bin in bin_parameters:
        bin_moffat_parameters, extraction_limits, bin_data = (
            tracing.moffat_fitting(
                each_bin,
                frame_dict["data"],
                axes_dict,
                data_scaling_factor,
                header_parameters,
                motes_parameters.sky_fwhm_multiplier,
                extraction_limits
            )
        )
        
        moffat_parameters_all_sky_bins.append(bin_moffat_parameters)

        # DIAGNOSTICS - Plot computed moffat profile over data for each bin.
        if motes_parameters.diag_plot_moffat:
            png_name = diagnostics.get_png_name(
                each_bin, axes_dict["wavelength_start"], input_file_name, "sky"
            )
            diagnostics.plot_fitted_spatial_profile(
                axes_dict["spatial_axis"],
                bin_data,
                axes_dict["hi_resolution_spatial_axis"],
                bin_moffat_parameters,
                axes_dict["data_spatial_floor"],
                header_parameters,
                png_name
            )

    moffat_parameters_all_sky_bins = np.array(moffat_parameters_all_sky_bins)
    logger.info("Fitting complete. Drawing target/sky boundaries.")
    extraction_limits = np.array(extraction_limits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D
    # spectrum.
    if motes_parameters.diag_plot_localisation:
        draw_lines = diagnostics.get_draw_lines(axes_dict, extraction_limits)
        title = "2D Spectrum Overplotted with Target/Sky Boundaries"
        data_frame = frame_dict["data"]
        diagnostics.show_img(data_frame, axes_dict, header_parameters, draw_lines, title)

    # Interpolate the extraction limits calculated for each median bin
    # such that each wavelength element across the entire unbinned
    # wavelength axis of the entire 2D spectrum has its own extraction
    # limits.
    sky_extraction_limits = tracing.interpolate_extraction_lims(
        extraction_limits, axes_dict["dispersion_axis_length"]
    )

    # DIAGNOSTICS - Plot the final extraction limits including the
    # extrapolated sections at the ends of the wavelength axis.
    if motes_parameters.diag_plot_localisation:
        draw_lines = diagnostics.get_draw_lines(axes_dict, sky_extraction_limits)
        title = "2D Spectrum Overplotted with Full Target/Sky Boundaries"
        data = frame_dict["data"]
        diagnostics.show_img(data, axes_dict, header_parameters, draw_lines, title)

    logger.info("Subtracting sky.")

    frame_dict = subtract_sky(
        sky_extraction_limits[0],
        sky_extraction_limits[1],
        frame_dict,
        axes_dict,
        motes_parameters,
        header_parameters,
        input_file_name
    )

    # DIAGNOSTICS - Plot the final extraction limits including the 
    # extrapolated sections at the ends of the wavelength axis.
    if motes_parameters.diag_plot_skysub:
        draw_lines = diagnostics.get_draw_lines(axes_dict, sky_extraction_limits)
        title = "Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries"
        data = frame_dict["data"]
        diagnostics.show_img(data, axes_dict, header_parameters, draw_lines, title)
    
    return frame_dict, moffat_parameters_all_sky_bins, sky_extraction_limits


def sky_polynomial(
    good_sky_pixels,
    good_sky_errs,
    good_sky_axis,
    column_axis,
    order,
    col,
    plot_fit,
    file_name,
    flux_unit
):
    """
    This function bootstraps a polynomial fit to sky pixels in a given
    pixel column to determine the best fit model and its uncertainties
    as input to sky subtraction.
    
    Args:
     -- good_sky_pixels (numpy.ndarray)
          Numpy array/column of pixels sampling the sky background.
     -- good_sky_errs (numpy.ndarray)
          Numpy array/column of pixels sampling the sky background
          uncertainties.
     -- good_sky_axis (numpy.ndarray)
          Spatial axis for good_sky_pixels.
     -- column_axis (numpy.ndarray)
          Full spatial axis for all column pixels including
          good_sky_pixels.
     -- order (int)
          Order of the polynomial fitted to the sky background.
     -- file_name (str)
          Name of the current file being processed.
     
    Return:
     -- column_sky_model (numpy.ndarray)
          Best fit sky model bootstrapped from the sky pixels via
          polynomial fitting.
     -- column_sky_model_err (numpy.ndarray)
          Uncertainty of the best fit sky model.
    """
    sample_size = 100
    bootstrapped_sky_pixels = []
    for j in range(len(good_sky_pixels)):
	    bootstapped_pixel = np.random.normal(
	        loc=good_sky_pixels[j],
	        scale=good_sky_errs[j], 
	        size=sample_size
	    )
	    bootstrapped_sky_pixels.append(bootstapped_pixel)
    bootstrapped_sky_pixels = np.array(bootstrapped_sky_pixels).T
    
    fit_coefficients = np.zeros((sample_size, order + 1))
    for j in range(sample_size):
        p_coeffs = np.polyfit(good_sky_axis, bootstrapped_sky_pixels[j], order)
        fit_coefficients[j] += np.flip(p_coeffs)
    coefficient_orders = np.arange(order + 1)
    coefficient_axes = [column_axis**x for x in coefficient_orders]
    coefficient_axes = np.array(coefficient_axes).T
    all_coefficient_axes = np.array([coefficient_axes]*sample_size)
    all_coefficients  = [[x]*len(column_axis) for x in fit_coefficients]
    all_coefficients = np.array(all_coefficients)
    sky_samples = np.sum((all_coefficient_axes * all_coefficients), axis=2).T
    column_sky_model = np.median(sky_samples, axis=1)
    column_sky_model_err = np.std(sky_samples ,axis=1)
    
    if plot_fit:
        diagnostics.plot_sky_fit(
            good_sky_axis,
            good_sky_pixels,
            good_sky_errs, 
            column_axis,
            column_sky_model,
            column_sky_model_err,
            col,
            file_name,
            flux_unit
        )
    
    return column_sky_model, column_sky_model_err


def subtract_sky(
    background_spatial_lo_limit,
    background_spatial_hi_limit,
    frame_dict,
    axes_dict,
    parameters,
    header_parameters,
    input_file_name
):
    """
    Subtracts the sky background from the 2D image by defining bg
    regions using limits input to the function and then fitting a
    profile to the background column by column while masking cosmic
    rays. The background level or profile for each column is subtracted
    from the full column to produce a background subtracted 2D image.

    Args:
     -- background_spatial_lo_limit (numpy.ndarray)
          Lower limits of the background regions.
     -- background_spatial_hi_limit (numpy.ndarray)
          Upper limits of the background regions.
     -- frame_dict (dict)
          A dictionary containing the 2D image.
     -- axes_dict (dict)
          A dictionary containing the axis information.
     -- parameters (dict)
          A dictionary containing MOTES parameters.
     -- header_parameters (dict)
          A dictionary containing header information.
     -- input_file_name (str)
          Name of the current file being processed.

    Returns:
     -- frame_dict (dict)
          Dictionary containing the background subtracted 2D image.
    """

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    sky_model = []
    sky_model_err = []
    number_of_columns = len(background_spatial_lo_limit)

    # Makes sure the limits are within the image.
    for ii in range(number_of_columns):
        if background_spatial_lo_limit[ii] < 0:
            background_spatial_lo_limit[ii] = 0
        data_section_upper_limit = axes_dict["spatial_axis"] + axes_dict["data_spatial_floor"]
        if background_spatial_hi_limit[ii] > data_section_upper_limit[-1]:
            background_spatial_hi_limit[ii] = data_section_upper_limit[-1]

        data_column = frame_dict["data"][ii]
        errs_column = frame_dict["errs"][ii]
        
        column_axis = np.array(range(len(data_column))) + axes_dict["data_spatial_floor"]
        
        in_sky_region = np.logical_or(
            column_axis < background_spatial_lo_limit[ii] + axes_dict["data_spatial_floor"],
            column_axis > background_spatial_hi_limit[ii] + axes_dict["data_spatial_floor"]
        )
        sky_pixels = data_column[np.where(in_sky_region)]
        sky_errors = errs_column[np.where(in_sky_region)]

        # Kill MOTES if there isn't enough background sky to perform sky
        # subtraction. Should probably be made more nuanced later.
        if len(sky_pixels) == 0:
            sky_fwhm_mult_limit = (
                round(np.min(np.shape(frame_dict["data"])) / (2 * header_parameters["seeing"]), 1)
            )
            logger.critical("No pixels contained inside sky region.")
            logger.critical("SKY_FWHM_MULTIPLIER is probably too large.")
            logger.critical(
                "SKY_FWHM_MULTIPLIER < %s is recommended in this case.", 
                sky_fwhm_mult_limit
            )
            logger.critical("Enlarging the 2D spectrum region in reg.txt is also a viable solution")
            logger.critical("Raising RuntimeError")
            raise RuntimeError

        sky_axis = column_axis[in_sky_region]
        if len(set(sky_pixels)) == 1:
            continue

        sky_pixels_median = np.nanmedian(sky_pixels)
        sky_pixels_sigma = np.nanstd(sky_pixels)
        good_sky_pixels_bool = np.logical_and(
            sky_pixels > sky_pixels_median - (10 * sky_pixels_sigma),
            sky_pixels < sky_pixels_median + (10 * sky_pixels_sigma),
        )
        good_sky_pixels_location = np.where(good_sky_pixels_bool)
        good_sky_pixels = sky_pixels[good_sky_pixels_location]
        good_sky_axis = sky_axis[good_sky_pixels_location]
        good_sky_errs = sky_errors[good_sky_pixels_location]
        order = int(parameters.sky_order)

        if order == 0:
            column_sky_model, column_sky_model_err = sky_median(good_sky_pixels)
            sky_model.append(column_sky_model)
            sky_model_err.append(column_sky_model_err)
        
        if order > 0:
            column_sky_model, column_sky_model_err = (
                sky_polynomial(
                    good_sky_pixels,
                    good_sky_errs,
                    good_sky_axis,
                    column_axis,
                    order,
                    ii + axes_dict["wavelength_start"],
                    parameters.diag_plot_skyfit,
                    input_file_name,
                    header_parameters["flux_unit"]
                )
            )
            sky_model.append(column_sky_model)
            sky_model_err.append(column_sky_model_err)

        frame_dict["data"][ii] -= column_sky_model
        frame_dict["errs"][ii] = (
            ((frame_dict["errs"][ii] ** 2) + (column_sky_model_err**2)) ** 0.5
        )

    sky_model = np.array(sky_model)
    sky_model_err = np.array(sky_model_err)
    if parameters.sky_order == 0:
        frame_dict["sky_model"] = (
            np.tile(sky_model, (np.shape(frame_dict["data"])[1], 1))
        )
        frame_dict["sky_uncertainty"] = (
            np.tile(sky_model_err, (np.shape(frame_dict["data"])[1], 1))
        )
    else:
        frame_dict["sky_model"] = sky_model.T
        frame_dict["sky_uncertainty"] = sky_model_err.T

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    chip_gaps = np.where(np.nanmedian(frame_dict["data"], axis=0) == 0)
    frame_dict["data"][:, chip_gaps[0]] += 1

    return frame_dict
