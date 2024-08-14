#!/home/tom/miniforge3/envs/work/bin/python

"""
sky.py - Contains all functions required for the sky subtraction
         capability of MOTES.
"""


import sys
import motes.common as common
import numpy as np

from scipy.optimize import least_squares


# Takes a data column, spatial axis and seeing of the observation and
# fits a Moffat function to the column using a Levenberg-Marquardt least
# squares method. Returns the best fit parameters of the Moffat function
def linear_least_squares(r, col):
    """
    Fit a linear profile to a data column using a Levenberg-Marquardt
    least squares method.

    Args:
     -- r (numpy.ndarray)
          Spatial axis of the data column being fitted.
     -- col (numpy.ndarray)
          Data column being fitted.

    Returns:
     -- [res_lsq.x[0], res_lsq.x[1]] (list)
          list containing the best fit parameters of the linear profile.
    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(linear_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1]]


def linear_resid(x, data_range, data):
    """
    Calculate residuals of fitted linear profile and the data for the
    Levenberg-Marquardt least squares method.

    Args:
     -- x (list)
          A list containing the best fit parameters of the linear
          profile.
     -- data_range (numpy.ndarray)
          the spatial axis of the data column.
     -- data (numpy.ndarray)
          the data column.

    Returns:
     -- residual (numpy.ndarray)
          The residuals of the fitted line and the data.
    """
    residual = (x[0] * data_range) + x[1] - data
    return residual


def poly2_least_squares(r, col):
    """
    Fits a second-order polynomial to a data column using a
    Levenberg-Marquardt least squares method.

    Args:
     -- r (numpy.ndarray)
          Spatial axis of the data column being fitted.
     -- col (numpy.ndarray)
          Data column being fitted.

    Returns:
     -- params (list)
          the parameters of the fitted polynomial.
    """

    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly2_resid, x0, args=(r, col), method="trf")

    params = [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2]]
    return params


def poly2_resid(x, data_range, data):
    """
    Calculates residuals of a fitted second-order polynomial and the
    data for the Levenberg-Marquardt least squares method.

    Args:
     -- x (list)
          A list containing the best fit parameters of the polynomial.
     -- data_range (numpy.ndarray)
          the spatial axis of the data column.
     -- data (numpy.ndarray)
          the data column.

    Returns:
     -- residual (numpy.ndarray)
          The residuals of the fitted polynominal and the data.
    """
    residual = (
        (x[0] * data_range * data_range) + (x[1] * data_range) + x[2] - data
    )

    return residual


def poly3_least_squares(r, col):
    """
    Fits a third-order polynomial to a data column using a 
    Levenberg-Marquardt least squares method.

    Args:
     -- r (numpy.ndarray)
          Spatial axis of the data column being fitted.
     -- col (numpy.ndarray)
          Data column being fitted.

    Returns:
     -- params (list)
          the best-fit parameters of the fitted polynomial.
    """
    # Set up initial conditions for the least squares fit.
    x0 = [0.0, 0.0, 0.0, np.median(col)]

    # Run the least squares fit.
    res_lsq = least_squares(poly3_resid, x0, args=(r, col), method="trf")

    return [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2], res_lsq.x[3]]


def poly3_resid(x, data_range, data):
    """
    Calculates residuals of a fitted third-order polynomial and the data
    for the Levenberg-Marquardt least squares method.

    Args:
     -- x (list)
          A list containing the best fit parameters of the polynomial.
     -- data_range (numpy.ndarray)
          the spatial axis of the data column.
     -- data (numpy.ndarray)
          the data column.

    Returns:
     -- residual (numpy.ndarray)
          The residuals of the fitted polynominal and the data.
    """
    residual = (
        (x[0] * data_range * data_range * data_range)
        + (x[1] * data_range * data_range)
        + (x[2] * data_range)
        + x[3]
        - data
    )
    return residual


def sky_locator(
    frame_dict,
    axes_dict,
    data_scaling_factor,
    header_parameters,
    bin_parameters,
    motes_parameters
):
    """
    Perform sky subtraction on the 2D spectrum. Locaalise the spectrum
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

    sys.stdout.write(
        " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
    )

    moffat_parameters_all_sky_bins = []
    extraction_limits = []
    for each_bin in bin_parameters:
        # Take the median spatial profile of the dispersion bin, and
        # leave out pixel columns in the chip gaps if this is a GMOS
        # spectrum.
        raw_bin_data = frame_dict["data"][:, each_bin[0] : each_bin[1]]
        chip_gap_location = np.where(np.median(raw_bin_data, axis=0) != 1)
        bin_data = np.nanmedian(raw_bin_data[:, chip_gap_location[0]], axis=1)

        # Use a Levenberg-Marquardt Least Squares method to fit a Moffat
        # function to the median spatial profile and return its
        # parameters.
        bin_moffat_parameters = common.moffat_least_squares(
            axes_dict["spatial_axis"],
            bin_data * data_scaling_factor,
            header_parameters["seeing"],
            header_parameters["pixel_resolution"],
        )

        # Scipy least squares doesn't like really tiny numbers like
        # fluxes in erg/s/cm^2/Angstrom, so it's necessary to scale the
        # data to a size that least squares can handle. The shape of the
        # profile fitted to the scaled spatial profile is the same as
        # the unscaled, but to to get a model profile that matches the
        # original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to
        # be scaled down again after the fitting.
        bin_moffat_parameters[0] /= data_scaling_factor
        bin_moffat_parameters[1] += axes_dict["data_spatial_floor"]
        bin_moffat_parameters[4] /= data_scaling_factor
        bin_moffat_parameters[5] /= data_scaling_factor

        # Define the extraction limits of the current dispersion bin
        # based on the parameters of the Moffat profile previously
        # fitted to it.
        limits_and_center = common.set_extraction_limits(
            bin_moffat_parameters,
            width_multiplier=motes_parameters["-BG_FWHM_MULTIPLIER"],
        )
        
        bin_lower_extraction_limit = limits_and_center[0]
        bin_upper_extraction_limit = limits_and_center[1]
        moffat_center = limits_and_center[2]

        extraction_limits.append(
            [
                (each_bin[0] + each_bin[1]) * 0.5,
                bin_lower_extraction_limit,
                bin_upper_extraction_limit,
                moffat_center
            ]
        )

        # Record the Moffat function parameters for each dispersion bin
        # and add the wavstart offset to the bin locations so they can
        # be saved as metadata along with the extracted spectrum.
        bin_moffat_parameters.append(
            each_bin[0] + axes_dict["wavelength_start"]
        )
        bin_moffat_parameters.append(
            each_bin[1] + axes_dict["wavelength_start"]
        )
        moffat_parameters_all_sky_bins.append(bin_moffat_parameters)

        # DIAGNOSTICS - Plot computed moffat profile over data for each
        # bin
        if motes_parameters["-DIAG_PLOT_MOFFAT"]:
            common.plot_fitted_spatial_profile(
                axes_dict["spatial_axis"],
                bin_data,
                axes_dict["hi_resolution_spatial_axis"],
                bin_moffat_parameters,
                axes_dict["data_spatial_floor"],
                header_parameters,
            )

    moffat_parameters_all_sky_bins = (
        np.array(moffat_parameters_all_sky_bins)
    )

    sys.stdout.write("     Fitting complete.\n")

    sys.stdout.write(" >>> Drawing target/sky boundaries. ")
    sys.stdout.flush()
    extraction_limits = np.array(extraction_limits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D
    # spectrum.
    if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        draw_lines = [
            extraction_limits[0] + axes_dict["wavelength_start"],
            (extraction_limits[1]) + axes_dict["data_spatial_floor"],
            extraction_limits[0] + axes_dict["wavelength_start"],
            (extraction_limits[2]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            draw_lines,
            "2D Spectrum Overplotted with Target/Sky Boundaries",
        )

    # Interpolate the extraction limits calculated for each median bin
    # such that each wavelength element across the entire unbinned
    # wavelength axis of the entire 2D spectrum has its own extraction
    # limits.
    sky_extraction_limits = common.interpolate_extraction_lims(
        extraction_limits, axes_dict["dispersion_axis_length"]
    )
    sys.stdout.write("DONE.\n")

    # DIAGNOSTICS - Plot the final extraction limits including the
    # extrapolated sections at the ends of the wavelength axis.
    if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        draw_lines = [
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
            (sky_extraction_limits[0]) + axes_dict["data_spatial_floor"],
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
            (sky_extraction_limits[1]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            draw_lines,
            "2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )

    sys.stdout.write(" >>> Subtracting sky.\n")
    sys.stdout.flush()

    frame_dict = subtract_sky(
        sky_extraction_limits[0],
        sky_extraction_limits[1],
        frame_dict,
        axes_dict,
        motes_parameters,
        header_parameters,
    )

    sys.stdout.write("\n     DONE.\n")
    sys.stdout.flush()

    # DIAGNOSTICS - Plot the final extraction limits including the 
    # extrapolated sections at the ends of the wavelength axis.
    if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        draw_lines = [
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
            (sky_extraction_limits[0]) + axes_dict["data_spatial_floor"],
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
            (sky_extraction_limits[1]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            draw_lines,
            "Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )
    
    return frame_dict, moffat_parameters_all_sky_bins, sky_extraction_limits


def subtract_sky(
    background_spatial_lo_limit,
    background_spatial_hi_limit,
    frame_dict,
    axes_dict,
    parameters,
    header_parameters
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

    Returns:
     -- frame_dict (dict)
          Dictionary containing the background subtracted 2D image.
    """

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    sky_model = []
    number_of_columns = len(background_spatial_lo_limit)

    # Makes sure the limits are within the image.
    for ii in range(number_of_columns):
        if background_spatial_lo_limit[ii] < 0:
            background_spatial_lo_limit[ii] = 0
        if (
            background_spatial_hi_limit[ii] > 
            (axes_dict["spatial_axis"] + axes_dict["data_spatial_floor"])[-1]
        ):
            background_spatial_hi_limit[ii] = (
                (axes_dict["spatial_axis"]
                + axes_dict["data_spatial_floor"])[-1]
            )

        data_column = frame_dict["data"][ii]
        column_axis = np.array(range(len(data_column)))
        sky_pixels = data_column[
            np.where(
                np.logical_or(
                    column_axis < background_spatial_lo_limit[ii],
                    column_axis > background_spatial_hi_limit[ii]
                )
            )
        ]

        # Kill MOTES if there isn't enough background sky to perform sky
        # subtraction. Should probably be made more nuanced later.
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
            np.where(
                np.logical_or(
                    column_axis < background_spatial_lo_limit[ii],
                    column_axis > background_spatial_hi_limit[ii]
                )
            )
        ]
        if len(set(sky_pixels)) == 1:
            continue

        sky_pixels_median = np.nanmedian(sky_pixels)
        sky_pixels_sigma = np.nanstd(sky_pixels)
        good_sky_pixels_location = np.where(
            np.logical_and(
                sky_pixels > sky_pixels_median - (10 * sky_pixels_sigma),
                sky_pixels < sky_pixels_median + (10 * sky_pixels_sigma),
            )
        )
        good_sky_pixels = sky_pixels[good_sky_pixels_location]
        good_sky_axis = sky_axis[good_sky_pixels_location]

        if parameters["-SKYSUB_MODE"] == "MEDIAN":
            bootstrapped_sky_pixels = (
                np.random.choice(
                    good_sky_pixels, (len(good_sky_pixels), 100), replace=True
                )
            )
            median_sky_sample = np.nanmedian(bootstrapped_sky_pixels, axis=0)
            column_sky_model = np.nanmean(median_sky_sample)
            sky_model.append(column_sky_model)
            sky_model_err = np.std(median_sky_sample) / (99**0.5)

        if parameters["-SKYSUB_MODE"] == "LINEAR":
            bootstrapped_sky_pixels = (
                np.random.choice(
                    good_sky_pixels, (len(good_sky_pixels), 100), replace=True
                )
            )
            linear_terms = []
            intercepts = []
            for jj in bootstrapped_sky_pixels.T:
                sky_model_parameters = linear_least_squares(good_sky_axis, jj)
                linear_terms.append(sky_model_parameters[0])
                intercepts.append(sky_model_parameters[1])
            intercepts = np.array(intercepts)
            linear_terms = np.array(linear_terms)
            sky_linear_term = np.mean(linear_terms)
            sky_linear_term_err = np.std(linear_terms) / (99**0.5)
            sky_intercept = np.mean(intercepts)
            sky_intercept_err = np.std(intercepts) / (99**0.5)
            column_sky_model = (sky_linear_term * column_axis) + sky_intercept
            sky_model.append(column_sky_model)
            sky_model_err = (
                (
                    sky_linear_term_err * sky_linear_term_err
                    * column_axis * column_axis
                )
                + (sky_intercept_err * sky_intercept_err)
            ) ** 0.5

        if parameters["-SKYSUB_MODE"] == "POLY2":
            bootstrapped_sky_pixels = (
                np.random.choice(
                    good_sky_pixels, (len(good_sky_pixels), 100), replace=True
                )
            )
            linear_terms = []
            intercepts = []
            quadratic_terms = []
            for jj in bootstrapped_sky_pixels.T:
                sky_model_parameters = poly2_least_squares(good_sky_axis, jj)
                quadratic_terms.append(sky_model_parameters[0])
                linear_terms.append(sky_model_parameters[1])
                intercepts.append(sky_model_parameters[2])
            quadratic_terms = np.array(quadratic_terms)
            intercepts = np.array(intercepts)
            linear_terms = np.array(linear_terms)
            sky_quadratic_term = np.mean(quadratic_terms)
            sky_quadratic_term_err = np.std(quadratic_terms) / (99**0.5)
            sky_linear_term = np.mean(linear_terms)
            sky_linear_term_err = np.std(linear_terms) / (99**0.5)
            sky_intercept = np.mean(intercepts)
            sky_intercept_err = np.std(intercepts) / (99**0.5)
            column_sky_model = (
                (sky_quadratic_term * column_axis * column_axis)
                + (sky_linear_term * column_axis)
                + sky_intercept
            )
            sky_model.append(column_sky_model)
            sky_model_err = (
                (
                    sky_quadratic_term_err * sky_quadratic_term_err 
                    * column_axis * column_axis * column_axis * column_axis
                )
                + (
                    sky_linear_term_err * sky_linear_term_err 
                    * column_axis * column_axis
                )
                + (sky_intercept_err * sky_intercept_err)
            ) ** 0.5

        if parameters["-SKYSUB_MODE"] == "POLY3":
            bootstrapped_sky_pixels = (
                np.random.choice(
                    good_sky_pixels, (len(good_sky_pixels), 100), replace=True
                )
            )
            linear_terms = []
            intercepts = []
            quadratic_terms = []
            cubic_terms = []
            for jj in bootstrapped_sky_pixels.T:
                sky_model_parameters = poly3_least_squares(good_sky_axis, jj)
                cubic_terms.append(sky_model_parameters[0])
                quadratic_terms.append(sky_model_parameters[1])
                linear_terms.append(sky_model_parameters[2])
                intercepts.append(sky_model_parameters[3])
            cubic_terms = np.array(cubic_terms)
            quadratic_terms = np.array(quadratic_terms)
            intercepts = np.array(intercepts)
            linear_terms = np.array(linear_terms)
            sky_cubic_term = np.mean(cubic_terms)
            sky_cubic_term_err = np.std(cubic_terms) / (99**0.5)
            sky_quadratic_term = np.mean(quadratic_terms)
            sky_quadratic_term_err = np.std(quadratic_terms) / (99**0.5)
            sky_linear_term = np.mean(linear_terms)
            sky_linear_term_err = np.std(linear_terms) / (99**0.5)
            sky_intercept = np.mean(intercepts)
            sky_intercept_err = np.std(intercepts) / (99**0.5)
            column_sky_model = (
                (sky_cubic_term * column_axis * column_axis * column_axis)
                + (sky_quadratic_term * column_axis * column_axis)
                + (sky_linear_term * column_axis)
                + sky_intercept
            )
            sky_model.append(column_sky_model)
            sky_model_err = (
                (
                    sky_cubic_term_err
                    * sky_cubic_term_err
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                    * column_axis
                )
                + (
                    sky_quadratic_term_err * sky_quadratic_term_err
                    * column_axis * column_axis * column_axis * column_axis
                )
                + (
                    sky_linear_term_err * sky_linear_term_err 
                    * column_axis * column_axis
                )
                + (sky_intercept_err * sky_intercept_err)
            ) ** 0.5

        frame_dict["data"][ii] -= column_sky_model
        frame_dict["errs"][ii] = (
            ((frame_dict["errs"][ii] ** 2) + (sky_model_err**2)) ** 0.5
        )

        sys.stdout.write(
            "     " 
            + str(ii + 1)
            + "/"
            + str(number_of_columns)
            + " columns completed.\r"
        )

    sky_model = np.array(sky_model)
    if parameters["-SKYSUB_MODE"] == "MEDIAN":
        frame_dict["sky_model"] = (
            np.tile(sky_model, (np.shape(frame_dict["data"])[1], 1))
        )
    else:
        frame_dict["sky_model"] = sky_model.T

    frame_dict["data"] = frame_dict["data"].T
    frame_dict["errs"] = frame_dict["errs"].T

    chip_gaps = np.where(np.nanmedian(frame_dict["data"], axis=0) == 0)
    frame_dict["data"][:, chip_gaps[0]] += 1

    return frame_dict
