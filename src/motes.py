"""
MOTES Modular and Optimal Tracer and Extractor of Spectra.
Description: Modular and Optimal Tracer and Extractor of Specrtra (MOTES). A Python package for
extracting spectrum from astronomical 2D spectrograms.
Version: 0.4.5
Date: 2023-12-23
Authors: Tom Seccull, Dominik Kiersz
Licence: GNU General Public License v3.0
"""

import copy
import datetime
import glob
import sys

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

import motes.common as common
import motes.harvester as harvester
import motes.startup as startup


def motes():
    """
    Run MOTES. This function is the entrypoint for the MOTES. It loads .fits files from the input
    subdirectory and harvests the 2D data frames. It estimates bins for the spectrum. It then
    performs optional sky subtraction and optimal extraction of 1D spectra.
    """
    # Run startup functions
    motes_parameters = startup.read_motes_parameter_file()  # Import parameters from file to dict
    data_regions = startup.read_regions()  # Search for, and read in, reg.txt

    # Open and process each spectrum contained in the current directory.
    for i, input_file_path in enumerate(sorted(glob.glob("./inputs/*.fits"))):
        sys.stdout.write(("/" * (70 - len(input_file_path[:70]))) + " " + input_file_path[:70] + "\n")
        sys.stdout.write(" >>> Beginning MOTES Processing\n")

        # Gather header metadata and the image data from the 2D image file.
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file.\n"
        )
        header_parameters, frame_dict, axes_dict, input_file_primary_header = harvester.data_harvest(
            i, input_file_path, data_regions
        )
        # Make backup copies of the original data and error frames.
        frame_dict["original_data"] = copy.deepcopy(frame_dict["data"])
        frame_dict["original_errs"] = copy.deepcopy(frame_dict["errs"])
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file completed.\n"
        )

        # Perform initial least-squares Moffat fitting over the entire 2D spectrum collapsed along
        # the dispersion axis with a median.
        sys.stdout.write(
            " >>> Fitting Moffat profile to median spatial profile of entire spectrum. "
        )
        sys.stdout.flush()

        # Calculate median spatial profile of the spectrum.
        full_median_spatial_profile = np.nanmedian(frame_dict["data"], axis=1)
        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle.
        # The shape of the profile fitted to the scaled spatial profile is the same as the
        # unscaled, but to to get a model profile that matches the original profile, the profile
        # amplitude (A), background level (B), and background gradient (m) all need to be scaled
        # down again after the fitting.
        data_scaling_factor = 10 ** np.abs(
            np.floor(np.log10(np.abs(np.nanmedian(full_median_spatial_profile))))
        )
        # Fit the median spatial profile with a Moffat function.
        moffat_profile_parameters = common.moffat_least_squares(
            axes_dict["spatial_axis"],
            full_median_spatial_profile * data_scaling_factor,
            header_parameters["seeing"],
            header_parameters["pixel_resolution"],
        )
        # Get an improved estimate of the FWHM of the spectrum from the best fit Moffat profile.
        header_parameters["seeing"] = (
            2 * moffat_profile_parameters[2] * np.sqrt((2 ** (1 / moffat_profile_parameters[3])) - 1)
        )
        # Scale the amplitude, background gradient, and background level of the model Moffat
        # profile down.
        moffat_profile_parameters[0] /= data_scaling_factor
        moffat_profile_parameters[4] /= data_scaling_factor
        moffat_profile_parameters[5] /= data_scaling_factor

        sys.stdout.write("DONE.\n")
        sys.stdout.write(
            " >>> FWHM of median spatial profile is "
            + str(round(header_parameters["seeing"], 2))
            + " Pixels, or "
            + str(round(header_parameters["seeing"] * header_parameters["pixel_resolution"], 2))
            + '"\n'
        )

        # Use the parameters of the Moffat profile fitted to the median spatial profile of the
        # entire spectrum to determine spatial limits that are used to bound the region of the
        # spectrum used by the common.get_bins function to bin the 2D spectrum while taking account
        # of its S/N.
        binning_region_spatial_floor, binning_region_spatial_ceiling, moffat_fwhm, moffat_center = common.set_extraction_limits()moffat_profile_parameters)
        sys.stdout.write(
            " >>> Spectrum localised to aperture in range of spatial pixel rows "
            + str(int(binning_region_spatial_floor + axes_dict["data_spatial_floor"]))
            + "-"
            + str(int(binning_region_spatial_ceiling + axes_dict["data_spatial_floor"]))
            + "\n"
        )

        # DIAGNOSTICS -  Plot fitted Moffat profile over collapsed 2D spectrum and print the
        # parameters of the fitted Moffat profile.
        if motes_parameters["-DIAG_PLOT_COLLAPSED_2D_SPEC"]:
            common.print_moffat_parameters(moffat_profile_parameters, axes_dict["data_spatial_floor"], data_scaling_factor)
            common.plot_fitted_spatial_profile(
                axes_dict["spatial_axis"],
                full_median_spatial_profile,
                axes_dict["hi_resolution_spatial_axis"],
                moffat_profile_parameters,
                axes_dict["data_spatial_floor"],
                header_parameters,
            )

        # Determine the location of bins on the dispersion axis within which to measure the spatial
        # profile.
        bin_parameters, frame_dict = common.get_bins(
            frame_dict,
            int(np.floor(binning_region_spatial_floor)),
            int(np.ceil(binning_region_spatial_ceiling)),
            axes_dict["dispersion_axis_length"],
            motes_parameters,
            has_sky=True,
        )

        # Will plot the location of the bins determined by get_bins if -DIAG_PLOT_BIN_LOC=1 in
        # motesparams.txt
        common.get_bins_output(
            bin_parameters,
            motes_parameters,
            binning_region_spatial_floor,
            binning_region_spatial_ceiling,
            frame_dict["data"],
            header_parameters,
            axes_dict,
        )
        sys.stdout.write(" >>> Bad pixels replaced.\n")
        # Subtract the sky spectrum if requested by the user.
        if motes_parameters["-SUBTRACT_SKY"]:
            frame_dict, moffat_parameters_all_sky_bins, sky_extraction_limits = sky_locator(
                frame_dict, axes_dict, data_scaling_factor, header_parameters, bin_parameters, motes_parameters
            )
        # Will plot the location of the bins determined by get_bins if -DIAG_PLOT_BIN_LOC=1 in
        # motesparams.txt
        bin_parameters, frame_dict = common.get_bins(
            frame_dict,
            int(np.floor(binning_region_spatial_floor)),
            int(np.ceil(binning_region_spatial_ceiling)),
            axes_dict["dispersion_axis_length"],
            motes_parameters,
        )

        common.get_bins_output(
            bin_parameters,
            motes_parameters,
            binning_region_spatial_floor,
            binning_region_spatial_ceiling,
            frame_dict["data"],
            header_parameters,
            axes_dict,
        )

        # For each dispersion bin determined by the get_bins function, median the bin along the
        # dispersion axis, fit a moffat profile to the median data and then use the parameters of
        # the fitted Moffat function to localise the 2D spectrum.
        sys.stdout.write(
            " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
        )

        moffat_parameters_all_bins = []
        extraction_limits = []

        for each_bin in bin_parameters:
            # Take the median spatial profile of the dispersion bin, and leave out pixel columns in
            # the chip gaps if this is a GMOS spectrum.
            raw_bin_data = frame_dict["data"][:, each_bin[0] : each_bin[1]]
            chip_gap_location = np.where(np.median(raw_bin_data, axis=0) != 1)
            bin_data = np.nanmedian(raw_bin_data[:, chip_gap_location[0]], axis=1)

            # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median
            # spatial profile and return its parameters.
            bin_moffat_parameters = common.moffat_least_squares(
                axes_dict["spatial_axis"],
                bin_data * data_scaling_factor,
                header_parameters["seeing"],
                header_parameters["pixel_resolution"],
            )

            bin_moffat_parameters[0] /= data_scaling_factor
            bin_moffat_parameters[4] /= data_scaling_factor
            bin_moffat_parameters[5] /= data_scaling_factor

            # Define the extraction limits of the current dispersion bin based on the parameters of
            # the Moffat profile previously fitted to it.
            bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_fwhm, moffat_center = common.set_extraction_limits()
                bin_moffat_parameters,
                width_multiplier=motes_parameters["-FWHM_MULTIPLIER"],
            )

            extraction_limits.append([(each_bin[0] + each_bin[1]) * 0.5, bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_center])

            # Record the Moffat function parameters for each dispersion bin and add the wavstart
            # offset to the bin locations so they can be saved as metadata along with the extracted
            # spectrum.
            bin_moffat_parameters.append(each_bin[0] + axes_dict["wavelength_start"])
            bin_moffat_parameters.append(each_bin[1] + axes_dict["wavelength_start"])
            moffat_parameters_all_bins.append(bin_moffat_parameters)

            # DIAGNOSTICS - Plot computed moffat profile over data for each bin
            if motes_parameters["-DIAG_PLOT_MOFFAT"]:
                common.plot_fitted_spatial_profile(
                    axes_dict["spatial_axis"],
                    bin_data,
                    axes_dict["hi_resolution_spatial_axis"],
                    bin_moffat_parameters,
                    axes_dict["data_spatial_floor"],
                    header_parameters,
                )

        moffat_parameters_all_bins = np.array(moffat_parameters_all_bins)
        sys.stdout.write("     Fitting complete.\n")

        sys.stdout.write(" >>> Drawing extraction aperture limits. ")
        sys.stdout.flush()
        extraction_limits = np.array(extraction_limits).T

        # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum. All pixels
        # fully within the aperture are extracted.
        if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            draw_lines = [
                extraction_limits[0] + axes_dict["wavelength_start"],
                extraction_limits[1] + axes_dict["data_spatial_floor"] - 1,
                extraction_limits[0] + axes_dict["wavelength_start"],
                extraction_limits[2] + axes_dict["data_spatial_floor"] + 1,
            ]

            common.show_img(
                frame_dict["data"],
                axes_dict,
                header_parameters,
                draw_lines,
                "2D Spectrum Overplotted with Extraction Limits",
            )

        # Interpolate the extraction limits calculated for each median bin such that each
        # wavelength element across the entire unbinned wavelength axis of the entire 2D spectrum
        # has its own extraction limits.

        final_extraction_limits = common.interpolate_extraction_lims(
            extraction_limits, axes_dict["dispersion_axis_length"]
        )
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the
        # ends of the wavelength axis. All pixels fully within the aperture are extracted.
        if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            draw_lines = [
                np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
                final_extraction_limits[0] + axes_dict["data_spatial_floor"] - 1,
                np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavelength_start"],
                final_extraction_limits[1] + axes_dict["data_spatial_floor"] + 1,
            ]

            common.show_img(
                frame_dict["data"],
                axes_dict,
                header_parameters,
                draw_lines,
                "2D Spectrum Overplotted with Full Extraction Limits",
            )

        # Extract the spectrum from a supersampled version of the 2D image using the extraction
        # limits.
        frame_dict["data"] = frame_dict["data"].T
        frame_dict["errs"] = frame_dict["errs"].T
        frame_dict["qual"] = frame_dict["qual"].T

        optimal_1d_data, optimal_1d_errs, aperture_1d_data, aperture_1d_errs = common.optimal_extraction(
            frame_dict["data"],
            frame_dict["errs"],
            final_extraction_limits,
            moffat_parameters_all_bins,
            axes_dict,
        )

        final_extraction_limits = np.array(final_extraction_limits)
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot extracted spectrum.
        if motes_parameters["-PLOT_EXTRACTED_SPECTRUM"]:
            # DIAGNOSTICS, EXTRACTED SPECTRUM
            plt.figure(figsize=(9, 6))
            plt.errorbar(
                axes_dict["wavelength_axis"],
                aperture_1d_data,
                yerr=aperture_1d_errs,
                color="k",
                marker=".",
                label="Aperture Spectrum",
            )
            plt.errorbar(
                axes_dict["wavelength_axis"],
                optimal_1d_data,
                yerr=optimal_1d_errs,
                color="r",
                marker=".",
                label="Optimal Spectrum",
            )
            plt.grid(alpha=0.5, linestyle="dotted")
            plt.title("Extracted 1D Spectrum")
            plt.ylabel("Flux, " + header_parameters["flux_unit"])
            plt.xlabel("Wavelength, " + header_parameters["wavelength_unit"])
            plt.legend()
            plt.show()

        # Save the extracted spectrum to a new .fits file.
        if motes_parameters["-SAVE"]:
            sys.stdout.write(" >>> Saving 1D spectrum and metadata.\n")
            sys.stdout.flush()
            save_fits(
                axes_dict,
                header_parameters,
                optimal_1d_data,
                optimal_1d_errs,
                aperture_1d_data,
                aperture_1d_errs,
                input_file_primary_header,
                motes_parameters,
                input_file_path,
                moffat_profile_parameters,
                frame_dict,
                moffat_parameters_all_bins,
                final_extraction_limits,
                moffat_parameters_all_sky_bins,
                sky_extraction_limits,
            )

        sys.stdout.write(" >>> Extraction of " + input_file_path + " completed.\n")

    sys.stdout.write(" >>> MOTES Processing Complete.\n\n")
    return None


def save_fits(
    axes_dict,
    header_parameters,
    optimal_1d_data,
    optimal_1d_errs,
    aperture_1d_data,
    aperture_1d_errs,
    input_file_primary_header,
    motes_parameters,
    input_file_path,
    moffat_profile_parameters,
    frame_dict,
    moffat_parameters_all_bins,
    extraction_limits,
    moffat_parameters_all_sky_bins,
    sky_extraction_limits,
):
    """
    This function saves the extracted spectrum and intermediate products in a single, newly
    constructed, FITS file.

    Args:
        axes_dict (dict)                        : A dictionary containing the axes information.
        header_parameters (dict)                       : A dictionary containing the header information.
        optimal_1d_data (numpy.ndarray)               : An array containing the flux values of the optimally
                                               extracted 1D spectrum.
        optimal_1d_errs (numpy.ndarray)               : An array containing the flux errors of the optimally
                                               extracted 1D spectrum.
        aperture_1d_data (numpy.ndarray)               : An array containing the flux values of the aperture
                                               extracted 1D spectrum.
        aperture_1d_errs (numpy.ndarray)               : An array containing the flux errors of the aperture
                                               extracted 1D spectrum.
        input_file_primary_header (astropy.io.fits.header.Header) : The original FITS header of the 2D spectrum.
        motes_parameters (dict)                          : A dictionary containing the MOTES parameters.
        input_file_path (str)                       : The filename of the 1D spectrum.
        moffat_profile_parameters (list)                      : A list containing the Moffat fit parameters.
        frame_dict (dict)                         : A dictionary containing the original 2D spectrum
                                               data and error frames.
        moffat_parameters_all_bins (numpy.ndarray)                : A dictionary containing the binning parameters.
        extraction_limits (numpy.ndarray)       : An array containing the extraction limits.
        moffat_parameters_all_sky_bins (numpy.ndarray)               : An array containing the binning parameters for the
                                               sky extraction.
        sky_extraction_limits (list)             : A list containing the extraction limits for the sky
                                               extraction.

    Returns:
        None
    """

    input_file_primary_header["MOTES"] = "######## Extracted 1D Spectrum Metadata ########"
    input_file_primary_header.add_blank("", before="MOTES")
    input_file_primary_header["HIERARCH UTC EXT DATE"] = (
        datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
        "file creation date",
    )
    input_file_primary_header["HIERARCH SPATPIXL"] = (
        axes_dict["data_spatial_floor"],
        "lower limit of spatial axis, pix",
    )
    input_file_primary_header["HIERARCH SPATPIXH"] = (
        axes_dict["data_spatial_ceiling"],
        "upper limit of spatial axis, pix",
    )
    input_file_primary_header["HIERARCH DISPPIXL"] = (
        axes_dict["wavelength_start"],
        "lower limit of dispersion axis, pix",
    )
    input_file_primary_header["HIERARCH DISPPIXH"] = (
        axes_dict["wavelength_end"],
        "upper limit of dispersion axis, pix",
    )
    input_file_primary_header["HIERARCH WAVL"] = (
        np.floor(axes_dict["wavelength_axis"][0]),
        "lower limit of wav range, " + header_parameters["wavelength_unit"],
    )
    input_file_primary_header["HIERARCH WAVH"] = (
        np.ceil(axes_dict["wavelength_axis"][-1]),
        "upper limit of wav range, " + header_parameters["wavelength_unit"],
    )
    input_file_primary_header["HIERARCH WAVU"] = header_parameters["wavelength_unit"], "Wavelength unit"

    input_file_primary_header["HIERARCH MOFF A"] = round(moffat_profile_parameters[0], 5), "moffat profile amplitude"
    input_file_primary_header.add_blank(
        "Parameters fit to the median spatial profile of the spectrum",
        before="HIERARCH MOFF A",
    )
    input_file_primary_header["HIERARCH MOFF C"] = (
        round(moffat_profile_parameters[1] + axes_dict["data_spatial_floor"], 5),
        "moffat profile center",
    )
    input_file_primary_header["HIERARCH MOFF ALPHA"] = (
        round(moffat_profile_parameters[2], 5),
        "moffat profile alpha value",
    )
    input_file_primary_header["HIERARCH MOFF BETA"] = (
        round(moffat_profile_parameters[3], 5),
        "moffat profile beta value",
    )
    input_file_primary_header["HIERARCH MOFF BACK"] = (
        round(moffat_profile_parameters[4], 5),
        "moffat profile background level",
    )
    input_file_primary_header["HIERARCH MOFF GRAD"] = (
        round(moffat_profile_parameters[5], 5),
        "moffat profile background slope",
    )
    input_file_primary_header["HIERARCH IQ"] = (
        round(header_parameters["seeing"] * header_parameters["pixel_resolution"], 2),
        'IQ measured from median profile, "',
    )

    input_file_primary_header["HIERARCH SNR BIN LIMIT"] = motes_parameters["-SNR_BIN_LIM"], "maximum SNR per bin"
    input_file_primary_header.add_blank(
        "Dispersion Binning and Spectrum Extraction",
        before="HIERARCH SNR BIN LIMIT",
    )
    input_file_primary_header["HIERARCH COL BIN LIMIT"] = (
        int(motes_parameters["-COL_BIN_LIM"]),
        "minimum number of columns per bin",
    )
    input_file_primary_header["HIERARCH FWHM MULTIPLIER"] = (
        motes_parameters["-FWHM_MULTIPLIER"],
        "FWHM used to define the extraction limits",
    )

    if motes_parameters["-SUBTRACT_SKY"]:
        input_file_primary_header["HIERARCH SKYSUB FWHM MULT"] = (
            motes_parameters["-BG_FWHM_MULTIPLIER"],
            "FWHM multiplier for defining background",
        )
        input_file_primary_header.add_blank("Sky Subtraction", before="HIERARCH SKYSUB FWHM MULT")
        input_file_primary_header["HIERARCH SKYSUB SNR BIN LIM"] = (
            motes_parameters["-SKY_SNR_BIN_LIM"],
            "max SNR per bin for sky subtraction",
        )
        sky_model_hdu = fits.ImageHDU(frame_dict["sky_model"])
        sky_model_hdu.header["EXTNAME"] = "2D_SKY"
        sky_bin_hdu = fits.ImageHDU(moffat_parameters_all_sky_bins)
        sky_bin_hdu.header["EXTNAME"] = "SKY_BIN_PARS"
        sky_extraction_limits = fits.ImageHDU(sky_extraction_limits)
        sky_extraction_limits.header["EXTNAME"] = "SKY_EXT_LIMS"

    input_file_primary_header["HIERARCH EXTRACTED HDU ROW 0"] = "Wavelength Axis, " + header_parameters["wavelength_unit"]
    input_file_primary_header.add_blank(
        "Data Saved in the Extracted Spectrum HDU",
        before="HIERARCH EXTRACTED HDU ROW 0",
    )
    input_file_primary_header["HIERARCH EXTRACTED HDU ROW 1"] = "Flux, " + header_parameters["flux_unit"]
    input_file_primary_header["HIERARCH EXTRACTED HDU ROW 2"] = "Flux Uncertainty, " + header_parameters["flux_unit"]
    input_file_primary_header["EXTNAME"] = "OPTI_1D_SPEC"

    optimal_1d_datahdu = fits.PrimaryHDU([axes_dict["wavelength_axis"], optimal_1d_data, optimal_1d_errs], header=input_file_primary_header)
    aperture_1d_datahdu = fits.ImageHDU([axes_dict["wavelength_axis"], aperture_1d_data, aperture_1d_errs], header=input_file_primary_header)
    aperture_1d_datahdu.header["EXTNAME"] = "APER_1D_SPEC"
    orig_2d_spec_hdu = fits.ImageHDU(frame_dict["original_data"])
    orig_2d_spec_hdu.header["EXTNAME"] = "ORIG_2D_SPEC"
    orig_2d_errs_hdu = fits.ImageHDU(frame_dict["original_errs"])
    orig_2d_errs_hdu.header["EXTNAME"] = "ORIG_2D_ERRS"
    orig_2d_qual_hdu = fits.ImageHDU(frame_dict["original_qual"])
    orig_2d_qual_hdu.header["EXTNAME"] = "ORIG_2D_QUAL"
    bins_moffat_parameters_hdu = fits.ImageHDU(moffat_parameters_all_bins)
    bins_moffat_parameters_hdu.header["EXTNAME"] = "EXT_BIN_PARS"
    extraction_limits = fits.ImageHDU(extraction_limits)
    extraction_limits.header["EXTNAME"] = "EXT_LIMS"
    hdu_list = [
        optimal_1d_datahdu,
        aperture_1d_datahdu,
        orig_2d_spec_hdu,
        orig_2d_errs_hdu,
        orig_2d_qual_hdu,
        bins_moffat_parameters_hdu,
        extraction_limits,
    ]

    if motes_parameters["-SUBTRACT_SKY"]:
        hdu_list.append(sky_model_hdu)
        hdu_list.append(sky_bin_hdu)
        hdu_list.append(sky_extraction_limits)

    fits_hdu_list = fits.HDUList(hdu_list)
    fits_hdu_list.writeto("m" + input_file_path.split("/")[-1])
    fits_hdu_list.close()

    sys.stdout.write(" >>> Spectrum extracted and saved:\n")
    sys.stdout.write(
        "     " + "/".join(input_file_path.split("/")[0:-1]) + "/m" + input_file_path.split("/")[-1] + "\n"
    )
    return None


def sky_locator(frame_dict, axes_dict, data_scaling_factor, header_parameters, bin_parameters, motes_parameters):
    """
    Perform sky subtraction on the 2D spectrum. Locaalise the spectrum in the same way done for the
    extraction, and then use the regions outside the boundaries defined by that process to
    characterise and subtract background sky emission.

    Args:
        frame_dict (dict)  : A dictionary containing the 2D spectrum and its associated errors and
                            quality arrays.
        axes_dict (dict)   : A dictionary containing the wavelength and spatial axes of the 2D
                            spectrum.
        data_scaling_factor (float) : A flux scale factor to convert the flux units of the 2D spectrum to the
                            same units as the sky.
        header_parameters (dict) : A dictionary containing the header parameters of the 2D spectrum.
        bin_parameters (dict)  : A dictionary containing the bin parameters of the 2D spectrum.
        motes_parameters (dict)  : A dictionary containing the parameters of the extraction.

    Returns:
        frame_dict (dict)         : A dictionary containing the 2D spectrum and its associated
                                   errors and quality arrays.
        moffat_parameters_all_sky_bins (list)            : A list containing bins for the sky background.
        sky_extraction_limits (list) : A list containing the extraction limits for the sky background.
    """

    sys.stdout.write(
        " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
    )

    moffat_parameters_all_sky_bins = []
    extraction_limits = []
    for each_bin in bin_parameters:
        # Take the median spatial profile of the dispersion bin, and leave out pixel columns in the
        # chip gaps if this is a GMOS spectrum.
        raw_bin_data = frame_dict["data"][:, each_bin[0] : each_bin[1]]
        chip_gap_location = np.where(np.median(raw_bin_data, axis=0) != 1)
        bin_data = np.nanmedian(raw_bin_data[:, chip_gap_location[0]], axis=1)

        # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median
        # spatial profile and return its parameters.
        bin_moffat_parameters = common.moffat_least_squares(
            axes_dict["spatial_axis"],
            bin_data * data_scaling_factor,
            header_parameters["seeing"],
            header_parameters["pixel_resolution"],
        )

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle. The shape of
        # the profile fitted to the scaled spatial profile is the same as the unscaled, but to to
        # get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled down again after
        # the fitting.
        bin_moffat_parameters[0] /= data_scaling_factor
        bin_moffat_parameters[4] /= data_scaling_factor
        bin_moffat_parameters[5] /= data_scaling_factor

        # Define the extraction limits of the current dispersion bin based on the parameters of the
        # Moffat profile previously fitted to it.
        bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_fwhm, moffat_center = common.set_extraction_limits()
            bin_moffat_parameters,
            width_multiplier=motes_parameters["-BG_FWHM_MULTIPLIER"],
        )

        extraction_limits.append([(each_bin[0] + each_bin[1]) * 0.5, bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_center])

        # Record the Moffat function parameters for each dispersion bin and add the wavstart offset
        # to the bin locations so they can be saved as metadata along with the extracted spectrum.
        bin_moffat_parameters.append(each_bin[0] + axes_dict["wavelength_start"])
        bin_moffat_parameters.append(each_bin[1] + axes_dict["wavelength_start"])
        moffat_parameters_all_sky_bins.append(bin_moffat_parameters)

        # DIAGNOSTICS - Plot computed moffat profile over data for each bin
        if motes_parameters["-DIAG_PLOT_MOFFAT"]:
            common.plot_fitted_spatial_profile(
                axes_dict["spatial_axis"],
                bin_data,
                axes_dict["hi_resolution_spatial_axis"],
                bin_moffat_parameters,
                axes_dict["data_spatial_floor"],
                header_parameters,
            )

    moffat_parameters_all_sky_bins = np.array(moffat_parameters_all_sky_bins)

    sys.stdout.write("     Fitting complete.\n")

    sys.stdout.write(" >>> Drawing target/sky boundaries. ")
    sys.stdout.flush()
    extraction_limits = np.array(extraction_limits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
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

    # Interpolate the extraction limits calculated for each median bin such that each wavelength
    # element across the entire unbinned wavelength axis of the entire 2D spectrum has its own
    # extraction limits.
    sky_extraction_limits = common.interpolate_extraction_lims(
        extraction_limits, axes_dict["dispersion_axis_length"]
    )

    sys.stdout.write("DONE.\n")

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the
    # ends of the wavelength axis.
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

    frame_dict = common.subtract_sky(
        sky_extraction_limits[0],
        sky_extraction_limits[1],
        frame_dict,
        axes_dict,
        motes_parameters,
        header_parameters,
    )

    sys.stdout.write("\n     DONE.\n")
    sys.stdout.flush()

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the
    # ends of the wavelength axis.
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


if __name__ == "__main__":
    motes()
