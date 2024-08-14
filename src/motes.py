#!/home/tom/miniforge3/envs/work/bin/python

"""
MOTES Modular and Optimal Tracer and Extractor of Spectra.

Description: Modular and Optimal Tracer and Extractor of Specrtra 
(MOTES). A Python package for extracting spectrum from astronomical
2D spectrograms.

Version: 0.5.0
Date: 2024-08-12
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
import motesio.motesio as motesio
import motessky.sky as sky


def motes():
    """
    Run MOTES. This function is the entrypoint for the MOTES. It loads .fits files from the input
    subdirectory and harvests the 2D data frames. It estimates bins for the spectrum. It then
    performs optional sky subtraction and optimal extraction of 1D spectra.
    """
    # Run startup functions
    motes_parameters = motesio.read_motes_parameter_file()  # Import parameters from file to dict
    data_regions = motesio.read_regions()  # Search for, and read in, reg.txt

    # Open and process each spectrum contained in the current directory.
    for i, input_file_path in enumerate(sorted(glob.glob("./inputs/*.fits"))):
        sys.stdout.write(("/" * (70 - len(input_file_path[:70]))) + " " + input_file_path[:70] + "\n")
        sys.stdout.write(" >>> Beginning MOTES Processing\n")

        # Gather header metadata and the image data from the 2D image file.
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file.\n"
        )
        header_parameters, frame_dict, axes_dict, original_hdu_list = (
            motesio.data_harvest(i, input_file_path, data_regions)
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
        limits_and_center = common.set_extraction_limits(moffat_profile_parameters)
        binning_region_spatial_floor = limits_and_center[0]
        binning_region_spatial_ceiling = limits_and_center[1]
        moffat_center = limits_and_center[2]
        
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
            frame_dict, moffat_parameters_all_sky_bins, sky_extraction_limits = sky.sky_locator(
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
            bin_moffat_parameters[1] += axes_dict["data_spatial_floor"]
            bin_moffat_parameters[4] /= data_scaling_factor
            bin_moffat_parameters[5] /= data_scaling_factor

            # Define the extraction limits of the current dispersion bin based on the parameters of
            # the Moffat profile previously fitted to it.
            limits_and_center = common.set_extraction_limits(
                bin_moffat_parameters,
                width_multiplier=motes_parameters["-FWHM_MULTIPLIER"],
            )
            bin_lower_extraction_limit = limits_and_center[0]
            bin_upper_extraction_limit = limits_and_center[1]
            moffat_center = limits_and_center[2]

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
            motesio.save_fits(
                original_hdu_list,
                axes_dict,
                header_parameters,
                optimal_1d_data,
                optimal_1d_errs,
                aperture_1d_data,
                aperture_1d_errs,
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


if __name__ == "__main__":
    motes()
