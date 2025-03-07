#!/usr/bin/env python3

"""
	core.py

	Copyright (C) 2025 Tom Seccull & Dominik Kiersz
	
	This module is part of the MOTES package hosted at 
	https://github.com/tseccull/motes
	https://doi.org/####################################################
	
	If used, please cite the MOTES DOI above.
	
	This script is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	Last updated - 2025-03-07

	Description---------------------------------------------------------
	`core.py` contains the core `motes()` function that runs MOTES. 
"""

import argparse
import copy
import importlib.metadata
import logging
import logging.config
import motes.diagnostics as diagnostics
import motes.extraction as extraction
import motes.logs as logs
import motes.sky as sky
import motes.tracing as tracing
import motes.io.motesio as motesio
import numpy as np
import os


def motes(motes_args):
    """
    Run MOTES. This function is the entrypoint for MOTES. It loads .fits
    files from the input subdirectory and harvests the 2D data frames.
    It estimates bins for the spectrum. It then performs optional sky
    subtraction and optimal extraction of 1D spectra.
    
    Args:
     -- motes_args (Namespace)
          Namespace containing the keyword args required by MOTES. These
          may be collected from the command line with argparse or may be
          fed directly to MOTES via a script.
          
    Returns: None
    """

    cwd = os.getcwd()
    
    logging.config.dictConfig(logs.motes_logs(cwd, motes_args.verbose))
    logger = logging.getLogger("motes")
    
    logger.info("MOTES function and logger initialised.")
    logger.info("MOTES version is v" + importlib.metadata.version("motes"))
    logger.info("Current working directory is " + cwd)
    [logger.info("motes_args.%s = %s", key, val) for key, val in vars(motes_args).items()]
    
    files_and_regions = motesio.read_regions(cwd)

    for regions_and_file in files_and_regions:
        input_file_name = regions_and_file[-1][:-5]
        
        logger.info("MOTES processing file: ./" + input_file_name + ".fits")
        logger.info("Switching to log file: ./" + input_file_name + "_motes.log")
        logger.handlers[1].setStream(open(cwd + "/" + input_file_name + "_motes.log", "w"))
        logger.info("Begin MOTES processing for ./inputs/" + input_file_name + ".fits")
        
        input_file_path = cwd + "/inputs/" + regions_and_file[-1]
        data_regions = regions_and_file[:4]

        header_parameters, frame_dict, axes_dict, original_hdu_list = (
            motesio.data_harvest(input_file_path, data_regions)
        )
        
        frame_dict["original_data"] = copy.deepcopy(frame_dict["data"])
        frame_dict["original_errs"] = copy.deepcopy(frame_dict["errs"])
        logger.info("Relevant image frames and header data harvested from input file.")

        # Perform initial least-squares Moffat fitting over the entire
        # 2D spectrum collapsed along the dispersion axis with a median.
        logger.info("Fitting Moffat profile to median spatial profile of entire spectrum.")
        full_median_spatial_profile = np.nanmedian(frame_dict["data"], axis=1)
        moffat_profile_parameters, data_scaling_factor = (
            tracing.median_moffat_fitting(
                full_median_spatial_profile,
                header_parameters, 
                axes_dict["spatial_axis"]
            ) 
        )
        iq_pixel = round(header_parameters["seeing"], 2)
        iq = round(header_parameters["seeing"] * header_parameters["pixel_resolution"], 2)
        logger.info("FWHM of median spatial profile is %s Pixels, or %s''.", iq_pixel, iq)

        # Use the parameters of the Moffat profile fitted to the median
        # spatial profile of the entire spectrum to determine spatial
        # limits that are used to bound the region of the spectrum used
        # by the common.get_bins function to bin the 2D spectrum while
        # taking account of its S/N.
        limits_and_center = tracing.set_extraction_limits(moffat_profile_parameters)
        binning_region_spatial_floor = limits_and_center[0]
        binning_region_spatial_ceiling = limits_and_center[1]
        moffat_center = limits_and_center[2]
        
        logger.info(
            "Spectrum localised to aperture covering spatial pixel rows %s-%s",
            int(binning_region_spatial_floor + axes_dict["data_spatial_floor"]),
            int(binning_region_spatial_ceiling + axes_dict["data_spatial_floor"])
        )
        
        logger.info("Median Moffat function parameters:")
        logger.info("     A = %s", moffat_profile_parameters[0])
        logger.info("     c = %s", moffat_profile_parameters[1] + axes_dict["data_spatial_floor"])
        logger.info(" alpha = %s", moffat_profile_parameters[2])
        logger.info("  beta = %s", moffat_profile_parameters[3])
        logger.info("     B = %s", moffat_profile_parameters[4])
        logger.info("     m = %s", moffat_profile_parameters[5])

        if motes_args.diag_plot_collapsed:
            png_name = "m" + input_file_name + "_moffat_median.png"
            diagnostics.plot_fitted_spatial_profile(
                axes_dict["spatial_axis"],
                full_median_spatial_profile,
                axes_dict["hi_resolution_spatial_axis"],
                moffat_profile_parameters,
                axes_dict["data_spatial_floor"],
                header_parameters,
                png_name
            )

        if motes_args.subtract_sky:
            logger.info("Beginning sky location and subtraction.")
            # Determine the location of bins on the dispersion axis
            # within which to measure the spatial profile.
            bin_parameters, frame_dict = tracing.get_bins(
                frame_dict,
                int(np.floor(binning_region_spatial_floor)),
                int(np.ceil(binning_region_spatial_ceiling)),
                axes_dict["dispersion_axis_length"],
                motes_args,
                has_sky=True,
            )
            logger.info("%s sky localisation bins defined on dispersion axis", len(bin_parameters))
            
            # Will plot the location of the bins determined by get_bins
            # if -DIAG_PLOT_BIN_LOC=1 in motesparams.txt
            if motes_args.diag_plot_bins:
                diagnostics.get_bins_output(
                    bin_parameters,
                    motes_args,
                    binning_region_spatial_floor,
                    binning_region_spatial_ceiling,
                    frame_dict["data"],
                    header_parameters,
                    axes_dict,
                )
            
            # Subtract the sky spectrum if requested by the user.
            frame_dict, moffat_parameters_all_sky_bins, sky_extraction_limits = (
                sky.sky_locator(
                    frame_dict,
                    axes_dict,
                    data_scaling_factor,
                    header_parameters,
                    bin_parameters,
                    motes_args,
                    input_file_name
                )
            )
            logger.info("Sky subtraction complete.")

        logger.info("Beginning spectrum localisation and extraction.")
        bin_parameters, frame_dict = tracing.get_bins(
            frame_dict,
            int(np.floor(binning_region_spatial_floor)),
            int(np.ceil(binning_region_spatial_ceiling)),
            axes_dict["dispersion_axis_length"],
            motes_args,
        )

        # Will plot the location of the bins determined by get_bins if
        # -DIAG_PLOT_BIN_LOC=1 in motesparams.txt
        if motes_args.diag_plot_bins:
            diagnostics.get_bins_output(
                bin_parameters,
                motes_args,
                binning_region_spatial_floor,
                binning_region_spatial_ceiling,
                frame_dict["data"],
                header_parameters,
                axes_dict,
            )

        # For each dispersion bin determined by the get_bins function,
        # median the bin along the dispersion axis, fit a moffat profile
        # to the median data and then use the parameters of the fitted
        # Moffat function to localise the 2D spectrum.
        logger.info("Fitting Moffat functions to 2D data to localise the spectrum.")

        moffat_parameters_all_bins = []
        extraction_limits = []

        for each_bin in bin_parameters:
            bin_moffat_parameters, extraction_limits, bin_data = (
                tracing.moffat_fitting(
                    each_bin,
                    frame_dict["data"],
                    axes_dict,
                    data_scaling_factor,
                    header_parameters,
                    motes_args.sky_fwhm_multiplier,
                    extraction_limits
                )
            )
            moffat_parameters_all_bins.append(bin_moffat_parameters)

            # DIAGNOSTICS - Plot computed moffat profile over data for
            # each bin
            if motes_args.diag_plot_moffat:
                png_name = diagnostics.get_png_name(
                    each_bin, axes_dict["wavelength_start"], input_file_name, "ext"
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

        moffat_parameters_all_bins = np.array(moffat_parameters_all_bins)
        logger.info("2D spectrum fitting completed.")
        logger.info("Drawing extraction aperture limits.")
        extraction_limits = np.array(extraction_limits).T

        # DIAGNOSTICS - Plot the determined extraction limits over the
        # 2D spectrum. All pixels fully within the aperture are
        # extracted.
        if motes_args.diag_plot_localisation:
            draw_lines = diagnostics.get_draw_lines(axes_dict, extraction_limits)
            title = "2D Spectrum Overplotted with Extraction Limits"
            data_frame = frame_dict["data"]
            diagnostics.show_img(data_frame, axes_dict, header_parameters, draw_lines, title)

        # Interpolate the extraction limits calculated for each median
        # bin such that each wavelength element across the entire
        # unbinned wavelength axis of the entire 2D spectrum has its own
        # extraction limits.
        final_extraction_limits = (
            tracing.interpolate_extraction_lims(
                extraction_limits, axes_dict["dispersion_axis_length"]
            )
        )

        # DIAGNOSTICS - Plot the final extraction limits including the
        # extrapolated sections at the ends of the wavelength axis. All
        # pixels fully within the aperture are extracted.
        if motes_args.diag_plot_localisation:
            draw_lines = diagnostics.get_draw_lines(axes_dict, final_extraction_limits)
            title = "2D Spectrum Overplotted with Full Extraction Limits"
            data_frame = frame_dict["data"]
            diagnostics.show_img(data_frame, axes_dict, header_parameters, draw_lines, title)

        # Extract the spectrum from a supersampled version of the 2D
        # image using the extraction limits.
        frame_dict["data"] = frame_dict["data"].T
        frame_dict["errs"] = frame_dict["errs"].T
        frame_dict["qual"] = frame_dict["qual"].T
        
        extracted = extraction.optimal_extraction(
            frame_dict["data"],
            frame_dict["errs"],
            final_extraction_limits,
            moffat_parameters_all_bins,
            axes_dict,
        )
        logger.info("Extraction completed.")

        final_extraction_limits = np.array(final_extraction_limits)
       

        # DIAGNOSTICS - Plot extracted spectrum.
        if motes_args.diag_plot_spectrum:
            diagnostics.plot_spectrum(axes_dict["wavelength_axis"], extracted, header_parameters)

        # Save the extracted spectrum to a new .fits file.
        if motes_args.save:
            motesio.save_fits(
                original_hdu_list,
                axes_dict,
                header_parameters,
                extracted,
                motes_args,
                input_file_name,
                moffat_profile_parameters,
                frame_dict,
                moffat_parameters_all_bins,
                final_extraction_limits,
                moffat_parameters_all_sky_bins,
                sky_extraction_limits,
            )

        logger.info("Processing %s completed.", input_file_path)
        logger.handlers[1].setStream(open(cwd + "/" + "motes.log", "a"))

    logger.info("MOTES Processing completed.")
    return None

