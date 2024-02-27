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
    data_region = startup.read_regions()  # Search for, and read in, reg.txt

    # Open and process each spectrum contained in the current directory.
    for i, input_file_path in enumerate(sorted(glob.glob("./inputs/*.fits"))):
        sys.stdout.write(("/" * (70 - len(input_file_path[:70]))) + " " + input_file_path[:70] + "\n")
        sys.stdout.write(" >>> Beginning MOTES Processing\n")

        # Gather header metadata and the image data from the 2D image file.
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file.\n"
        )
        header_parameters, frame_dict, axes_dict, input_file_primary_header = harvester.data_harvest(
            i, input_file_path, data_region
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
        binning_region_spatial_floor, binning_region_spatial_ceiling, moffat_fwhm, moffat_center = common.extraction_limits(moffat_profile_parameters)
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
            sky=True,
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
            frame_dict, sky_bin_parameters, sky_extraction_limits = sky_locator(
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
            bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_fwhm, moffat_center = common.extraction_limits(
                bin_moffat_parameters,
                width_multiplier=motes_parameters["-FWHM_MULTIPLIER"],
            )

            extraction_limits.append([(each_bin[0] + each_bin[1]) * 0.5, bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_center])

            # Record the Moffat function parameters for each dispersion bin and add the wavstart
            # offset to the bin locations so they can be saved as metadata along with the extracted
            # spectrum.
            bin_moffat_parameters.append(bin[0] + axes_dict["wavstart"])
            bin_moffat_parameters.append(bin[1] + axes_dict["wavstart"])
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

        binpars = np.array(moffat_parameters_all_bins)
        sys.stdout.write("     Fitting complete.\n")

        sys.stdout.write(" >>> Drawing extraction aperture limits. ")
        sys.stdout.flush()
        extraction_limits = np.array(extraction_limits).T

        # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum. All pixels
        # fully within the aperture are extracted.
        if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            drawlines = [
                extraction_limits[0] + axes_dict["wavstart"],
                extraction_limits[1] + axes_dict["data_spatial_floor"] - 1,
                extraction_limits[0] + axes_dict["wavstart"],
                extraction_limits[2] + axes_dict["data_spatial_floor"] + 1,
            ]

            common.show_img(
                frame_dict["data"],
                axes_dict,
                header_parameters,
                drawlines,
                "2D Spectrum Overplotted with Extraction Limits",
            )

        # Interpolate the extraction limits calculated for each median bin such that each
        # wavelength element across the entire unbinned wavelength axis of the entire 2D spectrum
        # has its own extraction limits.

        finalextractionlims = common.interpolate_extraction_lims(
            extraction_limits, axes_dict["dispersion_axis_length"]
        )
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the
        # ends of the wavelength axis. All pixels fully within the aperture are extracted.
        if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            drawlines = [
                np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
                finalextractionlims[0] + axes_dict["data_spatial_floor"] - 1,
                np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
                finalextractionlims[1] + axes_dict["data_spatial_floor"] + 1,
            ]

            common.show_img(
                frame_dict["data"],
                axes_dict,
                header_parameters,
                drawlines,
                "2D Spectrum Overplotted with Full Extraction Limits",
            )

        # Extract the spectrum from a supersampled version of the 2D image using the extraction
        # limits.
        sys.stdout.write(" >>> Extracting 1D spectrum. ")
        sys.stdout.flush()
        frame_dict["data"] = frame_dict["data"].T
        frame_dict["errs"] = frame_dict["errs"].T
        frame_dict["qual"] = frame_dict["qual"].T

        opdata1D, operrs1D, apdata1D, aperrs1D = common.optimal_extraction(
            frame_dict["data"],
            frame_dict["errs"],
            finalextractionlims,
            binpars,
            axes_dict,
        )

        finalextractionlims = np.array(finalextractionlims)
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot extracted spectrum.
        if motes_parameters["-PLOT_EXTRACTED_SPECTRUM"]:
            # DIAGNOSTICS, EXTRACTED SPECTRUM
            plt.figure(figsize=(9, 6))
            plt.errorbar(
                axes_dict["waxis"],
                apdata1D,
                yerr=aperrs1D,
                color="k",
                marker=".",
                label="Aperture Spectrum",
            )
            plt.errorbar(
                axes_dict["waxis"],
                opdata1D,
                yerr=operrs1D,
                color="r",
                marker=".",
                label="Optimal Spectrum",
            )
            plt.grid(alpha=0.5, linestyle="dotted")
            plt.title("Extracted 1D Spectrum")
            plt.ylabel("Flux, " + header_parameters["fluxunit"])
            plt.xlabel("Wavelength, " + header_parameters["wavunit"])
            plt.legend()
            plt.show()

        # Save the extracted spectrum to a new .fits file.
        if motes_parameters["-SAVE"]:
            sys.stdout.write(" >>> Saving 1D spectrum and metadata.\n")
            sys.stdout.flush()
            save_fits(
                axes_dict,
                header_parameters,
                opdata1D,
                operrs1D,
                apdata1D,
                aperrs1D,
                input_file_primary_header,
                motes_parameters,
                input_file_path,
                moffat_profile_parameters,
                frame_dict,
                binpars,
                finalextractionlims,
                sky_bin_parameters,
                sky_extraction_limits,
            )

        sys.stdout.write(" >>> Extraction of " + input_file_path + " completed.\n")

    sys.stdout.write(" >>> MOTES Processing Complete.\n\n")
    return None


def save_fits(
    axdict,
    hparams,
    opflux,
    operrs,
    apflux,
    aperrs,
    head,
    pars,
    filename,
    moffpars,
    fdict,
    bpars,
    extractionlims,
    sbpars,
    skyextractionlims,
):
    """
    This function saves the extracted spectrum and intermediate products in a single, newly
    constructed, FITS file.

    Args:
        axdict (dict)                        : A dictionary containing the axes information.
        hparams (dict)                       : A dictionary containing the header information.
        opflux (numpy.ndarray)               : An array containing the flux values of the optimally
                                               extracted 1D spectrum.
        operrs (numpy.ndarray)               : An array containing the flux errors of the optimally
                                               extracted 1D spectrum.
        apflux (numpy.ndarray)               : An array containing the flux values of the aperture
                                               extracted 1D spectrum.
        aperrs (numpy.ndarray)               : An array containing the flux errors of the aperture
                                               extracted 1D spectrum.
        head (astropy.io.fits.header.Header) : The original FITS header of the 2D spectrum.
        pars (dict)                          : A dictionary containing the MOTES parameters.
        filename (str)                       : The filename of the 1D spectrum.
        moffpars (list)                      : A list containing the Moffat fit parameters.
        fdict (dict)                         : A dictionary containing the original 2D spectrum
                                               data and error frames.
        bpars (numpy.ndarray)                : A dictionary containing the binning parameters.
        extractionlims (numpy.ndarray)       : An array containing the extraction limits.
        sbpars (numpy.ndarray)               : An array containing the binning parameters for the
                                               sky extraction.
        skyextractionlims (list)             : A list containing the extraction limits for the sky
                                               extraction.

    Returns:
        None
    """

    head["MOTES"] = "######## Extracted 1D Spectrum Metadata ########"
    head.add_blank("", before="MOTES")
    head["HIERARCH UTC EXT DATE"] = (
        datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
        "file creation date",
    )
    head["HIERARCH SPATPIXL"] = (
        axdict["data_spatial_floor"],
        "lower limit of spatial axis, pix",
    )
    head["HIERARCH SPATPIXH"] = (
        axdict["imgend"],
        "upper limit of spatial axis, pix",
    )
    head["HIERARCH DISPPIXL"] = (
        axdict["wavstart"],
        "lower limit of dispersion axis, pix",
    )
    head["HIERARCH DISPPIXH"] = (
        axdict["wavend"],
        "upper limit of dispersion axis, pix",
    )
    head["HIERARCH WAVL"] = (
        np.floor(axdict["waxis"][0]),
        "lower limit of wav range, " + hparams["wavunit"],
    )
    head["HIERARCH WAVH"] = (
        np.ceil(axdict["waxis"][-1]),
        "upper limit of wav range, " + hparams["wavunit"],
    )
    head["HIERARCH WAVU"] = hparams["wavunit"], "Wavelength unit"

    head["HIERARCH MOFF A"] = round(moffpars[0], 5), "moffat profile amplitude"
    head.add_blank(
        "Parameters fit to the median spatial profile of the spectrum",
        before="HIERARCH MOFF A",
    )
    head["HIERARCH MOFF C"] = (
        round(moffpars[1] + axdict["data_spatial_floor"], 5),
        "moffat profile center",
    )
    head["HIERARCH MOFF ALPHA"] = (
        round(moffpars[2], 5),
        "moffat profile alpha value",
    )
    head["HIERARCH MOFF BETA"] = (
        round(moffpars[3], 5),
        "moffat profile beta value",
    )
    head["HIERARCH MOFF BACK"] = (
        round(moffpars[4], 5),
        "moffat profile background level",
    )
    head["HIERARCH MOFF GRAD"] = (
        round(moffpars[5], 5),
        "moffat profile background slope",
    )
    head["HIERARCH IQ"] = (
        round(hparams["seeing"] * hparams["pixel_resolution"], 2),
        'IQ measured from median profile, "',
    )

    head["HIERARCH SNR BIN LIMIT"] = pars["-SNR_BIN_LIM"], "maximum SNR per bin"
    head.add_blank(
        "Dispersion Binning and Spectrum Extraction",
        before="HIERARCH SNR BIN LIMIT",
    )
    head["HIERARCH COL BIN LIMIT"] = (
        int(pars["-COL_BIN_LIM"]),
        "minimum number of columns per bin",
    )
    head["HIERARCH FWHM MULTIPLIER"] = (
        pars["-FWHM_MULTIPLIER"],
        "FWHM used to define the extraction limits",
    )
    head["HIERARCH INTERP KIND"] = (
        pars["-INTERP_KIND"],
        "interpolation mode used",
    )

    if pars["-SUBTRACT_SKY"]:
        head["HIERARCH SKYSUB FWHM MULT"] = (
            pars["-BG_FWHM_MULTIPLIER"],
            "FWHM multiplier for defining background",
        )
        head.add_blank("Sky Subtraction", before="HIERARCH SKYSUB FWHM MULT")
        head["HIERARCH SKYSUB SNR BIN LIM"] = (
            pars["-SKY_SNR_BIN_LIM"],
            "max SNR per bin for sky subtraction",
        )
        skymodhdu = fits.ImageHDU(fdict["skymod"])
        skymodhdu.header["EXTNAME"] = "2D_SKY"
        skybinhdu = fits.ImageHDU(sbpars)
        skybinhdu.header["EXTNAME"] = "SKY_BIN_PARS"
        skyextractionlims = fits.ImageHDU(skyextractionlims)
        skyextractionlims.header["EXTNAME"] = "SKY_EXT_LIMS"

    head["HIERARCH EXTRACTED HDU ROW 0"] = "Wavelength Axis, " + hparams["wavunit"]
    head.add_blank(
        "Data Saved in the Extracted Spectrum HDU",
        before="HIERARCH EXTRACTED HDU ROW 0",
    )
    head["HIERARCH EXTRACTED HDU ROW 1"] = "Flux, " + hparams["fluxunit"]
    head["HIERARCH EXTRACTED HDU ROW 2"] = "Flux Uncertainty, " + hparams["fluxunit"]
    head["EXTNAME"] = "OPTI_1D_SPEC"

    opfluxhdu = fits.PrimaryHDU([axdict["waxis"], opflux, operrs], header=head)
    apfluxhdu = fits.ImageHDU([axdict["waxis"], apflux, aperrs], header=head)
    apfluxhdu.header["EXTNAME"] = "APER_1D_SPEC"
    spec2Dhdu = fits.ImageHDU(fdict["original_data"])
    spec2Dhdu.header["EXTNAME"] = "ORIG_2D_SPEC"
    errs2Dhdu = fits.ImageHDU(fdict["original_errs"])
    errs2Dhdu.header["EXTNAME"] = "ORIG_2D_ERRS"
    qual2Dhdu = fits.ImageHDU(fdict["ogqual"])
    qual2Dhdu.header["EXTNAME"] = "ORIG_2D_QUAL"
    binhdu = fits.ImageHDU(bpars)
    binhdu.header["EXTNAME"] = "EXT_BIN_PARS"
    extractionlims = fits.ImageHDU(extractionlims)
    extractionlims.header["EXTNAME"] = "EXT_LIMS"
    hdu_list = [
        opfluxhdu,
        apfluxhdu,
        spec2Dhdu,
        errs2Dhdu,
        qual2Dhdu,
        binhdu,
        extractionlims,
    ]

    if pars["-SUBTRACT_SKY"]:
        hdu_list.append(skymodhdu)
        hdu_list.append(skybinhdu)
        hdu_list.append(skyextractionlims)

    hdulist = fits.HDUList(hdu_list)
    filenamelist = filename.split("_")
    hdulist.writeto("_".join(filenamelist[0:-2]) + "_" + "1D" + "_" + filenamelist[-1])
    hdulist.close()

    sys.stdout.write(" >>> Spectrum extracted and saved:\n")
    sys.stdout.write(
        "_".join(filenamelist[0:-2]) + "_" + "1D" + "_" + filenamelist[-1] + "\n"
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
        skybin (list)            : A list containing bins for the sky background.
        skyextractionlims (list) : A list containing the extraction limits for the sky background.
    """

    sys.stdout.write(
        " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
    )

    skybin = []
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
        bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_fwhm, moffat_center = common.extraction_limits(
            bin_moffat_parameters,
            width_multiplier=motes_parameters["-BG_FWHM_MULTIPLIER"],
        )

        extraction_limits.append([(each_bin[0] + each_bin[1]) * 0.5, bin_lower_extraction_limit, bin_upper_extraction_limit, moffat_center])

        # Record the Moffat function parameters for each dispersion bin and add the wavstart offset
        # to the bin locations so they can be saved as metadata along with the extracted spectrum.
        bin_moffat_parameters.append(bin[0] + axes_dict["wavstart"])
        bin_moffat_parameters.append(bin[1] + axes_dict["wavstart"])
        skybin.append(bin_moffat_parameters)

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

    skybin = np.array(skybin)

    sys.stdout.write("     Fitting complete.\n")

    sys.stdout.write(" >>> Drawing target/sky boundaries. ")
    sys.stdout.flush()
    extraction_limits = np.array(extraction_limits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
    if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        drawlines = [
            extraction_limits[0] + axes_dict["wavstart"],
            (extraction_limits[1]) + axes_dict["data_spatial_floor"],
            extraction_limits[0] + axes_dict["wavstart"],
            (extraction_limits[2]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            drawlines,
            "2D Spectrum Overplotted with Target/Sky Boundaries",
        )

    # Interpolate the extraction limits calculated for each median bin such that each wavelength
    # element across the entire unbinned wavelength axis of the entire 2D spectrum has its own
    # extraction limits.
    skyextractionlims = common.interpolate_extraction_lims(
        extraction_limits, axes_dict["dispersion_axis_length"]
    )

    sys.stdout.write("DONE.\n")

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the
    # ends of the wavelength axis.
    if motes_parameters["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        drawlines = [
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
            (skyextractionlims[0]) + axes_dict["data_spatial_floor"],
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
            (skyextractionlims[1]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            drawlines,
            "2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )

    sys.stdout.write(" >>> Subtracting sky.\n")
    sys.stdout.flush()

    frame_dict = common.subtract_sky(
        skyextractionlims[0],
        skyextractionlims[1],
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
        drawlines = [
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
            (skyextractionlims[0]) + axes_dict["data_spatial_floor"],
            np.array(range(axes_dict["dispersion_axis_length"])) + axes_dict["wavstart"],
            (skyextractionlims[1]) + axes_dict["data_spatial_floor"],
        ]

        common.show_img(
            frame_dict["data"],
            axes_dict,
            header_parameters,
            drawlines,
            "Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )

    return frame_dict, skybin, skyextractionlims


if __name__ == "__main__":
    motes()
