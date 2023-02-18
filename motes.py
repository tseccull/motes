"""
MOTES: Modular and Optimal Tracer and Extractor of Spectra.
Description: Modular and Optimal Tracer and Extractor of Specrtra (MOTES). A Python package for extracting spectrum from astronomical 2D spectrograms.
Version: 0.4.2
Date: 2023-02-16
Authors: Tom Seccull, Dominik Kiersz
Licence: GNU General Public License v3.0
"""
import astropy.io.fits as fits
import copy
import datetime
import glob
import motes.common as common
import motes.harvester as harvester
import motes.startup as startup
import matplotlib.pyplot as plt
import numpy as np
import sys


def motes():
    """Run MOTES.

    Description:
        This function is the entrypoint for the MOTES.
        It loads .fits files from the input subdirectory and harvest the data.
        Estimates bins for the spectrum.
        It then performs optional sky extraction and optimal extraction.
    """
    # Run startup functions
    params = startup.read_parfile()  # Import parameters from file to dict
    intreg = startup.read_regions()  # Search for, and read in, reg.txt

    # Open and process each spectrum contained in the current directory.
    for i, file_2D in enumerate(sorted(glob.glob("./inputs/*.fits"))):
        sys.stdout.write(
            ("/" * (70 - len(file_2D[:70]))) + " " + file_2D[:70] + "\n"
        )
        sys.stdout.write(" >>> Beginning MOTES Processing\n")

        # Gather header metadata and the image data from the 2D image file.
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file.\n"
        )
        headparams, framedict, axesdict, imghead = harvester.data_harvest(
            i, file_2D, intreg, params
        )
        # Make backup copies of the original data and error frames.
        framedict["ogdata"] = copy.deepcopy(framedict["data"])
        framedict["ogerrs"] = copy.deepcopy(framedict["errs"])
        sys.stdout.write(
            " >>> Gathering image frames and header data from input file completed.\n"
        )

        # Perform initial least-squares Moffat fitting over the entire 2D spectrum collapsed along the dispersion axis with a median.
        sys.stdout.write(
            " >>> Fitting Moffat profile to median spatial profile of entire spectrum. "
        )
        sys.stdout.flush()

        # Calculate median spatial profile of the spectrum.
        datadispcollapse = np.nanmedian(framedict["data"], axis=1)
        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle.
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled down again after the fitting.
        datascale = 10 ** np.abs(
            np.floor(np.log10(np.abs(np.nanmedian(datadispcollapse))))
        )
        # Fit the median spatial profile with a Moffat function.
        moffparams = common.moffat_least_squares(
            axesdict["saxis"],
            datadispcollapse * datascale,
            headparams["seeing"],
            headparams["pixresolution"],
        )
        # Get an improved estimate of the FWHM of the spectrum from the best fit Moffat profile.
        headparams["seeing"] = (
            2 * moffparams[2] * np.sqrt((2 ** (1 / moffparams[3])) - 1)
        )
        # Scale the amplitude, background gradient, and background level of the model Moffat profile down.
        moffparams[0] /= datascale
        moffparams[4] /= datascale
        moffparams[5] /= datascale

        sys.stdout.write("DONE.\n")
        sys.stdout.write(
            " >>> FWHM of median spatial profile is "
            + str(round(headparams["seeing"], 2))
            + " Pixels, or "
            + str(round(headparams["seeing"] * headparams["pixresolution"], 2))
            + '"\n'
        )

        # Use the parameters of the Moffat profile fitted to the median spatial profile of the entire spectrum to determine
        # spatial limits that are used to bound the region of the spectrum used by the common.get_bins function to
        # to bin the 2D spectrum while taking account of its S/N.
        lowext, highext, fwhm, cent = common.extraction_limits(
            moffparams, axesdict
        )
        sys.stdout.write(
            " >>> Spectrum localised to aperture in range of spatial pixel rows "
            + str(int(lowext + axesdict["imgstart"]))
            + "-"
            + str(int(highext + axesdict["imgstart"]))
            + "\n"
        )

        # DIAGNOSTICS -  Plot fitted Moffat profile over collapsed 2D spectrum and print the parameters of the fitted Moffat profile.
        if params["-DIAG_PLOT_COLLAPSED_2D_SPEC"]:
            common.printmoffparams(moffparams, axesdict["imgstart"], datascale)
            common.plot_fitted_spatial_profile(
                axesdict["saxis"],
                datadispcollapse,
                axesdict["hrsaxis"],
                moffparams,
                axesdict["imgstart"],
                headparams,
            )

        # Determine the location of bins on the dispersion axis within which to measure the spatial profile.
        binparams, framedict = common.get_bins(
            framedict,
            int(np.floor(lowext)),
            int(np.ceil(highext)),
            axesdict["dispaxislen"],
            params,
            sky=True,
        )

        # Will plot the location of the bins determined by get_bins if -DIAG_PLOT_BIN_LOC=1 in motesparams.txt
        common.get_bins_output(
            binparams,
            params,
            lowext,
            highext,
            framedict["data"],
            headparams,
            axesdict,
        )
        sys.stdout.write(" >>> Bad pixels replaced.\n")
        # Subtract the sky spectrum if requested by the user.
        if params["-SUBTRACT_SKY"]:
            framedict, skybinpars, skyextlims = skyloc(
                framedict, axesdict, datascale, headparams, binparams, params
            )
        # Will plot the location of the bins determined by get_bins if -DIAG_PLOT_BIN_LOC=1 in motesparams.txt
        binparams, framedict = common.get_bins(
            framedict,
            int(np.floor(lowext)),
            int(np.ceil(highext)),
            axesdict["dispaxislen"],
            params,
            replace_crbp=bool(params["-REPLACE_CRBP"]),
        )

        common.get_bins_output(
            binparams,
            params,
            lowext,
            highext,
            framedict["data"],
            headparams,
            axesdict,
        )

        # For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
        # moffat profile to the median data and then use the parameters of the fitted Moffat function to localise the 2D
        # spectrum.
        sys.stdout.write(
            " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
        )

        extbin = []
        extractionlimits = []

        for bin in binparams:

            # Take the median spatial profile of the dispersion
            # bin, and leave out pixel columns in the chip gaps if this is a GMOS spectrum.
            binimg = framedict["data"][:, bin[0] : bin[1]]
            chipgap = np.where(np.median(binimg, axis=0) != 1)
            bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)

            # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median spatial profile and
            # return its parameters.
            binmoffparams = common.moffat_least_squares(
                axesdict["saxis"],
                bindata * datascale,
                headparams["seeing"],
                headparams["pixresolution"],
            )

            binmoffparams[0] /= datascale
            binmoffparams[4] /= datascale
            binmoffparams[5] /= datascale

            # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
            # previously fitted to it.
            LowExt, HighExt, fwhm, centre = common.extraction_limits(
                binmoffparams,
                axesdict,
                width_multiplier=params["-FWHM_MULTIPLIER"],
            )

            extractionlimits.append(
                [(bin[0] + bin[1]) * 0.5, LowExt, HighExt, centre]
            )

            # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
            # locations so they can be saved as metadata along with the extracted spectrum.
            binmoffparams.append(bin[0] + axesdict["wavstart"])
            binmoffparams.append(bin[1] + axesdict["wavstart"])
            extbin.append(binmoffparams)

            # DIAGNOSTICS - Plot computed moffat profile over data for each bin
            if params["-DIAG_PLOT_MOFFAT"]:
                common.plot_fitted_spatial_profile(
                    axesdict["saxis"],
                    bindata,
                    axesdict["hrsaxis"],
                    binmoffparams,
                    axesdict["imgstart"],
                    headparams,
                )

        binpars = np.array(extbin)
        sys.stdout.write("     Fitting complete.\n")

        sys.stdout.write(" >>> Drawing extraction aperture limits. ")
        sys.stdout.flush()
        extractionlimits = np.array(extractionlimits).T

        # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum. All pixels fully within the aperture are extracted.
        if params["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            drawlines = [
                extractionlimits[0] + axesdict["wavstart"],
                extractionlimits[1] + axesdict["imgstart"] - 1,
                extractionlimits[0] + axesdict["wavstart"],
                extractionlimits[2] + axesdict["imgstart"] + 1,
            ]

            common.show_img(
                framedict["data"],
                axesdict,
                headparams,
                drawlines,
                "2D Spectrum Overplotted with Extraction Limits",
            )

        # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
        # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.

        finalextractionlims = common.interpolate_extraction_lims(
            extractionlimits, axesdict["dispaxislen"]
        )
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis. All pixels fully within the aperture are extracted.
        if params["-DIAG_PLOT_EXTRACTION_LIMITS"]:
            drawlines = [
                np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
                finalextractionlims[0] + axesdict["imgstart"] - 1,
                np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
                finalextractionlims[1] + axesdict["imgstart"] + 1,
            ]

            common.show_img(
                framedict["data"],
                axesdict,
                headparams,
                drawlines,
                "2D Spectrum Overplotted with Full Extraction Limits",
            )

        # Extract the spectrum from a supersampled version of the 2D image using the extraction limits.
        sys.stdout.write(" >>> Extracting 1D spectrum. ")
        sys.stdout.flush()
        framedict["data"] = framedict["data"].T
        framedict["errs"] = framedict["errs"].T
        framedict["qual"] = framedict["qual"].T

        opdata1D, operrs1D, apdata1D, aperrs1D = common.optimal_extraction(
            framedict["data"],
            framedict["errs"],
            finalextractionlims,
            binpars,
            axesdict,
        )

        finalextractionlims = np.array(finalextractionlims)
        sys.stdout.write("DONE.\n")

        # DIAGNOSTICS - Plot extracted spectrum.
        if params["-PLOT_EXTRACTED_SPECTRUM"]:
            # DIAGNOSTICS, EXTRACTED SPECTRUM
            plt.figure(figsize=(9, 6))
            plt.errorbar(
                axesdict["waxis"],
                apdata1D,
                yerr=aperrs1D,
                color="k",
                marker=".",
                label="Aperture Spectrum",
            )
            plt.errorbar(
                axesdict["waxis"],
                opdata1D,
                yerr=operrs1D,
                color="r",
                marker=".",
                label="Optimal Spectrum",
            )
            plt.grid(alpha=0.5, linestyle="dotted")
            plt.title("Extracted 1D Spectrum")
            plt.ylabel("Flux, " + headparams["fluxunit"])
            plt.xlabel("Wavelength, " + headparams["wavunit"])
            plt.legend()
            plt.show()

        # Save the extracted spectrum to a new .fits file.
        if params["-SAVE"]:
            sys.stdout.write(" >>> Saving 1D spectrum and metadata.\n")
            sys.stdout.flush()
            save_fits(
                axesdict,
                headparams,
                opdata1D,
                operrs1D,
                apdata1D,
                aperrs1D,
                imghead,
                params,
                file_2D,
                moffparams,
                framedict,
                binpars,
                finalextractionlims,
                skybinpars,
                skyextlims,
            )

        sys.stdout.write(" >>> Extraction of " + file_2D + " completed.\n")

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
    """This function saves the extracted spectrum and intermediate products in a single FITS file.

    Description:
    Constructs a fits file.

    Args:
        axdict (dict): _description_
        hparams (_type_): _description_
        opflux (_type_): _description_
        operrs (_type_): _description_
        apflux (_type_): _description_
        aperrs (_type_): _description_
        head (_type_): _description_
        pars (_type_): _description_
        filename (_type_): _description_
        moffpars (_type_): _description_
        fdict (_type_): _description_
        bpars (_type_): _description_
        extractionlims (_type_): _description_
        sbpars (_type_): _description_
        skyextractionlims (_type_): _description_

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
        axdict["imgstart"],
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
        round(moffpars[1] + axdict["imgstart"], 5),
        "moffat profile centre",
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
        round(hparams["seeing"] * hparams["pixresolution"], 2),
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

    head["HIERARCH EXTRACTED HDU ROW 0"] = (
        "Wavelength Axis, " + hparams["wavunit"]
    )
    head.add_blank(
        "Data Saved in the Extracted Spectrum HDU",
        before="HIERARCH EXTRACTED HDU ROW 0",
    )
    head["HIERARCH EXTRACTED HDU ROW 1"] = "Flux, " + hparams["fluxunit"]
    head["HIERARCH EXTRACTED HDU ROW 2"] = (
        "Flux Uncertainty, " + hparams["fluxunit"]
    )
    head["EXTNAME"] = "OPTI_1D_SPEC"

    opfluxhdu = fits.PrimaryHDU([axdict["waxis"], opflux, operrs], header=head)
    apfluxhdu = fits.ImageHDU([axdict["waxis"], apflux, aperrs], header=head)
    apfluxhdu.header["EXTNAME"] = "APER_1D_SPEC"
    spec2Dhdu = fits.ImageHDU(fdict["ogdata"])
    spec2Dhdu.header["EXTNAME"] = "ORIG_2D_SPEC"
    errs2Dhdu = fits.ImageHDU(fdict["ogerrs"])
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
    hdulist.writeto("m" + "_".join(filenamelist[0:-1]) + "_" + filenamelist[-1])
    hdulist.close()

    sys.stdout.write(" >>> Spectrum extracted and saved:\n")
    sys.stdout.write(
        "     m" + "_".join(filenamelist[0:-1]) + "_" + filenamelist[-1] + "\n"
    )
    return None


def skyloc(framedict, axesdict, datascale, headparams, binparams, params):
    """Perform sky subtraction on the 2D spectrum.

    Description:
    For each dispersion bin determined by the get_bins function, median the bin along the dispersion axis, fit a
    moffat profile to the median data and then use the parameters of the fitted Moffat function to localise the 2D spectrum.

    Args:
        framedict (_type_): _description_
        axesdict (_type_): _description_
        datascale (_type_): _description_
        headparams (_type_): _description_
        binparams (_type_): _description_
        params (_type_): _description_

    Returns:
        _type_: _description_
    """

    sys.stdout.write(
        " >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n"
    )

    skybin = []
    extractionlimits = []
    for bin in binparams:
        # Take the median spatial profile of the dispersion
        # bin, and leave out pixel columns in the chip gaps if this is a GMOS spectrum.
        binimg = framedict["data"][:, bin[0] : bin[1]]
        chipgap = np.where(np.median(binimg, axis=0) != 1)
        bindata = np.nanmedian(binimg[:, chipgap[0]], axis=1)

        # Use a Levenberg-Marquardt Least Squares method to fit a Moffat function to the median spatial profile and
        # return its parameters.
        binmoffparams = common.moffat_least_squares(
            axesdict["saxis"],
            bindata * datascale,
            headparams["seeing"],
            headparams["pixresolution"],
        )

        # Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom,
        # so it's necessary to scale the data to a size that least squares can handle (see *1E18 above).
        # The shape of the profile fitted to the scaled spatial profile is the same as the unscaled,
        # but to to get a model profile that matches the original profile, the profile amplitude (A),
        # background level (B), and background gradient (m) all need to be scaled down again after the fitting.

        binmoffparams[0] /= datascale
        binmoffparams[4] /= datascale
        binmoffparams[5] /= datascale

        # Define the extraction limits of the current dispersion bin based on the parameters of the Moffat profile
        # previously fitted to it.
        LowExt, HighExt, fwhm, centre = common.extraction_limits(
            binmoffparams,
            axesdict,
            width_multiplier=params["-BG_FWHM_MULTIPLIER"],
        )

        extractionlimits.append(
            [(bin[0] + bin[1]) * 0.5, LowExt, HighExt, centre]
        )

        # Record the Moffat function parameters for each dispersion bin and add the wavstart offset to the bin
        # locations so they can be saved as metadata along with the extracted spectrum.
        binmoffparams.append(bin[0] + axesdict["wavstart"])
        binmoffparams.append(bin[1] + axesdict["wavstart"])
        skybin.append(binmoffparams)

        # DIAGNOSTICS - Plot computed moffat profile over data for each bin
        if params["-DIAG_PLOT_MOFFAT"]:
            common.plot_fitted_spatial_profile(
                axesdict["saxis"],
                bindata,
                axesdict["hrsaxis"],
                binmoffparams,
                axesdict["imgstart"],
                headparams,
            )

    skybin = np.array(skybin)

    sys.stdout.write("     Fitting complete.\n")

    sys.stdout.write(" >>> Drawing target/sky boundaries. ")
    sys.stdout.flush()
    extractionlimits = np.array(extractionlimits).T

    # DIAGNOSTICS - Plot the determined extraction limits over the 2D spectrum.
    if params["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        drawlines = [
            extractionlimits[0] + axesdict["wavstart"],
            (extractionlimits[1]) + axesdict["imgstart"],
            extractionlimits[0] + axesdict["wavstart"],
            (extractionlimits[2]) + axesdict["imgstart"],
        ]

        common.show_img(
            framedict["data"],
            axesdict,
            headparams,
            drawlines,
            "2D Spectrum Overplotted with Target/Sky Boundaries",
        )

    # Interpolate the extraction limits calculated for each median bin such that each wavelength element across the
    # entire unbinned wavelength axis of the entire 2D spectrum has its own extraction limits.

    skyextractionlims = common.interpolate_extraction_lims(
        extractionlimits, axesdict["dispaxislen"]
    )

    sys.stdout.write("DONE.\n")

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis.
    if params["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        drawlines = [
            np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
            (skyextractionlims[0]) + axesdict["imgstart"],
            np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
            (skyextractionlims[1]) + axesdict["imgstart"],
        ]

        common.show_img(
            framedict["data"],
            axesdict,
            headparams,
            drawlines,
            "2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )

    sys.stdout.write(" >>> Subtracting sky.\n")
    sys.stdout.flush()

    framedict = common.subtract_sky(
        skyextractionlims[0],
        skyextractionlims[1],
        framedict,
        axesdict,
        params,
        headparams,
    )

    sys.stdout.write("\n     DONE.\n")
    sys.stdout.flush()

    # DIAGNOSTICS - Plot the final extraction limits including the extrapolated sections at the ends of the wavelength axis.
    if params["-DIAG_PLOT_EXTRACTION_LIMITS"]:
        drawlines = [
            np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
            (skyextractionlims[0]) + axesdict["imgstart"],
            np.array(range(axesdict["dispaxislen"])) + axesdict["wavstart"],
            (skyextractionlims[1]) + axesdict["imgstart"],
        ]

        common.show_img(
            framedict["data"],
            axesdict,
            headparams,
            drawlines,
            "Sky Subtracted 2D Spectrum Overplotted with Full Target/Sky Boundaries",
        )

    return framedict, skybin, skyextractionlims


if __name__ == "__main__":
    motes()
