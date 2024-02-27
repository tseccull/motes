"""
harvester.py - A collection of functions to read in data from the image file and repackage it into
               a dictionary for use in the rest of MOTES.
"""

import copy
import sys

import astropy.io.fits as fits
import numpy as np

from motes import common


def data_harvest(reg_counter, filename_2D, region):
    """Extract header metadata and data frames, and repackage them into dictionaries for use in the
       rest of MOTES.

    Args:
        reg_counter (int) : An integer noting which line in region is being called to define the
                            boundaries of the 2D data.
        filename_2D (str) : The name of the data file.
        region (list)     : A list of regions read in from reg.txt in startup.read_regions()

    Returns:
        head_dict (dict)  : A dictionary containing parameters and metadata read from the header of
                            the image file.
        frame_dict (dict) : A dictionary containing the 2D data frames read from the image file.
        axes_dict (dict)  : A dictionary containing the spatial and spectral axis arrays associated
                            with the data frames, along with metadata used to define the boundaries
                            of the 2D data.
        imghead (dict)    : A copy of the image file header; this is also a dictionary.
    """

    # Create dictionary to tell data_harvest which instrument specific function to call.
    instrument_dict = {
        "en06": harvest_floyds,
        "en12": harvest_floyds,
        "FORS2": harvest_fors2,
        "GMOS-N": harvest_gmos,
        "GMOS-S": harvest_gmos,
        "XSHOOTER": harvest_xshoo,
    }

    # Open file containing the spectral data, then extract the header, image frame, error frame,
    # and quality frame (if the file has one).
    with fits.open(filename_2D) as imgfile:
        imghead = imgfile[0].header
        inst = imghead["INSTRUME"]
        # Based on the value of inst, this calls one of the harvest_instrument functions.
        (
            imgdata,
            imgerrs,
            imgqual,
            ogimgqual,
            head_dict,
            wavaxis,
        ) = instrument_dict[
            inst
        ](imgfile, imghead)

    # Slice all dataframes based on the input from reg.txt
    imgshape = np.shape(imgdata)
    imgstart = int(0 + region[reg_counter][0])
    imgend = int(imgshape[0] - region[reg_counter][1])

    # Slice off the spatial rows outside the spatial region.
    datasliced = imgdata[imgstart : imgend + 1, :]
    errssliced = imgerrs[imgstart : imgend + 1, :]
    qualsliced = imgqual[imgstart : imgend + 1, :]
    dataslicedshape = np.shape(datasliced)
    sys.stdout.write(
        " >>> 2D spectrum sliced on spatial axis based on user defined limits:\n"
        "     New spatial axis covers pixel rows "
        + str(imgstart)
        + "-"
        + str(imgend)
        + ".\n"
    )

    # Create spatial axis for the 2D spectrum and a high resolution version (standard res * 5) for
    # the purposes of plotting
    spataxis = np.linspace(0.0, float(dataslicedshape[0] - 1), num=dataslicedshape[0])
    hiresspataxis = np.linspace(spataxis[0], spataxis[-1], num=len(spataxis) * 5)

    if region[reg_counter][2] < wavaxis[0] or region[reg_counter][3] > wavaxis[-1]:
        sys.stdout.write(
            " >>> User defined wavelength limit(s) are outside native wavelength range\n"
            '     Make sure "-LOW_WAV_SLICE" > lower limit of wavelength axis\n'
            '     Make sure "-HIGH_WAV_SLICE" < upper limit of wavelength axis\n'
            "     Terminating MOTES.\n\n"
        )
        sys.exit()

    wavslice = np.where(
        np.logical_and(
            wavaxis >= region[reg_counter][2], wavaxis <= region[reg_counter][3]
        )
    )
    wavstart = wavslice[0][0]
    wavend = wavslice[0][-1]
    wavaxis = wavaxis[wavslice]

    datasliced = np.squeeze(datasliced[:, wavslice])
    errssliced = np.squeeze(errssliced[:, wavslice])
    qualsliced = np.squeeze(qualsliced[:, wavslice])

    sys.stdout.write(
        " >>> 2D spectrum sliced on dispersion axis based on user defined limits:\n"
        "     New Wavelength range is "
        + str(region[reg_counter][2])
        + "-"
        + str(region[reg_counter][3])
        + " "
        + head_dict["wavunit"]
        + ".\n"
        "     This range is equivalent to pixel columns "
        + str(wavstart)
        + "-"
        + str(wavstart + len(wavaxis))
        + "\n"
    )

    frame_dict = {
        "data": datasliced,
        "errs": errssliced,
        "qual": qualsliced,
        "ogqual": ogimgqual,
    }

    axes_dict = {
        "spataxislen": len(spataxis),
        "spatial_axis": spataxis,
        "hi_resolution_spatial_axis": hiresspataxis,
        "data_spatial_floor": imgstart,
        "imgend": imgend,
        "dispersion_axis_length": len(wavaxis),
        "wavelength_axis": wavaxis,
        "wavelength_start": wavstart,
        "wavend": wavend,
    }

    return head_dict, frame_dict, axes_dict, imghead


def harvest_floyds(imgfilehdu, imgheader):
    """
    Harvest the header and data from a FLOYDS spectrum. Please note that this spectrum must not be
    flux calibrated, to ensure that a reliable ERR frame is made.

    Args:
        imgfilehdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        imgheader (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        imgdata (numpy.ndarray)   : the 2D data frame array
        imgerrs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5). In the
                                    case of FLOYDS, no variance or uncertainty frame is provided,
                                    so one is constructed using the data along with read noise and
                                    dark current metadata contained in the file header. This is why
                                    flux calibrated 2D FLOYDS spectra should not be extracted with
                                    MOTES, as the flux calibration spoils the construction of the
                                    error frame. Flux calibration should be applied after the
                                    spectrum is extracted if MOTES is used.
        imgqual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
                                    Since FLOYDS spectra are not provided with a qual frame, a
                                    blank one (flagging all pixels as good; i.e. ==1) is created to
                                    ensure compatibility with MOTES.
        ogimgqual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
                                    In the case of FLOYDS this frame is set to all ones (see
                                    above).
        headerdict (dict)         : a dictionary containing the header information required by
                                    MOTES.
    """

    # Retrieve the HDU and extract/construct the 2D data, err, and qual frames.
    imgdata = imgfilehdu[0].data
    imgerrs = np.sqrt(
        imgdata
        + (imgheader["RDNOISE"] * imgheader["RDNOISE"])
        + (imgheader["DARKCURR"] * imgheader["DARKCURR"])
    )
    imgqual = np.ones(np.shape(imgdata))
    ogimgqual = copy.deepcopy(imgqual) - 1

    # Determine the spatial pixel resolution of the image in arcsec.
    pixres = imgheader["PIXSCALE"]
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: " + str(pixres).strip("0") + '"\n'
    )

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    headerdict = {
        "object": imgheader["OBJECT"].replace(" ", "_"),
        "pixel_resolution": pixres,
        "exptime": imgheader["EXPTIME"],
        "inst": imgheader["INSTRUME"],
        "seeing": imgheader["AGFWHM"],  # Grabs estimated FWHM from autoguider.
        "flux_unit": "electrons",
        "wavunit": imgheader["WAT2_001"].split(" ")[2].split("=")[1],
    }
    sys.stdout.write("DONE.\n")

    # Create the wavelength axis of the spectrum.
    wavaxis = common.make_wav_axis(
        imgheader["CRVAL1"], imgheader["CD1_1"], imgheader["NAXIS1"]
    )

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis


def harvest_fors2(imgfilehdu, imgheader):
    """
    Harvest the header and data from a FORS2 spectrum.

    Args:
        imgfilehdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        imgheader (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        imgdata (numpy.ndarray)   : the 2D data frame
        imgerrs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5).
        imgqual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
                                    Since FORS2 spectra are not provided with a qual frame, a blank
                                    one (flagging all pixels as good; i.e. ==1) is created to
                                    ensure compatibility with MOTES.
        ogimgqual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
                                    In the case of FORS2 this frame is set to all zeros (see
                                    imgqual above).
        headerdict (dict)         : a dictionary containing the header information.
        wavaxis (numpy.ndarray)   : the 1D wavelength axis of the spectrum.
    """

    # Retrieve the data frame and error frame.
    imgdata = imgfilehdu[0].data
    imgerrs = imgfilehdu[1].data ** 0.5
    imgqual = np.ones(np.shape(imgdata))
    ogimgqual = copy.deepcopy(imgqual) - 1

    # Determine the spatial pixel resolution of the image in arcsec depending on the binning of the
    # detector and the configuration of the collimator (high resolution or standard resolution).
    # If the pixel resolution can't be determined, complain and quit MOTES.
    if (
        imgheader["HIERARCH ESO DET WIN1 BINY"] == 1
        and imgheader["HIERARCH ESO INS COLL NAME"] == "COLL_SR"
    ):
        pixres = 0.125
    elif (
        imgheader["HIERARCH ESO DET WIN1 BINY"] == 1
        and imgheader["HIERARCH ESO INS COLL NAME"] == "COLL_HR"
    ):
        pixres = 0.0632
    elif (
        imgheader["HIERARCH ESO DET WIN1 BINY"] == 2
        and imgheader["HIERARCH ESO INS COLL NAME"] == "COLL_SR"
    ):
        pixres = 0.25
    elif (
        imgheader["HIERARCH ESO DET WIN1 BINY"] == 2
        and imgheader["HIERARCH ESO INS COLL NAME"] == "COLL_HR"
    ):
        pixres = 0.125
    else:
        sys.stdout.write("FAILED.\n")
        sys.stdout.write(
            "     Non-standard binning used in image.\n"
            "     Spatial pixel resolution could not be determined.\n"
        )
        sys.stdout.write("     Terminating MOTES.\n\n")
        sys.exit()

    sys.stdout.write(" >>> Spatial pixel resolution determined: " + str(pixres) + '"\n')

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    headerdict = {
        "object": imgheader["OBJECT"].replace(" ", "_"),
        "pixel_resolution": pixres,
        "exptime": imgheader["HIERARCH ESO INS SHUT EXPTIME"],
        "inst": imgheader["INSTRUME"],
        "seeing": 0.5
        * (
            imgheader["HIERARCH ESO TEL AMBI FWHM START"]
            + imgheader["HIERARCH ESO TEL AMBI FWHM END"]
        ),
        "flux_unit": imgheader["BUNIT"],
        "wavunit": "Angstroms",
    }

    sys.stdout.write("DONE.\n")

    # Create the wavelength axis of the spectrum.
    wavaxis = common.make_wav_axis(
        imgheader["CRVAL1"], imgheader["CD1_1"], imgheader["NAXIS1"]
    )

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis


def harvest_gmos(imgfilehdu, imgheader):
    """
    Harvest the header and data from a GMOS spectrum.

    Args:
        imgfilehdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        imgheader (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        imgdata (numpy.ndarray)   : the 2D data frame
        imgerrs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5).
        imgqual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
        ogimgqual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
        headerdict (dict)         : a dictionary containing the header information.
        wavaxis (numpy.ndarray)   : the 1D wavelength axis of the spectrum.
    """

    # Retrieve the data frame, error frame, and qual frame. Also retrieve the header of the science
    # image frame, as some metadata is stored there instead of the primary header.
    scihead = imgfilehdu["SCI"].header
    imgdata = imgfilehdu["SCI"].data
    imgerrs = imgfilehdu["VAR"].data ** 0.5
    imgqual = imgfilehdu["DQ"].data
    ogimgqual = copy.deepcopy(imgqual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    # Pixels in chip gaps are kept as 1 to make sure they don't get flagged as bad.
    imgqual[np.where(imgdata + imgqual == 1)] = 0
    imgqual = 1 - imgqual

    # Set pixels with NaN value to 1 in the data frame, and flag them as bad pixels in the qual
    # frame.
    imgqual[~np.isfinite(imgdata)] = 0
    imgdata[~np.isfinite(imgdata)] = 1.0

    # All this is to get an initial estimate of the IQ. Tables below are based on the condition
    # constraints used by Gemini.
    # See https://www.gemini.edu/observing/telescopes-and-sites/sites#ImageQuality
    IQ_dict = {
        "20-percentile": 0,
        "70-percentile": 1,
        "85-percentile": 2,
        "100-percentile": 3,
        "Any": 3,
        "UNKNOWN": 3,
    }

    WavTab = np.array(
        [
            [0000.0, 4000.0, 0],
            [4000.0, 5500.0, 1],
            [5500.0, 7000.0, 2],
            [7000.0, 8500.0, 3],
            [8500.0, 9750.0, 4],
            [9750.0, 11000.0, 5],
        ]
    )

    IQTab = np.array(
        [
            [0.6, 0.90, 1.20, 2.00],
            [0.6, 0.85, 1.10, 1.90],
            [0.5, 0.75, 1.05, 1.80],
            [0.5, 0.75, 1.05, 1.70],
            [0.5, 0.70, 0.95, 1.70],
            [0.4, 0.70, 0.95, 1.65],
        ]
    )

    iq = imgheader["RAWIQ"]

    for i in WavTab:
        if scihead["CRVAL1"] > i[0] and scihead["CRVAL1"] < i[1]:
            seeing = float(IQTab[int(i[2])][int(IQ_dict[iq])])
            break

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    headerdict = {
        "object": imgheader["OBJECT"].replace(" ", "_"),
        "pixel_resolution": float(imgheader["PIXSCALE"]),
        "exptime": imgheader["EXPTIME"],
        "seeing": seeing,
        "inst": imgheader["INSTRUME"],
        "wavunit": scihead["WAT1_001"].split(" ")[2].split("=")[1],
    }

    # BUNIT only appears in the headers of GMOS spectra if they have been flux calibrated.
    if "BUNIT" in scihead:
        headerdict["flux_unit"] = scihead["BUNIT"]
    else:
        headerdict["flux_unit"] = "electrons"

    sys.stdout.write("DONE.\n")
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(headerdict["pixel_resolution"])
        + '"\n'
    )

    # Create the wavelength axis of the spectrum.
    wavaxis = common.make_wav_axis(
        scihead["CRVAL1"], scihead["CD1_1"], scihead["NAXIS1"]
    )

    # Sets all data and errs within the GMOS chip gaps to 1, so they don't get flagged as bad
    # pixels or trip up the bin definition stage. Chip gaps are identified as pixel columns which
    # are all zeros, and then three columns either side of the chip gaps are also flagged just to
    # be safe.
    zerorows = [1 if all(x == 0) else 0 for x in imgdata.T]
    boundary_cols = 3
    zerorows = np.concatenate(
        [np.zeros(boundary_cols), zerorows, np.zeros(boundary_cols)]
    )
    for n in reversed(range(boundary_cols)):
        zerorows = [
            1 if x == 0 and zerorows[y + 1] == 1 else x
            for y, x in enumerate(zerorows[:-1])
        ]
        zerorows = [
            1 if x == 0 and zerorows[y - 1] == 1 else x
            for y, x in enumerate(zerorows[1:])
        ]
    zerorows = [1 if x > 0 else x for x in zerorows]
    chipgapmap = np.tile(zerorows, (np.shape(imgdata)[0], 1))
    imgdata[chipgapmap == 1] = 1.0
    imgerrs[chipgapmap == 1] = 1.0

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis


def harvest_xshoo(imgfilehdu, imgheader):
    """
    Harvest the header and data from an X-Shooter spectrum.

    Args:
        imgfilehdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        imgheader (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        imgdata (numpy.ndarray)   : the 2D data frame
        imgerrs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5).
        imgqual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
        ogimgqual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
        headerdict (dict)         : a dictionary containing the header information.
        wavaxis (numpy.ndarray)   : the 1D wavelength axis of the spectrum.
    """

    print(type(imgfilehdu))
    print(type(imgheader))

    # Retrieve the data frame, error frame, and qual frame.
    imgdata = imgfilehdu[0].data
    imgerrs = imgfilehdu[1].data
    imgqual = imgfilehdu[2].data
    ogimgqual = copy.deepcopy(imgqual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    # Values flagged by the X-Shooter data reduction pipeline as interpolated are considered good.
    imgqual[imgqual == 4194304] = 0

    # If the bspline sky subtraction method has been used in the X-Shooter data reduction pipeline,
    # pixels flagged as outliers or inaccurate are considered good.
    imgqual[imgqual == 8388608] = 0
    imgqual[imgqual == 16777216] = 0
    imgqual[imgqual > 0] *= -1
    imgqual[imgqual == 0] = 1
    imgqual[imgqual < 0] = 0
    imgqual[~np.isfinite(imgdata)] = 0

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    headerdict = {
        "object": imgheader["OBJECT"].replace(" ", "_"),
        "pixel_resolution": imgheader["CDELT2"],
        "exptime": imgheader["EXPTIME"],
        "seeing": 0.5
        * (
            imgheader["HIERARCH ESO TEL AMBI FWHM START"]
            + imgheader["HIERARCH ESO TEL AMBI FWHM END"]
        ),
        "flux_unit": imgheader["BUNIT"],
        "wavunit": "nm",
    }
    sys.stdout.write("DONE.\n")
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(headerdict["pixel_resolution"])
        + '"\n'
    )

    # Create the wavelength axis of the spectrum.
    wavaxis = common.make_wav_axis(
        imgheader["CRVAL1"], imgheader["CDELT1"], imgheader["NAXIS1"]
    )

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis
