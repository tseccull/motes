#!/usr/bin/env python3

"""
floydsio.py - contains functions for reading and writing spectrum
              data files derived from FLOYDS longslit observations.
"""

import copy
import motes.common as common
import numpy as np
import sys


def harvest_floyds(input_fits_hdu, primary_header):
    """
    Harvest the header and data from a FLOYDS spectrum. Please note that this spectrum must not be
    flux calibrated, to ensure that a reliable ERR frame is made.

    Args:
        input_fits_hdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        primary_header (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        data (numpy.ndarray)   : the 2D data frame array
        errs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5). In the
                                    case of FLOYDS, no variance or uncertainty frame is provided,
                                    so one is constructed using the data along with read noise and
                                    dark current metadata contained in the file header. This is why
                                    flux calibrated 2D FLOYDS spectra should not be extracted with
                                    MOTES, as the flux calibration spoils the construction of the
                                    error frame. Flux calibration should be applied after the
                                    spectrum is extracted if MOTES is used.
        qual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
                                    Since FLOYDS spectra are not provided with a qual frame, a
                                    blank one (flagging all pixels as good; i.e. ==1) is created to
                                    ensure compatibility with MOTES.
        original_qual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
                                    In the case of FLOYDS this frame is set to all ones (see
                                    above).
        header_dict (dict)         : a dictionary containing the header information required by
                                    MOTES.
    """

    # Retrieve the HDU and extract/construct the 2D data, err, and qual frames.
    data = input_fits_hdu[0].data
    errs = np.sqrt(
        data
        + (primary_header["RDNOISE"] * primary_header["RDNOISE"])
        + (primary_header["DARKCURR"] * primary_header["DARKCURR"])
    )
    qual = np.ones(np.shape(data))
    original_qual = copy.deepcopy(qual) - 1

    # Determine the spatial pixel resolution of the image in arcsec.
    pixel_resolution = primary_header["PIXSCALE"]
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: " + str(pixel_resolution).strip("0") + '"\n'
    )

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "pixel_resolution": pixel_resolution,
        "exptime": primary_header["EXPTIME"],
        "instrument": primary_header["INSTRUME"],
        "seeing": primary_header["AGFWHM"],  # Grabs estimated FWHM from autoguider.
        "flux_unit": "electrons",
        "wavelength_unit": primary_header["WAT2_001"].split(" ")[2].split("=")[1],
    }
    sys.stdout.write("DONE.\n")

    # Create the wavelength axis of the spectrum.
    wavelength_axis = common.make_wav_axis(
        primary_header["CRVAL1"], primary_header["CD1_1"], primary_header["NAXIS1"]
    )

    return data, errs, qual, original_qual, header_dict, wavelength_axis


def save_floyds(hdu_list, original_hdu_list):

    return hdu_list
