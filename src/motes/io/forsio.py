#!/usr/bin/env python3

"""
forsio.py - contains functions for reading and writing spectrum
            data files derived from FORS2 longslit observations.
"""

import copy
import motes.common as common
import numpy as np
import sys


def harvest_fors2(input_fits_hdu, primary_header):
    """
    Harvest the header and data from a FORS2 spectrum.

    Args:
        input_fits_hdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        primary_header (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        data (numpy.ndarray)   : the 2D data frame
        errs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5).
        qual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
                                    Since FORS2 spectra are not provided with a qual frame, a blank
                                    one (flagging all pixels as good; i.e. ==1) is created to
                                    ensure compatibility with MOTES.
        original_qual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
                                    In the case of FORS2 this frame is set to all zeros (see
                                    qual above).
        header_dict (dict)         : a dictionary containing the header information.
        wavelength_axis (numpy.ndarray)   : the 1D wavelength axis of the spectrum.
    """

    # Retrieve the data frame and error frame.
    data = input_fits_hdu[0].data
    errs = input_fits_hdu[1].data ** 0.5
    qual = np.ones(np.shape(data))
    original_qual = copy.deepcopy(qual) - 1

    # Determine the spatial pixel resolution of the image in arcsec depending on the binning of the
    # detector and the configuration of the collimator (high resolution or standard resolution).
    # If the pixel resolution can't be determined, complain and quit MOTES.
    if (
        primary_header["HIERARCH ESO DET WIN1 BINY"] == 1
        and primary_header["HIERARCH ESO INS COLL NAME"] == "COLL_SR"
    ):
        pixel_resolution = 0.125
    elif (
        primary_header["HIERARCH ESO DET WIN1 BINY"] == 1
        and primary_header["HIERARCH ESO INS COLL NAME"] == "COLL_HR"
    ):
        pixel_resolution = 0.0632
    elif (
        primary_header["HIERARCH ESO DET WIN1 BINY"] == 2
        and primary_header["HIERARCH ESO INS COLL NAME"] == "COLL_SR"
    ):
        pixel_resolution = 0.25
    elif (
        primary_header["HIERARCH ESO DET WIN1 BINY"] == 2
        and primary_header["HIERARCH ESO INS COLL NAME"] == "COLL_HR"
    ):
        pixel_resolution = 0.125
    else:
        sys.stdout.write("FAILED.\n")
        sys.stdout.write(
            "     Non-standard binning used in image.\n"
            "     Spatial pixel resolution could not be determined.\n"
        )
        sys.stdout.write("     Terminating MOTES.\n\n")
        sys.exit()

    sys.stdout.write(" >>> Spatial pixel resolution determined: " + str(pixel_resolution) + '"\n')

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "pixel_resolution": pixel_resolution,
        "exptime": primary_header["HIERARCH ESO INS SHUT EXPTIME"],
        "instrument": primary_header["INSTRUME"],
        "seeing": 0.5
        * (
            primary_header["HIERARCH ESO TEL AMBI FWHM START"]
            + primary_header["HIERARCH ESO TEL AMBI FWHM END"]
        ),
        "flux_unit": primary_header["BUNIT"],
        "wavelength_unit": "Angstroms",
    }

    sys.stdout.write("DONE.\n")

    # Create the wavelength axis of the spectrum.
    wavelength_axis = common.make_wav_axis(
        primary_header["CRVAL1"], primary_header["CD1_1"], primary_header["NAXIS1"]
    )

    return data, errs, qual, original_qual, header_dict, wavelength_axis


def save_fors(hdu_list, original_hdu_list):

    return hdu_list
