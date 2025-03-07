#!/usr/bin/env python3

"""
xshooio.py - contains functions for reading and writing spectrum data
             files derived from X-Shooter longslit observations.
"""

import copy
import motes.common as common
import numpy as np
import sys


def harvest_xshoo(input_fits_hdu, primary_header):
    """
    Harvest the header and data from an X-Shooter spectrum.

    Args:
        input_fits_hdu (astropy.io.fits.hdu.image.PrimaryHDU) : the Header Data Unit (HDU) read in from
                                                            the data file.
        primary_header (astropy.io.fits.header.Header)         : the header read in from the data file.

    Returns:
        data (numpy.ndarray)   : the 2D data frame
        errs (numpy.ndarray)   : the 2D error/uncertainty frame (variance_frame^0.5).
        qual (numpy.ndarray)   : the 2D quality frame noting the locations of bad pixels etc.
        original_qual (numpy.ndarray) : the original 2D quality frame prior to manipulation by MOTES.
        header_dict (dict)         : a dictionary containing the header information.
        wavelength_axis (numpy.ndarray)   : the 1D wavelength axis of the spectrum.
    """

    print(type(input_fits_hdu))
    print(type(primary_header))

    # Retrieve the data frame, error frame, and qual frame.
    data = input_fits_hdu[0].data
    errs = input_fits_hdu[1].data
    qual = input_fits_hdu[2].data
    original_qual = copy.deepcopy(qual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    # Values flagged by the X-Shooter data reduction pipeline as interpolated are considered good.
    qual[qual == 4194304] = 0

    # If the bspline sky subtraction method has been used in the X-Shooter data reduction pipeline,
    # pixels flagged as outliers or inaccurate are considered good.
    qual[qual == 8388608] = 0
    qual[qual == 16777216] = 0
    qual[qual > 0] *= -1
    qual[qual == 0] = 1
    qual[qual < 0] = 0
    qual[~np.isfinite(data)] = 0

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "pixel_resolution": primary_header["CDELT2"],
        "exptime": primary_header["EXPTIME"],
        "seeing": 0.5
        * (
            primary_header["HIERARCH ESO TEL AMBI FWHM START"]
            + primary_header["HIERARCH ESO TEL AMBI FWHM END"]
        ),
        "flux_unit": primary_header["BUNIT"],
        "wavelength_unit": "nm",
    }
    sys.stdout.write("DONE.\n")
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(header_dict["pixel_resolution"])
        + '"\n'
    )

    # Create the wavelength axis of the spectrum.
    wavelength_axis = common.make_wav_axis(
        primary_header["CRVAL1"], primary_header["CDELT1"], primary_header["NAXIS1"]
    )

    return data, errs, qual, original_qual, header_dict, wavelength_axis


def save_xshoo(hdu_list, original_hdu_list):

    return hdu_list
