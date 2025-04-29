#!/usr/bin/env python3

"""
	gmosio.py

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
	gmosio.py contains functions for reading 2D spectrum data files
	derived from GMOS longslit observations and writing 1D spectrum
    data files produced by MOTES for the same datasets.
"""

import logging
import numpy as np

logger = logging.getLogger("motes")


def harvest_gmos(input_fits_hdu):
    """
    Harvest the header and data from a GMOS spectrum.

    Args:
     -- input_fits_hdu (astropy.io.fits.hdu.hdulist.HDUList)
          The list of Header Data Units (HDUs) read in from the data
          file.

    Returns:
     -- data (numpy.ndarray)
          The 2D data frame
     -- errs (numpy.ndarray)
          The 2D error/uncertainty frame (variance_frame^0.5).
     -- qual (numpy.ndarray)
          The 2D quality frame noting the locations of bad pixels etc.
     -- original_qual (numpy.ndarray)
          The original 2D quality frame prior to manipulation by MOTES.
     -- header_dict (dict)
          A dictionary containing the header information.
     -- wavelength_axis (numpy.ndarray)
          The 1D wavelength axis of the spectrum.
    """

    # Retrieve the data frame, error frame, and qual frame. Also
    # retrieve the header of the science image frame, as some metadata
    # is stored there instead of the primary header.
    primary_header = input_fits_hdu[0].header
    science_header = input_fits_hdu["SCI"].header
    data = input_fits_hdu["SCI"].data
    errs = input_fits_hdu["VAR"].data ** 0.5
    qual = input_fits_hdu["DQ"].data

	# Sets all data and errs within the GMOS chip gaps to 1, so they
	# don't get flagged as bad pixels or trip up the bin definition
	# stage.
	# Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    data[qual==16] = 1.
    errs[qual==16] = 1.
    new_qual = np.ones(np.shape(qual))
    new_qual[qual>0] = 0
    qual = new_qual

    # Make the wavelength axis.
    if science_header["CRPIX1"] > 0.5:
        reference_pixel = np.ceil(science_header["CRPIX1"])
        reference_shift = reference_pixel-science_header["CRPIX1"]
    else:
        reference_pixel = np.floor(science_header["CRPIX1"])
        reference_shift = science_header["CRPIX1"] - reference_pixel
    
    central_wavelength = science_header["CRVAL1"]
    delta_wavelength = science_header["CD1_1"]
    reference_wavelength = (reference_shift * delta_wavelength) + central_wavelength
    wavelen_pixel_axis = np.arange(-reference_pixel+1, science_header["NAXIS1"]-reference_pixel)
    wavelength_axis = (wavelen_pixel_axis * science_header["CD1_1"]) + reference_wavelength
    
    # If DRAGONS has been used to reduce the data, flip the wavelength
    # axis and the 2D data.
    #if science_header["CD1_1"] < 0:
    #    wavelength_axis = np.flip(wavelength_axis)
    #    data = np.flip(data, axis=1)
    #    errs = np.flip(errs, axis=1)
    #    qual = np.flip(qual, axis=1)

    # All this is to get an initial estimate of the IQ. Tables below are
    # based on the condition constraints used by Gemini.
    # See https://www.gemini.edu/observing/telescopes-and-sites/sites#ImageQuality
    iq_dict = {
        "20-percentile": 0,
        "70-percentile": 1,
        "85-percentile": 2,
        "100-percentile": 3,
        "Any": 3,
        "UNKNOWN": 3,
    }

    wavelength_table = np.array(
        [
            [000.0, 400.0, 0],
            [400.0, 550.0, 1],
            [550.0, 700.0, 2],
            [700.0, 850.0, 3],
            [850.0, 975.0, 4],
            [975.0, 1100.0, 5],
        ]
    )

    iq_table = np.array(
        [
            [0.6, 0.90, 1.20, 2.00],
            [0.6, 0.85, 1.10, 1.90],
            [0.5, 0.75, 1.05, 1.80],
            [0.5, 0.75, 1.05, 1.70],
            [0.5, 0.70, 0.95, 1.70],
            [0.4, 0.70, 0.95, 1.65],
        ]
    )

    iq = primary_header["RAWIQ"]
    
    short_wavelength = np.min(wavelength_axis)
    for i in wavelength_table:
        if short_wavelength > i[0] and short_wavelength < i[1]:
            seeing = float(iq_table[int(i[2])][int(iq_dict[iq])])
            break

    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "exptime": primary_header["EXPTIME"],
        "seeing": seeing,
        "instrument": primary_header["INSTRUME"],
        "flux_unit": science_header["BUNIT"],
        "wavelength_unit": science_header["CUNIT1"],
        "pixel_resolution": float(science_header["PIXSCALE"])
    }

    return data, errs, qual, header_dict, wavelength_axis


def save_gmos(hdu_list, original_hdu_list):
    """
    Prepare the original HDUs for inclusion into the output save file.
    
    Args:
     -- hdu_list (list)
          List of HDUs for inclusion in the output save file.
     -- original_hdu_list (class astropy.io.fits.hdu.hdulist.HDUList)
          HDUList read in from the original input file.
            
    Return:
     -- hdu_list (list)
          Same as hdu_list arg, but with relevant original HDUs added
          from original_hdu_list.
    """

    original_hdu_list["MDF"].header["EXTNAME"] = "ORIG_MDF"
    hdu_list.append(original_hdu_list["ORIG_MDF"])
    original_hdu_list["SCI"].header["EXTNAME"] = "ORIG_SCI"
    hdu_list.append(original_hdu_list["ORIG_SCI"])
    original_hdu_list["VAR"].header["EXTNAME"] = "ORIG_VAR"
    hdu_list.append(original_hdu_list["ORIG_VAR"])
    original_hdu_list["DQ"].header["EXTNAME"] = "ORIG_DQ"
    hdu_list.append(original_hdu_list["ORIG_DQ"])
    
    original_sci_keys = ["PIXSCALE", "CCDSUM", "GAIN", "GAINSET", "RDNOISE","CUNIT1"]
    original_sci_comments = [
        "Pixel scale in Y, ''/pixel",
        "Detector binning, pixels",
        "Amplifier gain",
        "Gain setting (low/high)",
        "Readout noise",
        "Wavelength units"]
    
    for i, key in enumerate(original_sci_keys):
        hdu_list[0].header[key] = (
            original_hdu_list["ORIG_SCI"].header[key],
            original_sci_comments[i]
        )
       
    return hdu_list
