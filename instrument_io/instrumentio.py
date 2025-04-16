#!/usr/bin/env python3

"""
	[instrumentio].py

	Copyright (C) [20XX Author Name(s)]
	
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

	Last updated - [YYYY-MM-DD]

	Description---------------------------------------------------------
	[instrumentio].py contains functions for reading 2D spectrum data files
	derived from [Instrument Name] longlist observations and writing 1D spectrum
    data files produced by MOTES for the same datasets.
"""

import logging
import numpy as np

logger = logging.getLogger("motes")


def harvest_instrument(input_fits_hdu):
    """
    Harvest the header and data from a[n Instrument Name] spectrum.

    Args:
     -- input_fits_hdu (astropy.io.fits.hdu.image.PrimaryHDU)
          The Header Data Unit (HDU) read in from the data file.
     -- primary_header (astropy.io.fits.header.Header)
          The header read in from the data file.

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

    # Retrieve the required header(s).
    header_one = input_fits_hdu[0].header
    optional_header = input_fits_hdu["HDU_EXTNAME"].header
    
    # Retrieve the 2D science data frame. MOTES assumes that the 
    # columns of the data frame align with the spectrum's spatial axis
    # and that its rows align with the dispersion axis.
    # If the spectrum does not need to be oriented, then use the
    # following.
    data = input_fits_hdu["DATA_EXTNAME"].data
    # If the data needs to be transposed, 
    # data = input_fits_hdu["DATA_EXTNAME"].data.T
    # If the data needs to be transposed, be careful to ensure that
    # the correct header keys are used for other steps in this
    # function (e.g. making the wavelength axis.
    
    # Retrieve the 2D uncertainty data frame. Is this frame is a 
    # variance frame, square root it to make sure errs represents the
    # uncertainty of the science data.
    errs = input_fits_hdu["UNCERTAINTY_EXTNAME"].data
    
    # If the input file has a data quality frame, retrieve it. MOTES
    # only accepts a boolean quality frame (0=GOOD; 1=BAD). If the
    # qual frame flags for multiple different pixel statuses with
    # different values, set those which can be processed to 0 and
    # those which should be avoided to 1. For example, the pixels
    # falling within chip gaps in GMOS spectra are flagged with the
    # value 16 in the associated quality frame. We want to ensure
    # MOTES skips them when localising the spectrum, so we set them
    # to a value of 1. 
    qual = input_fits_hdu["QUALITY_EXTNAME"].data
    data[qual==16] = 1.
    errs[qual==16] = 1.
    new_qual = np.ones(np.shape(qual))
    new_qual[qual>0] = 0
    qual = new_qual
    
    # If your input spectrum has no quality frame, make one now by
    # creating an array of zeros that is the same shape as the data
    # and uncertainty frames.
    # qual = np.zeros(np.shape(data))


    # Make the wavelength axis based on information contained in the
    # header(s). Often, the same header cards are used in different
    # instruments.
    # Reference (start or central) wavelength = header["CRVAL1"]
    # Reference (start or central) pixel index = header["CRPIX1"]
    # Wavelength increment per pixel = header["CD1_1"]
    # If the reference pixel is at the start of the wavelength axis,
    # the axis can be created quite simply. Be mindful of the
    # direction in which the wavelength increases when making the
    # wavelength axis
    wavelength_axis = (np.arange(header["NAXIS1"])*header["CD1_1"]) + header["CRVAL1"]
    
    # If the reference wavelength is the central wavelength, the
    # implementation is a little tricker. See the gmosio.py harvester
    # function as an example.

    # Next a dictionary is created to hold values taken or calculated
    # from the header(s). MOTES needs a number of specific values:
    # OBJECT NAME
    # FRAME EXPOSURE TIME
    # SEEING/IQ ESTIMATE.
    #   This might be provided in a keyword from the observatory, or
    #   it might be calculated based on the values of other keywords
    #   in the the header(s). See the gmosio.py harvester function for
    #   an example of something that's a bit more complex.
    # INSTRUMENT NAME.
    # BRIGHTNESS UNIT. 
    #   Unit used to measure brightness in the image (this could be
    #   electrons, ADUs, or something more elaborate if a flux
    #   calibration has been done).
    # WAVELENGTH UNIT.
    # PIXEL RESOLUTION / PLATE SCALE / PIXEL SCALE.
    #   This is the scale of the spatial axis of the spectrum in
    #   arcseconds/pixel.
    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "exptime": primary_header["EXPTIME"],
        "seeing": primary_header["FWHM"],
        "instrument": primary_header["INSTRUME"],
        "flux_unit": science_header["BUNIT"],
        "wavelength_unit": science_header["CUNIT1"],
        "pixel_resolution": float(science_header["PIXSCALE"])
    }

    return data, errs, qual, header_dict, wavelength_axis


def save_instrument(hdu_list, original_hdu_list):
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

    # Add any original unaltered HDUs from the input file to the
    # HDU list of the output file to store them as meta data frames.
    # Remember to set appropriate names for each HDU by updating
    # the EXTNAME card in each header.
    original_hdu_list["MDF"].header["EXTNAME"] = "ORIG_MDF"
    hdu_list.append(original_hdu_list["ORIG_MDF"])
    original_hdu_list["SCI"].header["EXTNAME"] = "ORIG_SCI"
    hdu_list.append(original_hdu_list["ORIG_SCI"])
    original_hdu_list["VAR"].header["EXTNAME"] = "ORIG_VAR"
    hdu_list.append(original_hdu_list["ORIG_VAR"])
    original_hdu_list["DQ"].header["EXTNAME"] = "ORIG_DQ"
    hdu_list.append(original_hdu_list["ORIG_DQ"])
    
    # Retrieve and store the values of any pertinent header cards
    # in the primary header of the output file. You may want to do
    # this to collate all necessary header keywords in the same header
    # if they were initially spread across the headers of more than
    # one HDU in the input file.
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
