#!/home/tom/miniforge3/envs/work/bin/python

import numpy as np
import sys


def harvest_gmos(input_fits_hdu):
    """
    Harvest the header and data from a GMOS spectrum.

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
            [0000.0, 4000.0, 0],
            [4000.0, 5500.0, 1],
            [5500.0, 7000.0, 2],
            [7000.0, 8500.0, 3],
            [8500.0, 9750.0, 4],
            [9750.0, 11000.0, 5],
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

    for i in wavelength_table:
        if science_header["CRVAL1"] > i[0] and science_header["CRVAL1"] < i[1]:
            seeing = float(iq_table[int(i[2])][int(iq_dict[iq])])
            break

    # Put header information into a dictionary
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()
    header_dict = {
        "object": primary_header["OBJECT"].replace(" ", "_"),
        "exptime": primary_header["EXPTIME"],
        "seeing": seeing,
        "instrument": primary_header["INSTRUME"]
    }

    # BUNIT only appears in the headers of GMOS spectra if they have
    # been flux calibrated of if they have been reduced with DRAGONS.
    if "BUNIT" in science_header:
        header_dict["flux_unit"] = science_header["BUNIT"]
    else:
        header_dict["flux_unit"] = "electron"
        
    if "CUNIT1" in science_header:
        header_dict["wavelength_unit"] = science_header["CUNIT1"]
    else:
        header_dict["wavelenght_unit"] = (
            science_header["WAT1_001"].split(" ")[2].split("=")[1]
        )
    
    if "PIXSCALE" in science_header:
        header_dict["pixel_resolution"] = float(science_header["PIXSCALE"])
    else:
        header_dict["pixel_resolution"] = float(primary_header["PIXSCALE"])

    sys.stdout.write("DONE.\n")
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(header_dict["pixel_resolution"])
        + '"\n'
    )

    # Make the wavelength axis.
    if science_header["CRPIX1"] > 0.5:
        reference_pixel = np.ceil(science_header["CRPIX1"])
        reference_difference = reference_pixel-science_header["CRPIX1"]
    else:
        reference_pixel = np.floor(science_header["CRPIX1"])
        reference_difference = science_header["CRPIX1"] - reference_pixel
    
    reference_wavelength = (
		(reference_difference * science_header["CD1_1"]) 
		+ science_header["CRVAL1"]
	)
	
    wavelength_pixel_axis = (
        np.arange(-reference_pixel+1, science_header["NAXIS1"]-reference_pixel)
    )
    
    wavelength_axis = (
        (wavelength_pixel_axis * science_header["CD1_1"])
        + reference_wavelength
    )

    # If DRAGONS has been used to reduce the data, flip the wavelength
    # axis and the 2D data.
    if science_header["CD1_1"] < 0:
        wavelength_axis = np.flip(wavelength_axis)
        data = np.flip(data, axis=1)
        errs = np.flip(errs, axis=1)
        qual = np.flip(qual, axis=1)

    return data, errs, qual, header_dict, wavelength_axis


def save_gmos(hdu_list, original_hdu_list):
	hdu_list.append(original_hdu_list["MDF"])
