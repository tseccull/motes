#!/home/tom/miniforge3/envs/work/bin/python

"""
motesio.py - Contains the main I/O functions needed for MOTES to
             function and calls the relevant harvester functions needed
             for each specific instrument.
"""

import astropy.io.fits as fits
import os
import motesio.floydsio as floydsio
import motesio.forsio as forsio
import motesio.gmosio as gmosio
import motesio.xshooio as xshooio
import numpy as np
import sys


def data_harvest(region_counter, input_file_path, data_regions):
    """
    Extract header metadata and data frames, and repackage them into 
    dictionaries for use in the rest of MOTES.

    Args:
     -- region_counter (int) 
            An integer noting which line in region is being called to
            define the boundaries of the 2D data.
     -- input_file_path (str) 
            The name of the data file.
     -- data_regions (list)
            A list of regions read in from reg.txt in
            startup.read_regions()

    Returns:
     -- header_dict (dict)
            A dictionary containing parameters and metadata read from
            the header of the image file.
     -- frame_dict (dict)
            A dictionary containing the 2D data frames read from the
            image file.
     -- axes_dict (dict)
            A dictionary containing the spatial and spectral axis arrays
            associated with the data frames, along with metadata used to
            define the boundaries of the 2D data.
     -- input_file_primary_header (dict)
            A copy of the image file header; this is also a dictionary.
    """

    # Create dictionary to tell data_harvest which instrument specific
    # function to call.
    instrument_dict = {
        "en06": floydsio.harvest_floyds,
        "en12": floydsio.harvest_floyds,
        "FORS2": forsio.harvest_fors2,
        "GMOS-N": gmosio.harvest_gmos,
        "GMOS-S": gmosio.harvest_gmos,
        "XSHOOTER": xshooio.harvest_xshoo
    }

    # Open file containing the spectral data, then extract the header,
    # image frame, error frame, and quality frame (if the file has one).
    with fits.open(input_file_path) as input_fits_hdu:
        input_file_primary_header = input_fits_hdu[0].header
        instrument = input_file_primary_header["INSTRUME"]
        # Based on the value of instrument, this calls one of the
        # harvest_instrument functions.
        (
            data,
            errs,
            qual,
            original_qual,
            header_dict,
            wavelength_axis,
        ) = instrument_dict[
            instrument
        ](input_fits_hdu, input_file_primary_header)

    # Slice all dataframes based on the input from reg.txt
    data_shape = np.shape(data)
    data_spatial_floor = int(0 + data_regions[region_counter][0])
    data_spatial_ceiling = int(data_shape[0] - data_regions[region_counter][1])

    # Slice off the spatial rows outside the spatial region.
    data_sliced = data[data_spatial_floor : data_spatial_ceiling + 1, :]
    errs_sliced = errs[data_spatial_floor : data_spatial_ceiling + 1, :]
    qual_sliced = qual[data_spatial_floor : data_spatial_ceiling + 1, :]
    data_sliced_shape = np.shape(data_sliced)
    sys.stdout.write(
        " >>> 2D spectrum sliced on spatial axis based on user defined limits:\n"
        "     New spatial axis covers pixel rows "
        + str(data_spatial_floor)
        + "-"
        + str(data_spatial_ceiling)
        + ".\n"
    )

    # Create spatial axis for the 2D spectrum and a high resolution
    # version (standard res * 5) for the purposes of plotting
    spatial_axis = (
        np.linspace(
            0.0, float(data_sliced_shape[0] - 1), num=data_sliced_shape[0]
        )
    )
    hi_resolution_spatial_axis = (
        np.linspace(
            spatial_axis[0], spatial_axis[-1], num=len(spatial_axis) * 5
        )
    )

    if (data_regions[region_counter][2] < wavelength_axis[0] or 
            data_regions[region_counter][3] > wavelength_axis[-1]):
        sys.stdout.write(
            " >>> User defined wavelength limit(s) are outside native wavelength range\n"
            '     Make sure "-LOW_WAV_SLICE" > lower limit of wavelength axis\n'
            '     Make sure "-HIGH_WAV_SLICE" < upper limit of wavelength axis\n'
            "     Terminating MOTES.\n\n"
        )
        sys.exit()

    wavelength_slice = np.where(
        np.logical_and(
            wavelength_axis >= data_regions[region_counter][2],
            wavelength_axis <= data_regions[region_counter][3]
        )
    )
    wavelength_start = wavelength_slice[0][0]
    wavelength_end = wavelength_slice[0][-1]
    wavelength_axis = wavelength_axis[wavelength_slice]

    data_sliced = np.squeeze(data_sliced[:, wavelength_slice])
    errs_sliced = np.squeeze(errs_sliced[:, wavelength_slice])
    qual_sliced = np.squeeze(qual_sliced[:, wavelength_slice])

    sys.stdout.write(
        " >>> 2D spectrum sliced on dispersion axis based on user defined limits:\n"
        "     New Wavelength range is "
        + str(data_regions[region_counter][2])
        + "-"
        + str(data_regions[region_counter][3])
        + " "
        + header_dict["wavelength_unit"]
        + ".\n"
        "     This range is equivalent to pixel columns "
        + str(wavelength_start)
        + "-"
        + str(wavelength_start + len(wavelength_axis))
        + "\n"
    )

    frame_dict = {
        "data": data_sliced,
        "errs": errs_sliced,
        "qual": qual_sliced,
        "original_qual": original_qual,
    }

    axes_dict = {
        "spataxislen": len(spatial_axis),
        "spatial_axis": spatial_axis,
        "hi_resolution_spatial_axis": hi_resolution_spatial_axis,
        "data_spatial_floor": data_spatial_floor,
        "data_spatial_ceiling": data_spatial_ceiling,
        "dispersion_axis_length": len(wavelength_axis),
        "wavelength_axis": wavelength_axis,
        "wavelength_start": wavelength_start,
        "wavelength_end": wavelength_end,
    }

    return header_dict, frame_dict, axes_dict, input_file_primary_header


def read_motes_parameter_file():
    """
    Import parameters from motesparams.txt parameter file into a 
    dictionary.

    Args:
     -- None

    Returns:
     -- parameter_dict (dict) 
         a dictionary containing the parameters read in from 
         motesparams.txt.
    """

    sys.stdout.write(" >>> Reading in parameters from motesparams.txt. ")
    sys.stdout.flush()

    # Read in MOTES parameter file line by line and filter out the empty lines.
    with open("motesparams.txt", "r", encoding="utf-8") as parameter_file:
        parameter_lines = parameter_file.read().splitlines()
        parameter_lines = filter(None, parameter_lines)

    # Flatten the 2D list of parameters and keywords into a 1D list where each
    # parameter's value follows its associated keyword.
    lumpy_parameter_list = [
	    x.split("=") for x in parameter_lines if x[0] == "-"
	]
    flat_parameter_list = [y for x in lumpy_parameter_list for y in x]

    # Convert all numerical values in the parameter list to floats.
    # If digit, convert to float. If not, leave as string.
    parameter_list = []
    for i in flat_parameter_list:
        if i.replace(".", "", 1).isdigit():
            parameter_list.append(float(i))
        else:
            parameter_list.append(i)

    # Assign parameters and their associated keywords to a dictionary.
    parameter_dict = dict(zip(parameter_list[::2], parameter_list[1::2]))

    sys.stdout.write("DONE.\n")

    return parameter_dict


def read_regions():
    """
    Search for, and read in, reg.txt file in root working directory. 
    Reads an input line from reg.txt and returns a list of integers 
    defining the boundaries of the region of the 2D data that contains
    the spectrum to be extracted. For each file, the first two integers
    are the number of pixel rows to remove from each end of the spatial
    axis of the spectrum, and the last two are the upper and lower
    wavelength bounds of the region on the dispersion axis.

    Args:
     -- None

    Returns:
     -- data_region (list)
         A list for each file contains a list of integers that define
         the boundaries of the region of the 2D spectum that will be
         used for the extraction.
    """

    # Search for reg.txt and read in the list that it contains.
    if os.path.exists("reg.txt") and os.path.isfile("reg.txt"):
        sys.stdout.write(" >>> reg.txt file found. Reading reg.txt. ")
        sys.stdout.flush()

        with open("reg.txt", "r", encoding="utf-8") as region_file:
            region_lines = region_file.read().splitlines()
            data_regions = [
                [int(limit) for limit in x.split(",")] for x in region_lines
            ]
        sys.stdout.write("DONE.\n")

    # Complain and quit MOTES if reg.txt isn't found.
    else:
        sys.stdout.write(" >>> reg.txt file not found in root working directory.\n")
        sys.stdout.write("     Root working directory is:\n")
        sys.stdout.write("     " + os.getcwd() + "\n")
        sys.stdout.write(
            "     Please provide a reg.txt file to tell MOTES the spatial extent of the 2D\n"
        )
        sys.stdout.write("     spectra. See the docs for further information.\n")
        sys.stdout.write("     Terminating MOTES.\n\n")
        sys.exit()

    return data_regions
