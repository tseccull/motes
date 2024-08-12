#!/home/tom/miniforge3/envs/work/bin/python

"""
motesio.py - Contains the main I/O functions needed for MOTES to
             function and calls the relevant harvester functions needed
             for each specific instrument.
"""

import astropy.io.fits as fits
import astropy.table as Table
import copy
import datetime
import os
import motesio.floydsio as floydsio
import motesio.forsio as forsio
import motesio.gmosio as gmosio
import motesio.xshooio as xshooio
import numpy as np
import sys

from astropy.table import Table


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
    instrument_harvest_dict = {
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
        instrument = input_fits_hdu[0].header["INSTRUME"]
        original_hdu_list = copy.deepcopy(input_fits_hdu)
        # Based on the value of instrument, this calls one of the
        # harvest_instrument functions.
        data, errs, qual, header_dict, wavelength_axis = (
            instrument_harvest_dict[instrument](input_fits_hdu)
        )

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

    return header_dict, frame_dict, axes_dict, original_hdu_list


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


def save_fits(
    original_hdu_list,
    axes_dict,
    header_parameters,
    optimal_1d_data,
    optimal_1d_errs,
    aperture_1d_data,
    aperture_1d_errs,
    motes_parameters,
    input_file_path,
    moffat_profile_parameters,
    frame_dict,
    moffat_parameters_all_bins,
    extraction_limits,
    moffat_parameters_all_sky_bins,
    sky_extraction_limits
):
    """
    This function saves the extracted spectrum and intermediate products
    in a single, newly constructed, FITS file.

    Args:
     -- original_hdu_list (astropy.io.fits.hdu.hdulist.HDUList)
            Astropy HDU List read in from the original input file.
     -- axes_dict (dict)
            A dictionary containing the axes information.
     -- header_parameters (dict)
            A dictionary containing the header information.
     -- optimal_1d_data (numpy.ndarray)
            An array containing the flux values of the optimally
            extracted 1D spectrum.
     -- optimal_1d_errs (numpy.ndarray)
            An array containing the flux errors of the optimally
            extracted 1D spectrum.
     -- aperture_1d_data (numpy.ndarray)
            An array containing the flux values of the aperture
            extracted 1D spectrum.
     -- aperture_1d_errs (numpy.ndarray)
            An array containing the flux errors of the aperture
            extracted 1D spectrum.
     -- motes_parameters (dict)
            A dictionary containing the MOTES parameters.
     -- input_file_path (str)
            The filename of the 1D spectrum.
     -- moffat_profile_parameters (list)
            A list containing the Moffat fit parameters.
     -- frame_dict (dict)
            A dictionary containing the original 2D spectrum data and
            error frames.
     -- moffat_parameters_all_bins (numpy.ndarray)
            A dictionary containing the binning parameters.
     -- extraction_limits (numpy.ndarray)
            An array containing the extraction limits.
     -- moffat_parameters_all_sky_bins (numpy.ndarray)
            An array containing the binning parameters for the sky
            extraction.
     -- sky_extraction_limits (list)
            A list containing the extraction limits for the sky
            extraction.

    Returns:
        None
    """

    primary_header = original_hdu_list[0].header
    primary_header["MOTES"] = ("motes.py", "Extraction script")
    primary_header["MOTESV"] = ("v0.5.0", "MOTES Version")
    primary_header["MOTESDOI"] = ("UNKNOWN", "MOTES DOI")
    primary_header["HIERARCH UTC EXT DATE"] = (
        datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
        "UT timestamp for extraction"
    )
    primary_header["HIERARCH SPATPIXL"] = (
        axes_dict["data_spatial_floor"],
        "lower limit of spatial axis, pix"
    )
    primary_header["HIERARCH SPATPIXH"] = (
        axes_dict["data_spatial_ceiling"],
        "upper limit of spatial axis, pix"
    )
    primary_header["HIERARCH DISPPIXL"] = (
        axes_dict["wavelength_start"],
        "lower limit of dispersion axis, pix"
    )
    primary_header["HIERARCH DISPPIXH"] = (
        axes_dict["wavelength_end"],
        "upper limit of dispersion axis, pix"
    )
    primary_header["HIERARCH WAVL"] = (
        np.floor(axes_dict["wavelength_axis"][0]),
        "lower limit of wav range, " + header_parameters["wavelength_unit"]
    )
    primary_header["HIERARCH WAVH"] = (
        np.ceil(axes_dict["wavelength_axis"][-1]),
        "upper limit of wav range, " + header_parameters["wavelength_unit"]
    )
    primary_header["HIERARCH WAVU"] = (
        header_parameters["wavelength_unit"],
        "Wavelength unit"
    )
    primary_header["HIERARCH MOFF A"] = (
        round(moffat_profile_parameters[0], 5),
        "moffat profile amplitude"
    )
    primary_header.add_blank(
        "Parameters fit to the median spatial profile of the spectrum",
        before="HIERARCH MOFF A"
    )
    primary_header["HIERARCH MOFF C"] = (
        round(
            moffat_profile_parameters[1] + axes_dict["data_spatial_floor"], 5
        ),
        "moffat profile center"
    )
    primary_header["HIERARCH MOFF ALPHA"] = (
        round(moffat_profile_parameters[2], 5),
        "moffat profile alpha value"
    )
    primary_header["HIERARCH MOFF BETA"] = (
        round(moffat_profile_parameters[3], 5),
        "moffat profile beta value"
    )
    primary_header["HIERARCH MOFF BACK"] = (
        round(moffat_profile_parameters[4], 5),
        "moffat profile background level"
    )
    primary_header["HIERARCH MOFF GRAD"] = (
        round(moffat_profile_parameters[5], 5),
        "moffat profile background slope"
    )
    primary_header["HIERARCH IQ"] = (
        round(
            header_parameters["seeing"]* header_parameters["pixel_resolution"],
            2
        ),
        'IQ measured from median profile, "'
    )
    primary_header["HIERARCH SNR BIN LIMIT"] = (
        motes_parameters["-SNR_BIN_LIM"],
        "maximum SNR per bin"
    )
    primary_header.add_blank(
        "Dispersion Binning and Spectrum Extraction",
        before="HIERARCH SNR BIN LIMIT"
    )
    primary_header["HIERARCH COL BIN LIMIT"] = (
        int(motes_parameters["-COL_BIN_LIM"]),
        "minimum number of columns per bin"
    )
    primary_header["HIERARCH FWHM MULTIPLIER"] = (
        motes_parameters["-FWHM_MULTIPLIER"],
        "FWHM used to define the extraction limits"
    )

    if motes_parameters["-SUBTRACT_SKY"]:
        primary_header["HIERARCH SKYSUB FWHM MULT"] = (
            motes_parameters["-BG_FWHM_MULTIPLIER"],
            "FWHM multiplier for defining background"
        )
        primary_header["HIERARCH SKYSUB SNR BIN LIM"] = (
            motes_parameters["-SKY_SNR_BIN_LIM"],
            "max SNR per bin for sky subtraction"
        )
        sky_model_hdu = fits.ImageHDU(frame_dict["sky_model"])
        sky_model_hdu.header["EXTNAME"] = "2D_SKY"
        sky_bin_hdu = fits.ImageHDU(moffat_parameters_all_sky_bins)
        sky_bin_hdu.header["EXTNAME"] = "SKY_BIN_PARS"
        sky_extraction_limits = fits.ImageHDU(sky_extraction_limits)
        sky_extraction_limits.header["EXTNAME"] = "SKY_EXT_LIMS"

    primary_header["HIERARCH EXTRACTED HDU ROW 0"] = (
        "Wavelength Axis, " + header_parameters["wavelength_unit"]
    )
    primary_header["HIERARCH EXTRACTED HDU ROW 1"] = (
        "Flux, " + header_parameters["flux_unit"]
    )
    primary_header["HIERARCH EXTRACTED HDU ROW 2"] = (
        "Flux Uncertainty, " + header_parameters["flux_unit"]
    )
    
    primary_header["EXTNAME"] = "OPTI_1D_SPEC"
    optimal_1d_datahdu = fits.PrimaryHDU(
        [axes_dict["wavelength_axis"], optimal_1d_data, optimal_1d_errs],
        header=primary_header
    )
    
    aperture_1d_datahdu = fits.ImageHDU(
        [axes_dict["wavelength_axis"], aperture_1d_data, aperture_1d_errs],
        header=primary_header
    )
    aperture_1d_datahdu.header["EXTNAME"] = "APER_1D_SPEC"
    
    ext_bin_moffat_parameters = Table(
        moffat_parameters_all_bins,
        names=(
            "amplitude",
            "center",
            "alpha",
            "beta",
            "background_grad",
            "background_level",
            "bin_limit_lo",
            "bin_limit_hi"
        ), 
        dtype = ("f8", "f8", "f8", "f8", "f8", "f8", "i4", "i4"),
        units = (
            header_parameters["flux_unit"] + ".",
            "spatial_pixel",
            "",
            "",
            header_parameters["flux_unit"] + " / spatial_pixel",
            header_parameters["flux_unit"] + ".",
            "spectral pixel",
            "spectral pixel"
        )
    )
    ext_bin_pars_tabhdu = fits.BinTableHDU(
        ext_bin_moffat_parameters, name="XBMPARAM"
    )
    
    print(header_parameters["flux_unit"])
    print(ext_bin_pars_tabhdu.header)
    exit()
    
    if motes_parameters["-SUBTRACT_SKY"]:
        sky_bin_moffat_parameters = Table(
            moffat_parameters_all_bins,
            names=(
                "amplitude",
                "center",
                "alpha",
                "beta",
                "background_grad",
                "background_level",
                "bin_limit_lo",
                "bin_limit_hi"
            ), 
            dtype=("f8", "f8", "f8", "f8", "f8", "f8", "i4", "i4")
        )
        sky_bin_pars_tabhdu = fits.BinTableHDU(
            ext_bin_moffat_parameters, name="SBMPARAM"
        )
    
    hdu_list = [
       optimal_1d_datahdu,
       aperture_1d_datahdu
       
    ]
    
    moffat_parameters_all_bins[:,7] -= 1
    moffat_parameters_all_sky_bins[:,7] -= 1
    print(Table(moffat_parameters_all_bins, names=("amplitude", "center", "alpha", "beta", "background_grad", "background_level", "bin_limit_lo", "bin_limit_hi"), dtype=("f8", "f8", "f8", "f8", "f8", "f8", "i4", "i4")))
    print(Table(moffat_parameters_all_sky_bins, names=("amplitude", "center", "alpha", "beta", "background_grad", "background_level", "bin limit lo", "bin limit hi"), dtype=("f8", "f8", "f8", "f8", "f8", "f8", "i4", "i4")))
    print(np.shape(moffat_parameters_all_bins))
    exit()
    
    bins_moffat_parameters_hdu = fits.ImageHDU(moffat_parameters_all_bins)
    bins_moffat_parameters_hdu.header["EXTNAME"] = "EXT_BIN_PARS"
    extraction_limits = fits.ImageHDU(extraction_limits)
    extraction_limits.header["EXTNAME"] = "EXT_LIMS"
   
    instrument_save_dict = {
        "en06": floydsio.save_floyds,
        "en12": floydsio.save_floyds,
        "FORS2": forsio.save_fors2,
        "GMOS-N": gmosio.save_gmos,
        "GMOS-S": gmosio.save_gmos,
        "XSHOOTER": xshooio.save_xshoo
    }
   
    hdu_list = instrument_save_dict[instrument](hdu_list, original_hdu_list)
    
    hdu_list = [
        optimal_1d_datahdu,
        aperture_1d_datahdu,
        orig_2d_spec_hdu,
        orig_2d_errs_hdu,
        orig_2d_qual_hdu,
        bins_moffat_parameters_hdu,
        extraction_limits,
    ]

    if motes_parameters["-SUBTRACT_SKY"]:
        hdu_list.append(sky_model_hdu)
        hdu_list.append(sky_bin_hdu)
        hdu_list.append(sky_extraction_limits)

    fits_hdu_list = fits.HDUList(hdu_list)
    fits_hdu_list.writeto("m" + input_file_path.split("/")[-1])
    fits_hdu_list.close()

    sys.stdout.write(" >>> Spectrum extracted and saved:\n")
    sys.stdout.write(
        "     " + "/".join(input_file_path.split("/")[0:-1]) + "/m" + input_file_path.split("/")[-1] + "\n"
    )
    return None
