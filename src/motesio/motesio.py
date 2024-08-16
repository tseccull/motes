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
import motes.notes as notes
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
          A list of regions read in from reg.txt in read_regions()

    Returns:
     -- header_dict (dict)
          A dictionary containing parameters and metadata read from the
          header of the image file.
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
    notes.data_harvest_1(data_spatial_floor, data_spatial_ceiling)

    # Create spatial axis for the 2D spectrum and a high resolution
    # version (standard res * 5) for the purposes of plotting
    spat_axis_len = data_sliced_shape[0]
    spatial_axis = np.linspace(0., float(spat_axis_len - 1), num=spat_axis_len)
    
    hi_resolution_spatial_axis = (
        np.linspace(spatial_axis[0], spatial_axis[-1], num=spat_axis_len * 5)
    )

    if (data_regions[region_counter][2] < wavelength_axis[0] or 
            data_regions[region_counter][3] > wavelength_axis[-1]):
        notes.data_harvest_2()
        sys.exit()

    wavelength_slice = (
        np.where(
            np.logical_and(
                wavelength_axis >= data_regions[region_counter][2],
                wavelength_axis <= data_regions[region_counter][3]
            )
        )
    )
    wavelength_start = wavelength_slice[0][0]
    wavelength_end = wavelength_slice[0][-1]
    wavelength_axis = wavelength_axis[wavelength_slice]

    data_sliced = np.squeeze(data_sliced[:, wavelength_slice])
    errs_sliced = np.squeeze(errs_sliced[:, wavelength_slice])
    qual_sliced = np.squeeze(qual_sliced[:, wavelength_slice])
    
    notes.data_harvest_3(
        data_regions, 
        region_counter, 
        header_dict["wavelength_unit"], 
        wavelength_start, 
        len(wavelength_axis)
    )

    frame_dict = {
        "data": data_sliced,
        "errs": errs_sliced,
        "qual": qual_sliced,
    }

    axes_dict = {
        "spataxislen": spat_axis_len,
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


def datasec_string(floor, ceiling):
    """
    Format a string containing the bounds for a section of the spatial
    or wavelength axes that will be saved to the header of MOTES output
    files.
    
    Args:
     -- floor (int)
          Lower bound of data section.
     -- ceiling (int)
          Upper bound of data section.
          
    Returns:
     -- "[" + str(floor) + ":" + str(ceiling) + "]" (str)
          Data section string for addition to a FITS header.
    """
    
    return "[" + str(floor) + ":" + str(ceiling) + "]"
    

def make_bin_table(bin_data, flux_unit, sky=False):
    """
    Take data about the Moffat profile fitted in each bin and convert it
    to FITS Table HDU.
    
    Args:
     -- bin_data (numpy.ndarray)
          array of bin fitting data
     -- flux_unit (str)
          unit (BUNIT) of the data in the science frame
    
    Returns:
     -- bin_table_hdu (astropy.io.fits.hdu.table.BinTableDHU)
          Table HDU containing bin fitting data
    """
    
    bin_data[:,7] -= 1
    bin_table = Table(
        bin_data,
        names=(
            "amplitude",
            "center",
            "alpha",
            "beta",
            "background_level",
            "background_grad",
            "bin_limit_lo",
            "bin_limit_hi"
        ), 
        dtype = ("f8", "f8", "f8", "f8", "f8", "f8", "i4", "i4")
    )
    if sky:
        bin_table_hdu = fits.BinTableHDU(bin_table, name="SBMPARAM")
    else:
        bin_table_hdu = fits.BinTableHDU(bin_table, name="XBMPARAM")
    bin_table_hdu.header["TABUNIT1"] = flux_unit
    bin_table_hdu.header["TABUNIT2"] = "spatial_pixel"
    bin_table_hdu.header["TABUNIT3"] = None
    bin_table_hdu.header["TABUNIT4"] = None
    bin_table_hdu.header["TABUNIT5"] = flux_unit + " / spatial_pixel"
    bin_table_hdu.header["TABUNIT6"] = flux_unit
    bin_table_hdu.header["TABUNIT7"] = "spectral_pixel"
    bin_table_hdu.header["TABUNIT8"] = "spectral_pixel"
	
    return bin_table_hdu


def make_ext_table(extraction_data, wavelength_unit, sky=False):
    """
    Take data about the spectrum extraction limits and convert it to
    FITS Table HDU.
    
    Args:
     -- bin_data (numpy.ndarray)
          array of bin fitting data
     -- wavelength_unit (str)
          unit of the wavelength axis
    
    Returns:
     -- ext_table_hdu (astropy.io.fits.hdu.table.BinTableDHU)
          Table HDU containing the extraction limit data
    """
    
    ext_table = Table(
        extraction_data,
        names=("wavelength", "ext_lim_lo", "ext_lim_hi"), 
        dtype = ("f8", "f8", "f8")
    )
    
    if sky:
        ext_table_hdu = fits.BinTableHDU(ext_table, name="SKYLMTAB")
    else:
        ext_table_hdu = fits.BinTableHDU(ext_table, name="XTLMTAB")
    ext_table_hdu.header["TABUNIT1"] = wavelength_unit
    ext_table_hdu.header["TABUNIT2"] = "pixel"
    ext_table_hdu.header["TABUNIT3"] = "pixel"
	
    return ext_table_hdu


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

    notes.read_motes_param_file_1()

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

    notes.done()

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
        notes.read_regions_1()

        with open("reg.txt", "r", encoding="utf-8") as region_file:
            region_lines = region_file.read().splitlines()
            data_regions = [
                [int(limit) for limit in x.split(",")] for x in region_lines
            ]
        notes.done()

    # Complain and quit MOTES if reg.txt isn't found.
    else:
        notes.read_regions_2(os.getcwd())
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
    moffat_parameters,
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
          An array containing the flux values of the optimally extracted
          1D spectrum.
     -- optimal_1d_errs (numpy.ndarray)
          An array containing the flux errors of the optimally extracted
          1D spectrum.
     -- aperture_1d_data (numpy.ndarray)
          An array containing the flux values of the aperture extracted
          1D spectrum.
     -- aperture_1d_errs (numpy.ndarray)
          An array containing the flux errors of the aperture extracted
          1D spectrum.
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
        round(moffat_parameters[0], 5),
        "moffat profile amplitude"
    )
    primary_header["HIERARCH MOFF C"] = (
        round(moffat_parameters[1] + axes_dict["data_spatial_floor"], 5),
        "moffat profile center"
    )
    primary_header["HIERARCH MOFF ALPHA"] = (
        round(moffat_parameters[2], 5),
        "moffat profile alpha value"
    )
    primary_header["HIERARCH MOFF BETA"] = (
        round(moffat_parameters[3], 5),
        "moffat profile beta value"
    )
    primary_header["HIERARCH MOFF BACK"] = (
        round(moffat_parameters[4], 5),
        "moffat profile background level"
    )
    primary_header["HIERARCH MOFF GRAD"] = (
        round(moffat_parameters[5], 5),
        "moffat profile background slope"
    )
    primary_header["HIERARCH IQ"] = (
        round(
            header_parameters["seeing"] * header_parameters["pixel_resolution"],
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
    
    ext_bin_pars_tabhdu = make_bin_table(
        moffat_parameters_all_bins, header_parameters["flux_unit"]
    )
    ext_limit_array = (
        np.vstack([axes_dict["wavelength_axis"], extraction_limits])
    )
    ext_limit_table_hdu = (
        make_ext_table(ext_limit_array.T, header_parameters["wavelength_unit"])
    )
    
    hdu_list = [
       optimal_1d_datahdu,
       aperture_1d_datahdu,
       ext_bin_pars_tabhdu,
       ext_limit_table_hdu,
    ]
    
    if motes_parameters["-SUBTRACT_SKY"]:
        sky_bin_pars_tabhdu = make_bin_table(
            moffat_parameters_all_sky_bins,
            header_parameters["flux_unit"],
            sky=True
        )
        sky_limit_array = (np.vstack([axes_dict["wavelength_axis"], sky_extraction_limits]))
        sky_limit_table_hdu = make_ext_table(
            sky_limit_array.T, header_parameters["wavelength_unit"], sky=True
        )
    
        sky_model_hdu = fits.ImageHDU(frame_dict["sky_model"])
        sky_model_hdu.header["EXTNAME"] = "2D_SKY"
        sky_model_hdu.header["DATATYPE"] = ("Model Intensity", "Type of Data")
        sky_model_hdu.header["DATASECS"] = (
            datasec_string(axes_dict["data_spatial_floor"], axes_dict["data_spatial_ceiling"]),
            "Data section; spatial axis, pixels"
        )
        sky_model_hdu.header["DATASECW"] = (
            datasec_string(axes_dict["wavelength_start"], axes_dict["wavelength_end"]),
            "Data section; wavelength axis, pixels"
        )
        sky_model_hdu.header["BUNIT"] = header_parameters["flux_unit"]
        sky_model_hdu.header["METHOD"] = (
            motes_parameters["-SKYSUB_MODE"], "MOTES sky subtraction method."
        )
        sky_model_hdu.header["COMMENT"] = (
            "2D_SKY data has same wavelength orientation as OPTI_1D_SPEC."
        )
        
        hdu_list.append(sky_model_hdu)
        hdu_list.append(sky_bin_pars_tabhdu)
        hdu_list.append(sky_limit_table_hdu)
   
    instrument_save_dict = {
        "en06": floydsio.save_floyds,
        "en12": floydsio.save_floyds,
        "FORS2": forsio.save_fors,
        "GMOS-N": gmosio.save_gmos,
        "GMOS-S": gmosio.save_gmos,
        "XSHOOTER": xshooio.save_xshoo
    }
   
    hdu_list = instrument_save_dict[primary_header["INSTRUME"]](
        hdu_list, original_hdu_list
    )

    fits_hdu_list = fits.HDUList(hdu_list)
    fits_hdu_list.writeto("m" + input_file_path.split("/")[-1])
    fits_hdu_list.close()

    notes.save_fits_1(input_file_path)
    
    return None
