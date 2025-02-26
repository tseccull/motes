#!/usr/bin/env python3

"""
motesio.py - Contains the main I/O functions needed for MOTES to
             function and calls the relevant harvester functions needed
             for each specific instrument.
"""

import astropy.io.fits as fits
import copy
import datetime
import logging
import motes.motesio.floydsio as floydsio
import motes.motesio.forsio as forsio
import motes.motesio.gmosio as gmosio
import motes.motesio.xshooio as xshooio
import numpy as np
import os

from astropy.table import Table

logger = logging.getLogger("motes")

def data_harvest(input_file_path, data_regions):
    """
    Extract header metadata and data frames, and repackage them into 
    dictionaries for use in the rest of MOTES.

    Args:
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
    
    [logger.info("header_dict[%s] = %s", key, val) for key, val in header_dict.items()]

    # Slice all dataframes based on the input from reg.txt
    data_shape = np.shape(data)
    data_spatial_floor = int(0 + data_regions[0])
    data_spatial_ceiling = int(data_shape[0] - data_regions[1])

    # Slice off the spatial rows outside the spatial region.
    data_sliced = data[data_spatial_floor : data_spatial_ceiling + 1, :]
    errs_sliced = errs[data_spatial_floor : data_spatial_ceiling + 1, :]
    qual_sliced = qual[data_spatial_floor : data_spatial_ceiling + 1, :]
    data_sliced_shape = np.shape(data_sliced)

    # Create spatial axis for the 2D spectrum and a high resolution
    # version (standard res * 5) for the purposes of plotting
    spat_axis_len = data_sliced_shape[0]
    spatial_axis = np.linspace(0., float(spat_axis_len - 1), num=spat_axis_len)
    
    hi_resolution_spatial_axis = (
        np.linspace(spatial_axis[0], spatial_axis[-1], num=spat_axis_len * 5)
    )
    
    wave_region_too_lo = data_regions[2] < wavelength_axis[0]
    wave_region_too_hi = data_regions[3] > wavelength_axis[-1]

    if wave_region_too_lo or wave_region_too_hi:
        logger.critical(
            "Wavelength limit(s) in 'reg.txt' are beyond the wavelength range of the input data."
        )
        logger.critical("Raising ValueError")
        raise ValueError
        
    in_wavelength_region = np.logical_and(
        wavelength_axis >= data_regions[2],
        wavelength_axis <= data_regions[3]
    )
    wavelength_slice = np.where(in_wavelength_region)
    
    wavelength_start = wavelength_slice[0][0]
    wavelength_end = wavelength_slice[0][-1]
    wavelength_axis = wavelength_axis[wavelength_slice]

    data_sliced = np.squeeze(data_sliced[:, wavelength_slice])
    errs_sliced = np.squeeze(errs_sliced[:, wavelength_slice])
    qual_sliced = np.squeeze(qual_sliced[:, wavelength_slice])
    
    logger.info(
        "2D data frame sliceed to [%s:%s,%s:%s]", wavelength_start, 
        wavelength_end + 1,data_spatial_floor, data_spatial_ceiling+1
    )
    
    logger.info(
        "Requested wavelength range is %s-%s %s.", 
        data_regions[2], data_regions[3], header_dict["wavelength_unit"]
    )
    logger.info(
        "Returned wavelength range is %s-%s %s.", wavelength_axis[0], 
        wavelength_axis[-1], header_dict["wavelength_unit"]
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


def read_regions(cwd):
    """
    Search for, and read in, reg.txt file in root working directory. 
    Reads an input line from reg.txt and returns a list of integers 
    defining the boundaries of the region of the 2D data that contains
    the spectrum to be extracted. For each file, the first two integers
    are the number of pixel rows to remove from each end of the spatial
    axis of the spectrum, and the last two are the upper and lower
    wavelength bounds of the region on the dispersion axis.

    Args:
     -- cwd (str)
          The absolute path to the current working directory for MOTES.

    Returns:
     -- data_region (list)
          A list for each file contains a list of integers that define
          the boundaries of the region of the 2D spectum that will be
          used for the extraction.
    """
    
    reg_path = "./inputs/reg.txt"
    if os.path.exists(reg_path) and os.path.isfile(reg_path):
        logger.info("./inputs/reg.txt file found.")

        with open(reg_path, "r", encoding="utf-8") as region_file:
            region_lines = region_file.read().splitlines()
            region_lines = [x.strip() for x in region_lines]
            unique_lines = set()
            duplicate_lines = [[x, i] for i, x in enumerate(region_lines) if x in unique_lines or unique_lines.add(x)]
            
            if len(duplicate_lines) > 0:
                [logger.warning("Ignoring duplicate region %s on line %s of reg.txt", x[0], x[1]) for x in duplicate_lines]
            
            data_regions = []
            for i, line in enumerate(unique_lines):
                if line[-5:] != ".fits":
                    logger.critical("File name string on line %s of reg.txt is not '.fits'.", i+1)
                    logger.critical(line)
                    logger.critical("MOTES accepts only .fits files. Perhaps this is a typo?")
                    logger.critical("Raising RuntimeError.")
                    raise RuntimeError
                if len(line.split(",")) != 5:
                    logger.critical("Number of comma-separated variables on line %s of reg.txt is not five.", i+1)
                    logger.critical(line)
                    logger.critical("See the section on reg.txt in the MOTES documentation for the required format.")
                    logger.critical("Raising RuntimeError.")
                    raise RuntimeError
                if not all([len(x) for x in line.split(",")]) > 0:
                    logger.critical("Variable(s) on line %s of reg.txt have zero length.", i+1)
                    logger.critical(line)
                    logger.critical("Raising RuntimeError.")
                    raise RuntimeError								
                try:
                    data_regions.append([int(limit) for limit in line.split(",")[:4]])
                except ValueError:
                    logger.critical("Non-numeric character found in numeric region variables on line %s of reg.txt.", i+1)
                    logger.critical(line)
                    logger.critical("Raising ValueError.")
                    raise ValueError					    					
            
            [data_regions[x].append(y.split(",")[-1]) for x, y in enumerate(unique_lines)]
        
        logger.info("%s regions read from ./inputs/reg.txt. Basic format checks passed.", len(unique_lines))

    else:
        logger.critical("reg.txt region file not found at " + cwd +  "/inputs/")
        logger.critical("reg.txt is needed to define which part(s) of each input 2D spectrum file MOTES will process.")
        logger.critical("See documentation for more details.")
        logger.critical("Raising FileNotFoundError.")
        raise FileNotFoundError

    return data_regions


def save_fits(
    original_hdu_list,
    axes_dict,
    header_parameters,
    extracted,
    motes_parameters,
    input_file_name,
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
     -- extracted (list)
	      List of numpy arrays containing the extracted spectrum data
	      and the associated uncertainties.
     -- motes_parameters (dict)
          A dictionary containing the MOTES parameters.
     -- input_file_name (str)
          The file name of the 1D spectrum.
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
    
    ut_now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
    spatial_floor = axes_dict["data_spatial_floor"]
    spatial_ceiling = axes_dict["data_spatial_ceiling"]
    wave_lo = np.floor(axes_dict["wavelength_axis"][0])
    wave_hi = np.ceil(axes_dict["wavelength_axis"][-1])
    
    save_file_name = (
        "m" 
        + input_file_name 
        + "_" + str(int(wave_lo)) 
        + ":" + str(int(wave_hi)) 
        + "_" + str(spatial_floor) 
        + ":" + str(spatial_ceiling) 
        + ".fits"
    )
    logger.info("Saving extracted spectrum and metadata to %s", save_file_name)
    
    
    wave_pix_lo = axes_dict["wavelength_start"]
    wave_pix_hi = axes_dict["wavelength_end"]
    wave_unit = header_parameters["wavelength_unit"]
    moffat_amplitude = round(moffat_parameters[0], 5)
    moffat_center = round(moffat_parameters[1] + axes_dict["data_spatial_floor"], 5)
    moffat_alpha = round(moffat_parameters[2], 5)
    moffat_beta = round(moffat_parameters[3], 5)
    moffat_bg_level = round(moffat_parameters[4], 5)
    moffat_bg_grad = round(moffat_parameters[5], 5)
    plate_scale = round(header_parameters["pixel_resolution"], 2)
    pixel_fwhm = round(header_parameters["seeing"], 2)
    motes_fwhm = round(header_parameters["seeing"] * header_parameters["pixel_resolution"], 2)
    col_bin_lim = int(motes_parameters.minimum_column_limit)
    fwhm_multiplier = motes_parameters.extraction_fwhm_multiplier
    spectrum_unit = header_parameters["flux_unit"]
    wave_axis = axes_dict["wavelength_axis"]
    
    head_card_input = [
        ["MOTES", "motes.py", "Extraction script"],
        ["MOTESV", "v1.0.0", "MOTES Version"],
        ["MOTESDOI", "UNKNOWN", "MOTES DOI"],
        ["UTXTIME", ut_now, "UT timestamp for MOTES"],
        ["SPATPIXL", spatial_floor, "lower limit of spatial axis, pix"],
        ["SPATPIXH", spatial_ceiling, "upper limit of spatial axis, pix"],
        ["DISPPIXL", wave_pix_lo, "lower limit of dispersion axis, pix"],
        ["DISPPIXH", wave_pix_hi, "upper limit of dispersion axis, pix"],
        ["WAVL", wave_lo, "lower limit of wav range, " + wave_unit],
        ["WAVH", wave_hi, "upper limit of wav range, " + wave_unit],
        ["WAVU", wave_unit, "wavelength unit"],
        ["MOFFA", moffat_amplitude, "median Moffat profile amplitude"],
        ["MOFFC", moffat_center, "median Moffat profile center"],
        ["MOFFALPH", moffat_alpha, "median Moffat profile alpha value"],
        ["MOFFBETA", moffat_beta, "median Moffat profile beta value"],
        ["MOFFBGLV", moffat_bg_level, "median Moffat profile background level"],
        ["MOFFBGSL", moffat_bg_grad, "median Moffat profile background slope"],
        ["PLATESCL", plate_scale, "spatial pixel plate scale, arcsec/pix"],
        ["PIXELIQ", pixel_fwhm, "IQ measured from median profile, pix"],
        ["ARCSECIQ", motes_fwhm, "IQ measured from median profile, arcsec"],
        ["SNRBNLIM", motes_parameters.extraction_snr_limit, "maximum SNR per bin"],
        ["COLBNLIM", col_bin_lim, "minimum number of columns per bin"],
        ["FWHMMULT", fwhm_multiplier, "FWHM used to define the extraction limits"]
    ]

    if motes_parameters.subtract_sky:
        sky_fwhm_mult = motes_parameters.sky_fwhm_multiplier
        sky_snr_bin_lim = motes_parameters.sky_snr_limit
        sky_order = motes_parameters.sky_order
        
        head_card_input.append(["SFWHMMLT", sky_fwhm_mult, "FWHM multiplier to define sky region"])
        head_card_input.append(["SSNRBNLM", sky_snr_bin_lim, "max SNR per bin for sky subtraction"])
        head_card_input.append(["SKYORDER", sky_order, "polynomial order of spatial sky model"])

    head_card_input.append(["HDUROW0", "Wavelength Axis, " + wave_unit, ""])
    head_card_input.append(["HDUROW1", "Spectrum data, " + spectrum_unit, ""])
    head_card_input.append(["HDUROW2", "Spectrum uncertainty, " + spectrum_unit, ""])
    head_card_input.append(["EXTNAME", "OPTI_1D_SPEC", ""])
    
    primary_header = original_hdu_list[0].header
    for card_deets in head_card_input:
        primary_header[card_deets[0]] = (card_deets[1], card_deets[2])
    
    optimal_1d_datahdu = fits.PrimaryHDU(
        [wave_axis, extracted[0], extracted[1]], header=primary_header
    )
    
    aperture_1d_datahdu = fits.ImageHDU(
        [wave_axis, extracted[2], extracted[3]], header=primary_header
    )
    aperture_1d_datahdu.header["EXTNAME"] = "APER_1D_SPEC"
    
    ext_bin_pars_tabhdu = make_bin_table(moffat_parameters_all_bins, spectrum_unit)
    ext_limit_array = np.vstack([wave_axis, extraction_limits])
    ext_limit_table_hdu = (make_ext_table(ext_limit_array.T, wave_unit))
    
    hdu_list = [optimal_1d_datahdu, aperture_1d_datahdu, ext_bin_pars_tabhdu, ext_limit_table_hdu]
    
    if motes_parameters.subtract_sky:
        spatial_datasec = datasec_string(spatial_floor, spatial_ceiling)
        wavelen_datasec = datasec_string(wave_pix_lo, wave_pix_hi)
        sky_order = motes_parameters.sky_order
		
        sky_head_card_info = [
            ["DATATYPE", "Model Intensity", "Type of data"],
            ["DATASECS", spatial_datasec, "Data section; spatial axis, pixels"],
            ["DATASECW", wavelen_datasec, "Data section; wavelength axis, pixels"],
            ["BUNIT", spectrum_unit, ""],
            ["SKYORDER", int(sky_order), "MOTES sky characterization order."],
            ["COMMENT", "2D_SKY data has same wavelength orientation as OPTI_1D_SPEC.", ""]
        ]
        
        sky_bin_pars_tabhdu = make_bin_table(
            moffat_parameters_all_sky_bins, spectrum_unit, sky=True
        )
        
        sky_limit_array = (np.vstack([wave_axis, sky_extraction_limits]))
        sky_limit_table_hdu = make_ext_table(sky_limit_array.T, wave_unit, sky=True)
    
        sky_model_hdu = fits.ImageHDU(frame_dict["sky_model"])
        sky_uncertainty_hdu = fits.ImageHDU(frame_dict["sky_uncertainty"])
        for card_deets in sky_head_card_info:
            sky_model_hdu.header[card_deets[0]] = (card_deets[1], card_deets[2])
            sky_uncertainty_hdu.header[card_deets[0]] = (card_deets[1], card_deets[2])
        
        sky_model_hdu.header["EXTNAME"] = ("2D_SKY", "")
        sky_uncertainty_hdu.header["EXTNAME"] = ("2D_SKY_UNC", "")
        
        hdu_list.append(sky_model_hdu)
        hdu_list.append(sky_uncertainty_hdu)
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
    fits_hdu_list.writeto(save_file_name)
    fits_hdu_list.close()

    logger.info("Data saved.")
    
    return None
