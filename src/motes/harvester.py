"""
harvester.py - A collection of functions to read in data from the image file and repackage it into
               a dictionary for use in the rest of MOTES.
"""

import copy
import sys

import astropy.io.fits as fits
import numpy as np

import motes.common as common


def data_harvest(region_counter, input_file_path, data_regions):
    """Extract header metadata and data frames, and repackage them into dictionaries for use in the
       rest of MOTES.

    Args:
        region_counter (int) : An integer noting which line in region is being called to define the
                            boundaries of the 2D data.
        input_file_path (str) : The name of the data file.
        data_regions (list)     : A list of regions read in from reg.txt in startup.read_regions()

    Returns:
        header_dict (dict)  : A dictionary containing parameters and metadata read from the header of
                            the image file.
        frame_dict (dict) : A dictionary containing the 2D data frames read from the image file.
        axes_dict (dict)  : A dictionary containing the spatial and spectral axis arrays associated
                            with the data frames, along with metadata used to define the boundaries
                            of the 2D data.
        input_file_primary_header (dict)    : A copy of the image file header; this is also a dictionary.
    """

    # Create dictionary to tell data_harvest which instrument specific function to call.
    instrument_dict = {
        "en06": harvest_floyds,
        "en12": harvest_floyds,
        "FORS2": harvest_fors2,
        "GMOS-N": harvest_gmos,
        "GMOS-S": harvest_gmos,
        "XSHOOTER": harvest_xshoo,
    }

    # Open file containing the spectral data, then extract the header, image frame, error frame,
    # and quality frame (if the file has one).
    with fits.open(input_file_path) as input_fits_hdu:
        input_file_primary_header = input_fits_hdu[0].header
        instrument = input_file_primary_header["INSTRUME"]
        # Based on the value of instrument, this calls one of the harvest_instrument functions.
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

    # Create spatial axis for the 2D spectrum and a high resolution version (standard res * 5) for
    # the purposes of plotting
    spatial_axis = np.linspace(0.0, float(data_sliced_shape[0] - 1), num=data_sliced_shape[0])
    hi_resolution_spatial_axis = np.linspace(spatial_axis[0], spatial_axis[-1], num=len(spatial_axis) * 5)

    if data_regions[region_counter][2] < wavelength_axis[0] or data_regions[region_counter][3] > wavelength_axis[-1]:
        sys.stdout.write(
            " >>> User defined wavelength limit(s) are outside native wavelength range\n"
            '     Make sure "-LOW_WAV_SLICE" > lower limit of wavelength axis\n'
            '     Make sure "-HIGH_WAV_SLICE" < upper limit of wavelength axis\n'
            "     Terminating MOTES.\n\n"
        )
        sys.exit()

    wavelength_slice = np.where(
        np.logical_and(
            wavelength_axis >= data_regions[region_counter][2], wavelength_axis <= data_regions[region_counter][3]
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


def harvest_gmos(input_fits_hdu, primary_header):
    """
    Harvest the header and data from a GMOS spectrum.

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

    # Retrieve the data frame, error frame, and qual frame. Also retrieve the header of the science
    # image frame, as some metadata is stored there instead of the primary header.
    science_header = input_fits_hdu["SCI"].header
    data = input_fits_hdu["SCI"].data
    errs = input_fits_hdu["VAR"].data ** 0.5
    qual = input_fits_hdu["DQ"].data
    original_qual = copy.deepcopy(qual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    # Pixels in chip gaps are kept as 1 to make sure they don't get flagged as bad.
    qual[np.where(data + qual == 1)] = 0
    qual = 1 - qual

    # Set pixels with NaN value to 1 in the data frame, and flag them as bad pixels in the qual
    # frame.
    qual[~np.isfinite(data)] = 0
    data[~np.isfinite(data)] = 1.0

    # All this is to get an initial estimate of the IQ. Tables below are based on the condition
    # constraints used by Gemini.
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
        "pixel_resolution": float(primary_header["PIXSCALE"]),
        "exptime": primary_header["EXPTIME"],
        "seeing": seeing,
        "instrument": primary_header["INSTRUME"],
        "wavelength_unit": science_header["WAT1_001"].split(" ")[2].split("=")[1],
    }

    # BUNIT only appears in the headers of GMOS spectra if they have been flux calibrated.
    if "BUNIT" in science_header:
        header_dict["flux_unit"] = science_header["BUNIT"]
    else:
        header_dict["flux_unit"] = "electrons"

    sys.stdout.write("DONE.\n")
    sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(header_dict["pixel_resolution"])
        + '"\n'
    )

    # Create the wavelength axis of the spectrum.
    wavelength_axis = common.make_wav_axis(
        science_header["CRVAL1"], science_header["CD1_1"], science_header["NAXIS1"]
    )

    # Sets all data and errs within the GMOS chip gaps to 1, so they don't get flagged as bad
    # pixels or trip up the bin definition stage. Chip gaps are identified as pixel columns which
    # are all zeros, and then three columns either side of the chip gaps are also flagged just to
    # be safe.
    zero_rows = [1 if all(x == 0) else 0 for x in data.T]
    boundary_columns = 3
    zero_rows = np.concatenate(
        [np.zeros(boundary_columns), zero_rows, np.zeros(boundary_columns)]
    )
    for n in reversed(range(boundary_columns)):
        zero_rows = [
            1 if x == 0 and zero_rows[y + 1] == 1 else x
            for y, x in enumerate(zero_rows[:-1])
        ]
        zero_rows = [
            1 if x == 0 and zero_rows[y - 1] == 1 else x
            for y, x in enumerate(zero_rows[1:])
        ]
    zero_rows = [1 if x > 0 else x for x in zero_rows]
    chip_gap_map = np.tile(zero_rows, (np.shape(data)[0], 1))
    data[chip_gap_map == 1] = 1.0
    errs[chip_gap_map == 1] = 1.0

    return data, errs, qual, original_qual, header_dict, wavelength_axis


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
