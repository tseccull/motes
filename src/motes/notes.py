#!/home/tom/miniforge3/envs/work/bin/python

"""
messages.py - Contains functions to print all stdout messages printed to
              terminal by MOTES during operation. Messages are presented
              here in no particular order.
"""


import numpy as np
import sys


def data_harvest_1(floor, ceiling):
	sys.stdout.write(
        " >>> 2D spectrum sliced on spatial axis based on user defined limits:\n"
        "     New spatial axis covers pixel rows "
        + str(floor)
        + "-"
        + str(ceiling)
        + ".\n"
    )

    
def data_harvest_2():
	sys.stdout.write(
        " >>> User defined wavelength limit(s) are outside native wavelength range\n"
        '     Make sure "-LOW_WAV_SLICE" > lower limit of wavelength axis\n'
        '     Make sure "-HIGH_WAV_SLICE" < upper limit of wavelength axis\n'
        "     Terminating MOTES.\n\n"
    )


def data_harvest_3(
    data_regions, region_counter, wavelength_unit, wavelength_start, wavelength_axis_length
):
	sys.stdout.write(
        " >>> 2D spectrum sliced on dispersion axis based on user defined limits:\n"
        "     New Wavelength range is "
        + str(data_regions[region_counter][2])
        + "-"
        + str(data_regions[region_counter][3])
        + " "
        + wavelength_unit
        + ".\n"
        "     This range is equivalent to pixel columns "
        + str(wavelength_start)
        + "-"
        + str(wavelength_start + wavelength_axis_length)
        + "\n"
    )


def done():
    sys.stdout.write("DONE.\n")


def harvest_gathering():
    sys.stdout.write(" >>> Gathering required information from FITS header. ")
    sys.stdout.flush()


def harvest_pixel_resolution(pixel_resolution):
	sys.stdout.write(
        " >>> Spatial pixel resolution determined: "
        + str(pixel_resolution)
        + '"\n'
    )


def read_motes_param_file_1():
    sys.stdout.write(" >>> Reading in parameters from motesparams.txt. ")
    sys.stdout.flush()


def read_regions_1():
    sys.stdout.write(" >>> reg.txt file found. Reading reg.txt. ")
    sys.stdout.flush()
    

def read_regions_2(working_directory):
    sys.stdout.write(" >>> reg.txt file not found in root working directory.\n")
    sys.stdout.write("     Root working directory is:\n")
    sys.stdout.write("     " + working_directory + "\n")
    sys.stdout.write(
        "     Please provide a reg.txt file to tell MOTES the spatial extent of the 2D\n"
    )
    sys.stdout.write("     spectra. See the docs for further information.\n")
    sys.stdout.write("     Terminating MOTES.\n\n")


def save_fits_1(file_path):
    sys.stdout.write(" >>> Spectrum extracted and saved:\n")
    sys.stdout.write(
        "     "
        + "/".join(file_path.split("/inputs")[0:-1])
        + "/m" + file_path.split("/")[-1]
        + "\n"
    )

def sky_locator_1():
    sys.stdout.write(" >>> Fitting Moffat Functions to each bin to localise 2D spectrum.\n")


def sky_locator_2():
    sys.stdout.write("     Fitting complete.\n")
    sys.stdout.write(" >>> Drawing target/sky boundaries. ")
    sys.stdout.flush()
    
   
def sky_locator_3():
    sys.stdout.write(" >>> Subtracting sky.\n")
    sys.stdout.flush()
    
    
def sky_locator_4():
    sys.stdout.write("\n")
    sys.stdout.flush()
    
    
def subtract_sky_1(data_shape, seeing):
    bg_fwhm_multiplier_limit = round(np.min(data_shape) / (2 * seeing), 1)
    sys.stdout.write(" >>> No pixels contained inside sky region.\n")
    sys.stdout.write("     -BG_FWHM_MULTIPLIER in motesparams.txt is probably too large.\n")
    sys.stdout.write("     Please reduce -BG_FWHM_MULTIPLIER and try again.\n")
    sys.stdout.write(" >>> No pixels contained inside sky region.\n")
    sys.stdout.write("     -BG_FWHM_MULTIPLIER in motesparams.txt is probably too large.\n")
    sys.stdout.write("     Please reduce -BG_FWHM_MULTIPLIER and try again.\n")
    sys.stdout.write(
        "     -BG_FWHM_MULTIPLIER < "
        + str(bg_fwhm_multiplier_limit)
        + " recommended in this case.\n"
    )
    sys.stdout.write(
        "     Enlarging the 2D spectrum region in reg.txt is also a viable solution.\n"
    )
    sys.stdout.write("     Terminating MOTES.\n")
    
    
def subtract_sky_2(ii, number_of_columns):
    sys.stdout.write(
        "     " 
        + str(ii + 1)
        + "/"
        + str(number_of_columns)
        + " columns completed.\r"
    )
