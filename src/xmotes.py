#!/usr/bin/env python3

"""
MOTES Modular and Optimal Tracer and Extractor of Spectra.

Description: Modular and Optimal Tracer and Extractor of Specrtra 
(MOTES). A Python package for extracting spectrum from astronomical
2D spectrograms.

Version: 1.0.0
Date: 2025-02-28
Authors: Tom Seccull, Dominik Kiersz
Licence: GNU General Public License v3.0

This file is the MOTES executable that can be called from the command
line.
"""


import argparse
import motes.core as core


parser = argparse.ArgumentParser(
	description="MOTES is the Modular and Optimized Tracer and\
	Extractor of Spectra. It is designed to locate the traces of\
	spectroscopically observed astronomical point sources in 2D\
	spectroscopic imaging data and then optimally extract those\
	traces to form 1D spectra. Optionally, sky subtraction can be\
	performed via subtraction of polynomial fits to the spatial\
	profiles pixels located outside the extraction aperture."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="[boolean] When set, spectra extracted by MOTES will be\
	saved to new fits files along with relevant metadata stored\
	in the header."
)
parser.add_argument(
    "-m", "--minimum_column_limit", type=int, default=15,
    help="[int; default=15] Defines the minimum number of good\
    pixels allowed in each spatial row of a dispersion bin when\
    defining the localisation bins on the dispersion axis with the\
    getbins function. If this threshold is not met getbins will\
    continuing adding dispersion pixel columns to the bin until it\
    is. Must be positive."
)
parser.add_argument(
    "-xs", "--extraction_snr_limit", type=int, default=10,
    help="[int; default=10] Sets the SNR Limit at which the code\
    stops adding pixel columns to the current bin. SNR is\
    calculated as the summed flux divided by the root sum square of\
    the error values associated with that flux. Must be positive."
)
parser.add_argument(
    "-xf", "--extraction_fwhm_multiplier", type=float, default=2.0,
    help="[float; default=2.0] Sets the multiple of the bin's\
    moffat FWHM along the spatial axis. This is used to define the\
    distance from the centre of the moffat profile at which the\
    extraction limits should be set. A value in the range 2.0-3.0\
    typically gives the best results. Must be positive."
)
parser.add_argument(
	"-k", "--subtract_sky", action="store_true",
	help="[boolean] When set, MOTES will attempt a background/sky\
	subtraction before extracting the spectrum"
)
parser.add_argument(
    "-ko", "--sky_order", type=int, default=0,
    help="[int; default=0] Sets the order of the polynomial fitting\
    to the background regions of the 2D spectrum when performing\
    the sky subtraction. If set to 0, a median sky subtraction is\
    performed. Value must be >= 0."
)
parser.add_argument(
    "-ks", "--sky_snr_limit", type=int, default=100,
    help="[int; default=100] Sets the SNR Limit at which the code\
    stops adding pixel columns to the current bin when defining the\
    binning for location of the background sky regions. SNR is\
    calculated as the summed flux divided by the root sum square of\
    the error values associated with that flux. Must be positive."
)
parser.add_argument(
    "-kf", "--sky_fwhm_multiplier", type=float, default=3.0,
    help="[float; default=3.0] Sets the multiple of the bin's\
    moffat FWHM along the spatial axis. This is used to define the\
    distance from the centre of the moffat profile at which the\
    sky limits should be set. A value >= 3.0 typically gives the\
    best results. Must be positive."
)
parser.add_argument(
	"-ds", "--diag_plot_spectrum", action="store_true",
	help="[boolean] When set, spectra extracted by MOTES will be\
	plotted on screen with matplotlib."
)
parser.add_argument(
	"-dc", "--diag_plot_collapsed", action="store_true",
	help="[boolean] When set, MOTES will make a plot of the median\
	collapsed spatial profile of the 2D spectrum data and the\
	Moffat profile that fitted to it. This plot will be saved\
	with a file name with the form m[input_file_name]_moffat_median.png"
)
parser.add_argument(
	"-dm", "--diag_plot_moffat", action="store_true",
	help="[boolean] When set, the median spatial profile of each\
	dispersion bin will be plotted against the Moffat profile\
	fitted to it. Each plot will be saved as a .png file."
)
parser.add_argument(
	"-db", "--diag_plot_bins", action="store_true",
	help="[boolean] When set, MOTES will plot the locations of bins\
	determined by the get_bins() function on screen."
)
parser.add_argument(
	"-dl", "--diag_plot_localisation", action="store_true",
	help="[boolean] When set, MOTES will plot the extraction limits\
	or sky boundary limits over the 2D data for comparison on\
	screen. The limits/boundaries will be plotted both before and\
	after calculation of the extrapolated ends of the the\
	limits/boundaries."
)
parser.add_argument(
	"-df", "--diag_plot_skyfit", action="store_true",
	help="[boolean] When set, the spatial profile of each\
	background (sky) data column will be plotted against the\
	polynomial fitted to it. Each plot will be saved as a .png file."
)
parser.add_argument(
	"-dk", "--diag_plot_skysub", action="store_true",
	help="[boolean] When set, the 2D sky-subtracted data frame will\
	be plotted on screen."
)
parser.add_argument(
	"-v", "--verbose", action="store_true",
	help="[boolean] When set logging information will be printed\
	in the console. If not set, only warnings or other critical\
	messages will be printed."
)

args = parser.parse_args()
core.motes(args)
