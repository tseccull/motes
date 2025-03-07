#!/usr/bin/env python3

"""
	common.py

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
	`common.py` contains common functions shared by at least two other 
	modules in the MOTES package.
	
"""

import numpy as np


def filter_data(data_2d, errs_2d):
    """
    This function takes in the data_2d and errs_2d and outputs frames
    where any NaN or Inf values are 0.0 (this can be later filtered to a
    median of the column). This is used before the extraction procedures
    to ensure the S/N in the bin is numerical (required for optimal
    extraction).

    Args:
     -- data_2d (numpy.ndarray)
          Original data_2d
     -- errs_2d (numpy.ndarray)
          Original errs_2d

    Returns:
     -- data_2d (numpy.ndarray)
          Filtered data_2d
     -- errs_2d (numpy.ndarray)
          Filtered errs_2d
    """

    errs_2d[~np.isfinite(errs_2d)] = 0.0
    data_2d[~np.isfinite(errs_2d)] = 0.0
    data_2d[~np.isfinite(data_2d)] = 0.0
    errs_2d[~np.isfinite(data_2d)] = 0.0

    return data_2d, errs_2d


def make_wav_axis(start, increment, length):
    """
    Returns a wavelength axis array using a start wavelength, the
    wavelength increment and the number of values required along the
    axis.

    Args:
     -- start (float)
          The start wavelength of the axis.
     -- increment (float)
          The wavelength increment between each value along the axis.
     -- length (int)
          Number of values along the axis required.

    Returns:
     -- wavaxis (numpy.ndarray)
          A wavelength axis array.
    """
    
    end = (length * increment) + start
    wavaxis = np.arange(start=start, step=increment, stop=end)
    
    return wavaxis


def moffat(parameters, data_range):
    """
    Creates moffat profile added to a linear sloped background based on
    input parameters.

    Args:
     -- parameters (list)
          List containing the parameters of the Moffat profile.
     -- data_range (numpy.ndarray)
          x axis of the Moffat profile.

    Returns:
     -- moffat_plus_background (numpy.ndarray)
          A moffat profile defined on the data_range axis using the
          parameters input to the function, with added background flux.
    """

    amplitude = parameters[0]
    center = parameters[1]
    alpha = parameters[2]
    beta = parameters[3]
    bg_level = parameters[4]
    bg_grad = parameters[5]

    centered_range = data_range - center    
    unit_moffat = (1 + ((centered_range*centered_range) / (alpha*alpha))) ** -beta
    scaled_moffat = amplitude * unit_moffat        
    moffat_plus_background = scaled_moffat + bg_level + (data_range * bg_grad)

    return moffat_plus_background


def moffat_fwhm(alpha, beta):
    """
    Calculate the fwhm of a Moffat function based on its alpha and beta.
    
    Args:
     -- alpha (float)
          alpha value of the Moffat function.
     -- beta (float)
          beta value of the Moffat function.
    
    Returns:
     -- fwhm (float)
          Full Width at Half Maximum of the Moffat function.
    """
    
    fwhm = 2 * alpha * np.sqrt((2**(1/beta)) - 1)
    
    return fwhm 
