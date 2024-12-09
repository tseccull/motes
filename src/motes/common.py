#!/home/tom/miniforge3/envs/work/bin/python

"""
common.py - Common functions for the MOTES pipeline.
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
