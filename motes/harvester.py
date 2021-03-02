###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import astropy.io.fits as fits
import copy
import datetime
import glob
import motes.common as common
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

from astropy.table import Table

###############################################################################
# FUNCTIONS ///////////////////////////////////////////////////////////////// #
###############################################################################
# ////////////////////////////////////////////////////////////////////////////#
# EXTRACT HEADER DATA AND DATAFRAMES, AND REPACKAGE THEM INTO DICTIONARIES

def data_harvest(intreg, filename_2D, chunk, param_dict):
    # Create dictionary to tell data_harvest which instrument specific
    # function to call.
    instrument_dict = {   'FORS2': harvest_fors2, 
                         'GMOS-N': harvest_gmos,
                         'GMOS-S': harvest_gmos, 
                       'XSHOOTER': harvest_xshoo
                       }

    # Open file containing the spectral data, open it up and extract the
    # header, image frame, error frame, and quality frame (if it has one).
    with fits.open(filename_2D) as imgfile:
        imghead = imgfile[0].header
        inst = imghead['INSTRUME']
        # Based on the value of inst, this calls one of the harvest_instrument functions.
        imgdata, imgerrs, imgqual, ogimgqual, head_dict, wavaxis = instrument_dict[inst](imgfile, imghead)

    # SLICE ALL DATAFRAMES BASED ON THE INPUT FROM reg.txt
    imgshape = np.shape(imgdata)
    imgstart = int(0 + chunk[intreg][0])
    imgend = int(imgshape[0] - chunk[intreg][1])

    # Slice off the spatial rows outside the spatial region.
    datasliced = imgdata[imgstart:imgend + 1, :]
    errssliced = imgerrs[imgstart:imgend + 1, :]
    qualsliced = imgqual[imgstart:imgend + 1, :]
    dataslicedshape = np.shape(datasliced)
    sys.stdout.write(' >>> 2D spectrum sliced on spatial axis based on user defined limits:\n'
                     '     New spatial axis covers pixel rows ' + str(imgstart) + '-' + str(imgend) + '.\n')

    # Create spatial axis for the 2D spectrum and a high resolution version (standard res * 5) for the purposes of plotting
    spataxis = np.linspace(0., float(dataslicedshape[0] - 1), num=dataslicedshape[0])
    hiresspataxis = np.linspace(spataxis[0], spataxis[-1], num=len(spataxis) * 5)

    if chunk[intreg][2]<wavaxis[0] or chunk[intreg][3]>wavaxis[-1]:
        sys.stdout.write(' >>> User defined wavelength limit(s) are outside native wavelength range\n'
                         '     Make sure "-LOW_WAV_SLICE" > lower limit of wavelength axis\n'
                         '     Make sure "-HIGH_WAV_SLICE" < upper limit of wavelength axis\n'
                         '     Terminating MOTES.\n\n')
        exit()

    wavslice = np.where(np.logical_and(wavaxis >= chunk[intreg][2], wavaxis <= chunk[intreg][3]))
    wavstart = wavslice[0][0]
    wavend = wavslice[0][-1]
    wavaxis = wavaxis[wavslice]

    datasliced = np.squeeze(datasliced[:, wavslice])
    errssliced = np.squeeze(errssliced[:, wavslice])
    qualsliced = np.squeeze(qualsliced[:, wavslice])

    sys.stdout.write(' >>> 2D spectrum sliced on dispersion axis based on user defined limits:\n'
                     '     New Wavelength range is ' + str(chunk[intreg][2]) + '-' + str(chunk[intreg][3]) + ' ' + head_dict['wavunit'] + '.\n'
                     '     This range is equivalent to pixel columns ' + str(wavstart) + '-' + str(wavstart + len(wavaxis)) + '\n')

    frame_dict = {'data': datasliced, 
                  'errs': errssliced, 
                  'qual': qualsliced,
                'ogqual': ogimgqual}

    axes_dict = {'spataxislen': len(spataxis),
                       'saxis': spataxis, 
                     'hrsaxis': hiresspataxis, 
                    'imgstart': imgstart,
                      'imgend': imgend, 
                 'dispaxislen': len(wavaxis), 
                       'waxis': wavaxis,
                    'wavstart': wavstart,
                      'wavend': wavend
                 }

    return head_dict, frame_dict, axes_dict, imghead


# ////////////////////////////////////////////////////////////////////////////#
# HARVEST THE HEADER AND DATA FROM A FORS2 SPECTRUM
def harvest_fors2(imgfilehdu, imgheader):
    # Retrieve the data frame and error frame.
    imgdata = imgfilehdu[0].data
    imgerrs = imgfilehdu[1].data**0.5
    imgqual = np.ones(np.shape(imgdata))
    ogimgqual = copy.deepcopy(imgqual)-1

    # Determine the spatial pixel resolution of the image in arcsec.
    if imgheader['HIERARCH ESO DET WIN1 BINY'] == 1 and imgheader['HIERARCH ESO INS COLL NAME'] == 'COLL_SR':
        pixres = 0.125
    elif imgheader['HIERARCH ESO DET WIN1 BINY'] == 1 and imgheader['HIERARCH ESO INS COLL NAME'] == 'COLL_HR':
        pixres = 0.0632
    elif imgheader['HIERARCH ESO DET WIN1 BINY'] == 2 and imgheader['HIERARCH ESO INS COLL NAME'] == 'COLL_SR':
        pixres = 0.25
    elif imgheader['HIERARCH ESO DET WIN1 BINY'] == 2 and imgheader['HIERARCH ESO INS COLL NAME'] == 'COLL_HR':
        pixres = 0.125
    else:
        sys.stdout.write('FAILED.\n')
        sys.stdout.write('     Non-standard binning used in image.\n'
                         '     Spatial pixel resolution could not be determined.\n')
        sys.stdout.write('     Extraction of ' + direc + ' could not be completed.\n')
        sys.stdout.write('     Terminating MOTES.\n\n')
        exit()

    sys.stdout.write(' >>> Spatial pixel resolution determined: ' + str(pixres) + '"\n')

    # Put header information into a dictionary
    sys.stdout.write(' >>> Gathering required information from FITS header. ')
    sys.stdout.flush()
    headerdict = {       'object': imgheader['OBJECT'].replace(' ', '_'), 
                  'pixresolution': pixres,
                        'exptime': imgheader['HIERARCH ESO INS SHUT EXPTIME'], 
                           'inst': imgheader['INSTRUME'],
                         'seeing': 0.5 * (imgheader['HIERARCH ESO TEL AMBI FWHM START'] + imgheader['HIERARCH ESO TEL AMBI FWHM END']),
                       'fluxunit': imgheader['BUNIT'], 
                        'wavunit': 'Angstrom'
                   }
    
    sys.stdout.write('DONE.\n')

    # SLICE IMAGE ON WAVELENGTH AXIS #
    # Create the wavelength axis of the spectrum, define the region of the spectrum in dispersion space to be
    # kept, and then slice the 2D image on wavelength axis accordingly.
    wavaxis = common.make_wav_axis(imgheader['CRVAL1'], imgheader['CD1_1'], imgheader['NAXIS1'])

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis


# ////////////////////////////////////////////////////////////////////////////#
# HARVEST THE HEADER AND DATA FROM A GMOS SPECTRUM
def harvest_gmos(imgfilehdu, imgheader):
    # Retrieve the data frame, error frame, and qual frame.
    scihead = imgfilehdu['SCI'].header
    imgdata = imgfilehdu['SCI'].data
    imgerrs = imgfilehdu['VAR'].data**0.5
    imgqual = imgfilehdu['DQ'].data
    ogimgqual = copy.deepcopy(imgqual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    # Pixels in chip gaps are kept as 1 to make sure they don't get flagged as CRs later on.
    imgqual[np.where(imgdata+imgqual==1)] = 0
    imgqual[imgqual==0] = 1

    # All this is to get an initial estimate of the IQ. Tables below are based on the condition constraints used by Gemini.
    # See https://www.gemini.edu/observing/telescopes-and-sites/sites#ImageQuality
    IQ_dict = { '20-percentile': 0, 
                '70-percentile': 1, 
                '85-percentile': 2, 
               '100-percentile': 3
               }
    
    WavTab = np.array([[0000.,4000.,0],
                       [4000.,5500.,1],
                       [5500.,7000.,2],
                       [7000.,8500.,3],
                       [8500.,9750.,4],
                       [9750.,11000.,5]])
    
    IQTab = np.array([[0.6, 0.90, 1.20, 2.00],
                      [0.6, 0.85, 1.10, 1.90],
                      [0.5, 0.75, 1.05, 1.80],
                      [0.5, 0.75, 1.05, 1.70],
                      [0.5, 0.70, 0.95, 1.70],
                      [0.4, 0.70, 0.95, 1.65]])

    iq = imgheader['RAWIQ']
    wav_min = scihead['CRVAL1']

    for i in WavTab:
        if scihead['CRVAL1'] > i[0] and scihead['CRVAL1'] < i[1]:
            seeing = IQTab[int(i[2])][int(IQ_dict[iq])]
            break

    # Put header information into a dictionary
    sys.stdout.write(' >>> Gathering required information from FITS header. ')
    sys.stdout.flush()
    headerdict = {       'object': imgheader['OBJECT'].replace(' ', '_'), 
                  'pixresolution': imgheader['PIXSCALE'],
                        'exptime': imgheader['EXPTIME'], 
                         'seeing': seeing,
                           'inst': imgheader['INSTRUME'], 
                        'wavunit': scihead['WAT1_001'].split(' ')[2].split('=')[1]
                  }
    
    # BUNIT only appears in the headers of GMOS spectra if they have been flux calibrated.
    if 'BUNIT' in scihead:
        headerdict['fluxunit'] = scihead['BUNIT']
    else:
        headerdict['fluxunit'] = 'electrons'
                  
    sys.stdout.write('DONE.\n')
    sys.stdout.write(' >>> Spatial pixel resolution determined: ' + str(headerdict['pixresolution']) + '"\n')

    # SLICE IMAGE ON WAVELENGTH AXIS #
    # Create the wavelength axis of the spectrum, define the region of the spectrum in dispersion space to be
    # kept, and then slice the 2D image on wavelength axis accordingly.
    wavaxis = common.make_wav_axis(scihead['CRVAL1'], scihead['CD1_1'], scihead['NAXIS1'])

    # Sets all values within the GMOS chip gaps to NaN, so they don't get flagged as bad pixels later on during CR masking.
    # This will be reversed later following cr_handling()
    loc = np.where(np.where(np.median(imgdata, axis=0)==0)[0]-np.roll(np.where(np.median(imgdata, axis=0)==0)[0], 1)!=1)
    locs = [0, loc[0][1]-1, loc[0][1], -1]
    
    vals = np.array([np.where(np.median(imgdata, axis=0)==0)[0][locs[0]]-2,
                     np.where(np.median(imgdata, axis=0)==0)[0][locs[1]]+2,
                     np.where(np.median(imgdata, axis=0)==0)[0][locs[2]]-2,
                     np.where(np.median(imgdata, axis=0)==0)[0][locs[3]]+2]
                     )
    imgdata[:,vals[0]:vals[1]] = 1
    imgdata[:,vals[2]:vals[3]] = 1
    
    # Set errors in the chip gaps to 1.0 so they don't trip up the Moffat fitting.
    imgerrs[:,vals[0]:vals[1]] = 1.
    imgerrs[:,vals[2]:vals[3]] = 1.
    
    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis


# ////////////////////////////////////////////////////////////////////////////#
# HARVEST THE HEADER AND DATA FROM A X-SHOOTER SPECTRUM
def harvest_xshoo(imgfilehdu, imgheader):
    # Retrieve the data frame, error frame, and qual frame.
    imgdata = imgfilehdu[0].data
    imgerrs = imgfilehdu[1].data
    imgqual = imgfilehdu[2].data
    ogimgqual = copy.deepcopy(imgqual)

    # Convert qual frame to boolean. Good pixels = 1; Bad pixels = 0
    imgqual[imgqual>0] *= -1
    imgqual[imgqual==0] = 1
    imgqual[imgqual<0] = 0

    # Put header information into a dictionary
    sys.stdout.write(' >>> Gathering required information from FITS header. ')
    sys.stdout.flush()
    headerdict = {       'object': imgheader['OBJECT'].replace(' ', '_'), 
                  'pixresolution': imgheader['CDELT2'],
                        'exptime': imgheader['EXPTIME'],
                         'seeing': 0.5 * (imgheader['HIERARCH ESO TEL AMBI FWHM START'] + imgheader['HIERARCH ESO TEL AMBI FWHM END']),
                       'fluxunit': imgheader['BUNIT'], 
                        'wavunit': 'nm'
                  }
    sys.stdout.write('DONE.\n')
    sys.stdout.write(' >>> Spatial pixel resolution determined: ' + str(headerdict['pixresolution']) + '"\n')

    # SLICE IMAGE ON WAVELENGTH AXIS #
    # Create the wavelength axis of the spectrum, define the region of the spectrum in dispersion space to be
    # kept, and then slice the 2D image on wavelength axis accordingly.
    wavaxis = common.make_wav_axis(imgheader['CRVAL1'], imgheader['CDELT1'], imgheader['NAXIS1'])

    return imgdata, imgerrs, imgqual, ogimgqual, headerdict, wavaxis
