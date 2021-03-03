###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import copy
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import sys
import warnings

from matplotlib.widgets import Slider
from scipy.optimize import least_squares

###############################################################################
# FUNCTIONS //////////////////////////////////////////////////////////////////#
###############################################################################

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Calculate the FWHM and hence the extraction limits from a Moffat profile based on the distance from the central peak
# as a multiple of FWHM. A high resolution profile can also be produced if necessary by setting highres=True.
def extraction_limits(moffparams, width_multiplier, axesdict, highres=False):

    # If requested create a Moffat profile based on the input parameters that has 50 times the resolution of the
    # original data. Else, use the native resolution of the data
    if highres:
        r = np.linspace(0, axesdict['spataxislen'] - 0.02, num=axesdict['spataxislen'] * 50)
        moffprof = moffat(moffparams[0], moffparams[1], moffparams[2], moffparams[3], moffparams[4], moffparams[5], r)
        fwhm = (2 * moffparams[2] * (((2**(1/moffparams[3]))-1)**0.5))*50

    else:
        r = np.linspace(0, axesdict['spataxislen'] - 1, num=axesdict['spataxislen'])
        moffprof = moffat(moffparams[0], moffparams[1], moffparams[2], moffparams[3], 0, 0, r)
        fwhm = 2 * moffparams[2] * (((2**(1/moffparams[3]))-1)**0.5)

    # Locate the peak of the Moffat profile
    indexmax = np.argmax(moffprof)

    # Respectively define the upper and lower extraction limits at a distance above and below the peak of the Moffat
    # profile that equals the FWHM of the Moffat profile times a multiplying factor.
    lower_extraction_limit = indexmax - (width_multiplier * fwhm)
    upper_extraction_limit = indexmax + (width_multiplier * fwhm)

    return lower_extraction_limit, upper_extraction_limit, fwhm, moffparams[1]


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Linearly extrapolate the extraction limits at the ends of the 2D spectrum.
def extrap_extraction_lims(extlims, dispaxislen, shortend, longend):
    short_extrap_grad1 = extrap_grad(extlims[0], [   0,  150, 300])
    short_extrap_grad2 = extrap_grad(extlims[1], [   0,  150, 300])
    long_extrap_grad1  = extrap_grad(extlims[0], [-300, -150,  -1])
    long_extrap_grad2  = extrap_grad(extlims[1], [-300, -150,  -1])

    short_extrap_lim1 = extlims[0][0] - (short_extrap_grad1 * shortend)
    short_extrap_lim2 = extlims[1][0] - (short_extrap_grad2 * shortend)
    long_extrap_lim1  = extlims[0][-1] + (long_extrap_grad1 * (dispaxislen - longend))
    long_extrap_lim2  = extlims[1][-1] + (long_extrap_grad2 * (dispaxislen - longend))

    return short_extrap_lim1, short_extrap_lim2, long_extrap_lim1, long_extrap_lim2


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Given a range of data and limits to define a region of that data, calculate the region's gradient.
def extrap_grad(intextlims, median_lims):
    median_of_y_points_x_to_y = np.median(intextlims[median_lims[0]:median_lims[1]])
    median_of_y_points_y_to_z = np.median(intextlims[median_lims[1]:median_lims[2]])
    median_of_x_points_x_to_y = np.median([median_lims[0], median_lims[1]])
    median_of_x_points_y_to_z = np.median([median_lims[1], median_lims[2]])

    return (median_of_y_points_x_to_y - median_of_y_points_y_to_z)/(median_of_x_points_x_to_y - median_of_x_points_y_to_z)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Define the bins of data over which Moffat profiles will be fitted. Each bin is defined such that when summed it will
# have a given signal to noise. So lower signal to noise regions will have larger bins. This function works from the
# centre of the input 2D spectrum outward to the ends in order to ensure a good start for the binning process.
# Flagged bad pixels will be replaced with appropriate values using a method described by the ESO X-Shooter
# pipeline. Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211
def get_bins(fdict, slow, shigh, dispaxislen, params, sky=False):
    # Take S/N threshold (minSNR) and minimum number of columns per dispersion bin (mincols)
    if params['-SUBTRACT_SKY'] and sky==True:
        minSNR = params['-SKY_SNR_BIN_LIM']
    elif params['-SUBTRACT_SKY'] and sky==False:
        minSNR = params['-SNR_BIN_LIM']
    else:
        minSNR = params['-SNR_BIN_LIM']

    mincols = params['-COL_BIN_LIM']

    # Start x at the central pixel column of the dispersion axis
    x = int(dispaxislen / 2)
    width = 0
    binlocations = []

    sys.stdout.write(' >>> Determining spectrum localisation bins on dispersion axis.\n')
    sys.stdout.write('     User-defined S/N threshold = ' + str(minSNR) + '\n')

    # Start at the centre of the dispersion axis and start binning towards the short wavelength end of the spectrum.
    while x - width > 0:
        snrestimate = 0.

        # If the S/N of the current bin has not yet reached the user-defined threshold (minSNR), add one more pixel
        # column to the bin.
        while snrestimate <= minSNR:
            width += 1

            # Stop the loop once the short wavelength end of the spectrum os reached.
            if x - width < 0:
                width = int(0 + x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin, continue to the next
            # iteration and add another pixel column to the bin.
            shortrows = len(list(filter(lambda x : x <= mincols, np.sum(fdict['qual'][:, int(x - width):int(x)], axis=1))))
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.
            datacol = fdict['data'][:, int(x - width):int(x)]
            errscol = fdict['errs'][:, int(x - width):int(x)]
            bindatacol = np.nansum(datacol, axis=1)
            binerrscol = np.sqrt(np.nansum(np.power(errscol, 2), axis=1))
            signal = np.nansum(bindatacol[int(slow):int(shigh)])
            rssnoise = np.sqrt(np.nansum(np.power(binerrscol[int(slow):int(shigh)], 2)))

            #Estimate the S/N
            snrestimate = signal / rssnoise

        if params['-REPLACE_CRBP']:
            # Replace bad pixels, and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211
		    
            # Here a median spatial profile for the current bin is determined by bootstrapping the good pixels
            # in each spatial pixel row with repeats to estimate a distribution of median values for each value in the
            # median spatial distribution for the bin. The mean of each of these distributions is then taken to be the
            # value of that median spatial pixel. The standard error of each of these distributions becomes the error
            # of the flux.
            meddatacol = []
            errmeddatacol = []
            for j in range(np.shape(fdict['data'])[0]):
                bootsamp = np.random.choice(fdict['data'][j, int(x - width):int(x)], size=(int(width), 100), replace=True)
                meddatacol.append(np.median(bootsamp))
                meanstddist = np.std(bootsamp, axis=0)
                meanstd = np.mean(meanstddist)
                errmeddatacol.append(meanstd / np.sqrt(width - 1))
		    
            meddatacol = np.array(meddatacol)
            errmeddatacol = np.array(errmeddatacol)
		    
            # Each pixel in the bin that was flagged to have be bad pixels or be contaminated by CRs are replaced by
            # the value at the corresponding spatial coordinate in the median spatial profile.
            for i in range(int(width)):
                if 0. in fdict['data'][:, int(x-i)]:
                    nocr = np.where(fdict['data'][:, int(x-i)] != 0.)
		    
                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    datanocr = fdict['data'][:, int(x-i)][nocr]
                    medcolnocr = meddatacol[nocr]
                    scale = np.sum(datanocr[slow:shigh])/np.sum(medcolnocr[slow:shigh])
                    if scale < 0:
                        scale = np.max(datanocr[slow:shigh])/np.max(medcolnocr[slow:shigh])
                    meddatacol *= scale
                    errmeddatacol *= scale
		    
                    # Substittue scaled median spatial profile values into the pixels where CRs or bad pixels have been
                    # identified.
                    multiplier = np.ones(len(meddatacol))
                    multiplier[nocr] *= 0
                    crsubstitute = multiplier * meddatacol
                    fdict['data'][:, int(x-i)] += crsubstitute
                    crerrsubstitute = multiplier * errmeddatacol
                    fdict['errs'][:, int(x - i)] += crerrsubstitute

        binlocations.append([int(x - width), int(x), snrestimate])

        x -= width
        width = 0

    x = int(dispaxislen / 2) + 1

    # Repeat the same process as above, starting at the centre of the dispersion axis, but moving outward toward the
    # longer wavelength end of the 2D spectrum.
    while x + width < dispaxislen:
        snrestimate = 0.

        # If the S/N of the current bin has not yet reached the user-defined threshold (minSNR), add one more pixel
        # column to the bin.
        while snrestimate <= minSNR:
            width += 1

            # Stop the loop once the long wavelength end of the spectrum is reached.
            if x + width > dispaxislen:
                width = int(dispaxislen - x)
                break

            # If there aren't enough good pixels in each spatial column of the current bin, continue to the next
            # iteration and add another pixel column to the bin.
            shortrows = len(list(filter(lambda x: x <= mincols, np.sum(fdict['qual'][:, int(x):int(x + width)], axis=1))))
            if shortrows > 0:
                continue

            # Sum the bin in the dispersion direction and determine the S/N where the signal is the sum of the flux in
            # in the bin and the noise is the root sum square of the errors on the flux in the bin. Errors are taken
            # from the spectrum's error frame.
            bindatacol = np.nansum(fdict['data'][:, int(x):int(x + width)], axis=1)
            binerrscol = np.sqrt(np.nansum(np.power(fdict['errs'][:, int(x):int(x + width)], 2), axis=1))
            signal = np.nansum(bindatacol[int(slow):int(shigh)])
            rssnoise = np.sqrt(np.nansum(np.power(binerrscol[int(slow):int(shigh)], 2)))

            # Estimate the S/N
            snrestimate = signal / rssnoise

        if params['-REPLACE_CRBP']:
            # Replace bad pixels, and cosmic rays if requested by the user. The method employed here
            # is the same as that used by ESO in their X-Shooter data reduction pipeline.
            # Modigliani et al. (2010)  Proc. SPIE, 7737, 28 https://doi.org/10.1117/12.857211
            meddatacol = []
            errmeddatacol = []
		    
            # Here a median spatial profile for the current bin is determined by bootstrapping the good pixels
            # in each spatial pixel row with repeats to estimate a distribution of median values for each value in the
            # median spatial distribution for the bin. The mean of each of these distributions is then taken to be the
            # value of that median spatial pixel. The standard error of each of these distributions becomes the error
            # of the flux.
            for j in range(np.shape(fdict['data'])[0]):
                bootsamp = np.random.choice(fdict['data'][j, int(x):int(x + width)], size=(int(width), 100), replace=True)
                meddatacol.append(np.median(bootsamp))
                meanstddist = np.std(bootsamp, axis=0)
                meanstd = np.mean(meanstddist)
                errmeddatacol.append(meanstd / np.sqrt(width - 1))
		    
            meddatacol = np.array(meddatacol)
            errmeddatacol = np.array(errmeddatacol)
		    
            # Each pixel in the bin that was flagged to have be bad pixels or be contaminated by CRs are replaced by
            # the value at the corresponding spatial coordinate in the median spatial profile.
            for i in range(int(width)):
                if 0. in fdict['data'][:, int(x + i)]:
                    nocr = np.where(fdict['data'][:, int(x + i)] != 0.)
		    
                    # Scale the median spatial profile so its summed flux within the spectrum aperture is equal to the
                    # same for the pixel column being fixed.
                    datanocr = fdict['data'][:, int(x + i)][nocr]
                    medcolnocr = meddatacol[nocr]
                    scale = np.sum(datanocr[slow:shigh]) / np.sum(medcolnocr[slow:shigh])
                    if scale < 0:
                        scale = np.max(datanocr[slow:shigh]) / np.max(medcolnocr[slow:shigh])
                    meddatacol *= scale
                    errmeddatacol *= scale
		    
                    # Substitue scaled median spatial profile values into the pixels where CRs or bad pixels have been
                    # identified.
                    multiplier = np.ones(len(meddatacol))
                    multiplier[nocr] *= 0
                    crsubstitute = multiplier * meddatacol
                    fdict['data'][:, int(x + i)] += crsubstitute
                    crerrsubstitute = multiplier * errmeddatacol
                    fdict['errs'][:, int(x + i)] += crerrsubstitute

        binlocations.append([int(x), int(x + width), snrestimate])

        x += width
        width = 0

    # Sort out the bins into the correct order on the dispersion axis.
    binlocations.sort(key=lambda x: x[0])

    return binlocations, fdict


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Print and plot the output of get_bins
def get_bins_output(binparams, params, lowext, highext, data2D, headparams, axdict):
    sys.stdout.write(' >>> ' + str(len(binparams)) + ' spectrum localisation bins determined on dispersion axis.\n')

    # DIAGNOSTICS - Plot boundary locations of localisation bins on the dispersion axis.
    if params['-DIAG_PLOT_BIN_LOC']:
        drawlines = []
        for bin in binparams:
            binlineloc = np.where(np.logical_and(axdict['saxis'] > lowext - ((highext - lowext) * 0.2), axdict['saxis'] < highext + ((highext - lowext) * 0.2)))
            drawlines.append(np.ones(len(axdict['saxis'][binlineloc])) * bin[0] + axdict['wavstart'])
            drawlines.append(axdict['saxis'][binlineloc] + axdict['imgstart'])
            drawlines.append(np.ones(len(axdict['saxis'][binlineloc])) * bin[1] + axdict['wavstart'])
            drawlines.append(axdict['saxis'][binlineloc] + axdict['imgstart'])

        show_img(data2D, axdict, headparams, drawlines, '2D Spectrum with Boundaries of Localisation Bins')

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Takes an input of extraction limits from the fitting of the binned data and interpolates the limits over the unbinned
# data. Limits are also linearly extrapolated towards the ends of the spectral range.
def interpolate_extraction_lims(extractionlims, dispaxislen, interpkind):
    # If the 2D spectrum was so faint that only 1 dispersion bin could be determined, set the extraction limits
    # across the unbinned spectrum to be a simple linear aperture that has the same extraction limits as was
    # determined for that one bin.
    if len(extractionlims[0]) == 1:
        finalextlims = [np.repeat(extractionlims[1][0], dispaxislen), np.repeat(extractionlims[2][0], dispaxislen)]

    # Otherwise, interpolate the extraction limits form the bins across the unbinned wavelength axis. Also
    # extrapolate the extraction limits in a linear fashion at the ends of the wavelength axis of the 2D spectrum
    # so that the full wavelength axis is covered.
    else:
        interpextract_1 = interp.interp1d(extractionlims[0], extractionlims[1], kind=interpkind)
        interpextract_2 = interp.interp1d(extractionlims[0], extractionlims[2], kind=interpkind)

        intermextlims = [interpextract_1(np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))).astype(int),
                         interpextract_2(np.array(np.arange(extractionlims[0][0], extractionlims[0][-1]))).astype(int)]

        shortextraplim1, shortextraplim2, longextraplim1, longextraplim2 = extrap_extraction_lims(intermextlims, dispaxislen, extractionlims[0][0], extractionlims[0][-1])

        extlim1 = np.insert(intermextlims[0], 0, shortextraplim1)
        extlim1 = np.append(extlim1, longextraplim1)

        extlim2 = np.insert(intermextlims[1], 0, shortextraplim2)
        extlim2 = np.append(extlim2, longextraplim2)

        nextextlims = [extlim1, extlim2]
        nextxaxis = np.insert(np.arange(extractionlims[0][0], extractionlims[0][-1]), 0, 0)
        nextxaxis = np.append(nextxaxis, dispaxislen)

        interpnextlims1 = interp.interp1d(nextxaxis, nextextlims[0], kind='linear', fill_value='extrapolate')
        interpnextlims2 = interp.interp1d(nextxaxis, nextextlims[1], kind='linear', fill_value='extrapolate')

        finalextlims = [interpnextlims1(np.array(range(dispaxislen))), interpnextlims2(np.array(range(dispaxislen)))]

    return finalextlims


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Returns a wavelength axis array using a start wavelength, the wavelength increment and the number of values along the
# axis required.
def make_wav_axis(start, increment, length):
    return (increment * np.array(range(length))) + start


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function detects cosmic rays in the 2D image and masks the pixels that contain them, by separating the spectrum
# region from the two background regions. Each region is sigma-clipped separately due to the different distribution of
# fluxes present in the background and the vicinity of the 2D spectrum.
def mask_cosmic_rays(data, lext, hext, multiplier=4.):
    # Divide 2D spectrum into spectrum and background sections
    sections = [data[:int(lext), :], data[int(lext):int(hext), :], data[int(hext):, :]]
    masks = []

    # Sequentially sigma-clip each section, and create a mask marking the locations of pixels to be kept.
    for sect in sections:

        datamed = np.median(sect)
        datastd = np.std(sect)

        upperlim = datamed + (multiplier * datastd)
        lowerlim = datamed - (multiplier * datastd)

        sectmask2D = (sect < upperlim) & (sect > lowerlim)
        masks.append(sectmask2D)

    # Assemble masks for each section into a full mask for the input 2D spectrum.
    mask2D = np.vstack((masks[0], masks[1], masks[2]))

    return mask2D


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Creates moffat profile added a linear sloped background based on input parameters
def moffat(amp, c, alpha, beta, bglevel, bggrad, datarange):
    moff = amp * ((1 + (((datarange - c) * (datarange - c)) / (alpha * alpha))) ** -beta)

    return moff + bglevel + (datarange * bggrad)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to the column using a
# Levenberg-Marquardt least squares method. Returns the best fit parameters of the Moffat function
def moffat_least_squares(r, col, seeing, pixres):
    # Set up initial conditions for the least squares fit.
    # amplitude, centre, alpha, beta, background gradient, background level
    # alpha value calculated by substituting seeing (FWHM) converted to pixels and beta = 4.765 into the equation:
    #               alpha = FWHM / (2 * ((2 ^ (1 / beta)) - 1) ^ 0.5)
    # Initial beta value comes from optimal value from atmospheric turbulence theory as described in
    # Trujillo I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [np.nanmedian(np.sort(col)[-3:]), np.argmax(col), seeing / 0.7914, 4.765, 0.,
          np.median(np.concatenate((col[:5], col[-5:])))]

    # Run the least squares fit.
    res_lsq = least_squares(moffat_resid, 
                            x0, 
                            bounds=([0., np.argmax(col)-1., 0., 0., 0., -np.inf], [np.inf, np.argmax(col)+1, (5*seeing/0.7914), 5., np.inf, np.inf]), 
                            args=(r, col), 
                            method='trf', 
                            ftol=1E-12
                            )

    return [res_lsq.x[0], res_lsq.x[1], res_lsq.x[2], res_lsq.x[3], res_lsq.x[4], res_lsq.x[5]]


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Calculates residuals of fitted moffat profile and the data for the Levenberg Marquardt least squares method.
# A = x[0]
# c = x[1]
# alpha = x[2]
# beta = x[3]
# B = x[4]
# m = x[5]
#@jit
def moffat_resid(x, datarange, data):
    moff = x[0] * ((1 + ((datarange - x[1]) * (datarange - x[1])) / (x[2] * x[2])) ** -x[3])

    return moff + x[4] + (datarange * x[5]) - data


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Plot the spatial profile of a collapsed spectrum or a collapsed bin therein, and plot the Moffat function fitted to
# the data on top.
def plot_fitted_spatial_profile(spataxis, bindata, hiresspataxis, binmoffparams, imgstart, headparams):
    plt.figure(figsize=(7,4.5))
    plt.plot(hiresspataxis + imgstart, 
             moffat(binmoffparams[0], binmoffparams[1], binmoffparams[2], binmoffparams[3], binmoffparams[4], binmoffparams[5], hiresspataxis), 
             color='r', 
             linewidth=3, 
             label='Fitted Moffat Profile'
             )
    
    plt.plot(spataxis + imgstart, bindata, color='k', label='Spatial Profile')
    plt.grid(linestyle='dashed', color='gray')
    plt.legend()
    plt.title('Spectrum Spatial Profile and Fitted Moffat Profile')
    plt.xlabel('Spatial Axis, Pixels')
    plt.ylabel('Median Flux, ' + headparams['fluxunit'])
    plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# This function prints, to the terminal, a set parameters describing a fitted Moffat function.
def printmoffparams(moffparams, imgstart, datascale):
    sys.stdout.write(' >>> Fitted Moffat function parameters:\n')
    sys.stdout.write('         A = ' + str(moffparams[0]) + '\n')
    sys.stdout.write('         c = ' + str(moffparams[1] + imgstart) + '\n')
    sys.stdout.write('     alpha = ' + str(moffparams[2]) + '\n')
    sys.stdout.write('      beta = ' + str(moffparams[3]) + '\n')
    sys.stdout.write('         B = ' + str(moffparams[4]) + '\n')
    sys.stdout.write('         m = ' + str(moffparams[5]) + '\n\n')
    sys.stdout.write(' >>> Profile scaling factor used for fitting: ' + str(datascale) + '\n')
    sys.stdout.write(' >>> Plot of median spatial profile presents the orginal unscaled profile.\n')

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Takes an input image and line data to be drawn on that image and creates a figure to be shown on screen.
def show_img(data2D, axdict, headparams, drawlines, title):

    # Catch and suppress the following UserWarning:
    # UserWarning: Warning: 'partition' will ignore the 'mask' of the MaskedArray.
    #     a.partition(kth, axis=axis, kind=kind, order=order)
    # Numpy masked array used only when displaying the spectrograms.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        warnings.warn("partition", UserWarning)
        
        power = int(np.floor(np.log10(np.mean(data2D))))
        data2D = copy.deepcopy(data2D)/10**power

        figwidth = 10.
        fig = plt.figure(figsize=(figwidth, figwidth / 1.9))
        gs = gridspec.GridSpec(18, 33)
        ax = plt.subplot(gs[:, :32])
        colax = plt.subplot(gs[1:, 32])
        masked_data2D = np.ma.masked_where(data2D == 0, data2D)
        cmap = matplotlib.cm.afmhot
        cmap.set_bad(color='red')
        s = ax.imshow(masked_data2D, 
                      aspect='auto', 
                      vmin=0, 
                      vmax=np.nanmedian(masked_data2D) + (0.5 * np.nanstd(masked_data2D)), 
                      origin='lower', 
                      cmap=cmap,
                      extent=[axdict['wavstart'], axdict['wavstart'] + len(axdict['waxis']), axdict['saxis'][0] + axdict['imgstart'], axdict['saxis'][-1] + axdict['imgstart']]
                      )
        
        for i in range(int(len(drawlines)/2)):
            ax.plot(drawlines[i*2], drawlines[(i*2)+1], color='white')
        cbar = fig.colorbar(s, cax=colax)
        cbar.ax.yaxis.set_offset_position('left')
        cbar.ax.set_ylabel('Pixel Flux, x10^' + str(power) + ' ' + headparams['fluxunit'])
        ax2 = ax.twiny()
        ax2.plot(axdict['waxis'], data2D[0, :], alpha=0)
        ax2.set_xlim(axdict['waxis'][0], axdict['waxis'][-1])
        ax2.set_xlabel('Wavelength, ' + headparams['wavunit'])
        ax.set_ylim(axdict['saxis'][0] + axdict['imgstart'], axdict['saxis'][-1] + axdict['imgstart'])
        ax.set_ylabel('Spatial Axis, Pixels')
        ax.set_xlim(axdict['wavstart'], axdict['wavstart'] + len(axdict['waxis']))

        ax.set_xlabel('Dispersion Axis, Pixels')

        plt.title(title, y=1.095)
       
        # Add interactive scaling bar to figures if the 2D spectrum isn't flux calibrated. 
        # Again for some reason teeny tiny numbers cause things to break.
        fig.subplots_adjust(bottom=0.2)
        axvmin = plt.axes([0.1, 0.05, 0.8, 0.03])
        axvmax = plt.axes([0.1, 0.01, 0.8, 0.03])
        smin = Slider(axvmin, 'LowCut', 0, np.nanmax(masked_data2D)-1, valinit=0., valstep=0.001*np.nanmax(masked_data2D))
        smax = Slider(axvmax, 'HighCut', 1., np.nanmax(masked_data2D), valinit=np.nanmedian(masked_data2D) + (0.5 * np.nanstd(masked_data2D)), valstep=0.001*np.nanmax(masked_data2D))
	    
        def update(val):
            vmax = smax.val
            smin.valmax = vmax-1
            vmin = smin.val
            smax.valmin = vmin+1
            
            s.set_clim(vmin, vmax)
            cbar.ax.set_ylabel('Pixel Flux, x10^' + str(power) + ' ' + headparams['fluxunit'])
            fig.canvas.draw_idle()
            
        smin.on_changed(update)
        smax.on_changed(update)

        plt.show()

    return None


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Subtracts the sky background from the 2D image by defining bg regions using limits input to the function and then
# fitting a profile to the background column by column while masking cosmic rays. The fitted linear for
# each column is subtracted from the full column to produce a background subtracted 2D image.
def subtract_sky(bglowext, bghighext, fdict, axdict, crmask_multiplier=1.):

    fdict['data'] = fdict['data'].T
    fdict['errs'] = fdict['errs'].T

    medsky = []

    for ii in range(len(bglowext)):
        if bglowext[ii] < 0:
            bglowext[ii] = 0
        if bghighext[ii] > axdict['saxis'][-1]:
            bghighext[ii] = axdict['saxis'][-1]

        datacol = fdict['data'][ii]
        colrange = range(len(datacol))
        skypix = datacol[np.where(np.logical_or(colrange<bglowext[ii], colrange>bghighext[ii]))]
        bootsky = np.random.choice(skypix, (len(skypix), 100), replace=True)
        skysamp = np.median(bootsky, axis=0)
        skylevel = np.mean(skysamp)
        medsky.append(skylevel)
        skyerr = np.std(skysamp)/(99**0.5)
        fdict['data'][ii] -= skylevel
        fdict['errs'][ii] = ((fdict['errs'][ii]**2)+(skyerr**2))**0.5

    medsky = np.array(medsky)
    fdict['skymod'] = np.tile(medsky, (np.shape(fdict['data'])[1], 1))
    fdict['data'] = fdict['data'].T
    fdict['errs'] = fdict['errs'].T
    
    chipgaps = np.where(np.median(fdict['data'], axis=0)==0)
    fdict['data'][:, chipgaps[0]] += 1

    return fdict
