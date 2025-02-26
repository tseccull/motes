#!/home/tom/miniforge3/envs/work/bin/python

import copy
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import motes.common as common
import numpy as np

from matplotlib.widgets import Slider

def get_bins_output(
    bin_parameters,
    parameters,
    spatial_lo_limit,
    spatial_hi_limit,
    data_2d,
    header_parameters,
    axes_dict
):
    """
    Plot the output of get_bins.

    Args:
     -- bin_parameters (list)
          A list containing the locations and S/N of each bin
     -- parameters (dict)
          A dictionary of parameters read in from the motesparams.txt
          configuration file.
     -- spatial_lo_limit (int)
          Lower limit of the spatial region measured for S/N in
          get_bins()
     -- spatial_hi_limit (int)
          Upper limit of the spatial region measured for S/N in
          get_bins()
     -- data_2d (numpy.ndarray)
          2D spectroscopic data
     -- header_parameters (dict)
          A dictionary of parameters full from the datafile header.
     -- axes_dict (dict)
          A dictionary containing axes and axis metadata for the
          current extraction.

    Returns:
     -- None
    """

    # DIAGNOSTICS - Plot boundary locations of localisation bins on the
    # dispersion axis.
    spatial_axis = axes_dict["spatial_axis"]
    wavelength_start = axes_dict["wavelength_start"]
    spatial_floor = axes_dict["data_spatial_floor"]
    
    bin_line_lo_limit = spatial_lo_limit - ((spatial_hi_limit - spatial_lo_limit) * 0.2)
    bin_line_hi_limit = spatial_hi_limit + ((spatial_hi_limit - spatial_lo_limit) * 0.2)
    bin_line_location = (
        np.where(
            np.logical_and(spatial_axis > spatial_lo_limit, spatial_axis < spatial_hi_limit)
        )
    )
    
    draw_lines = []
    for b in bin_parameters:
        draw_lines.append(np.ones(len(spatial_axis[bin_line_location])) * b[0] + wavelength_start)
        draw_lines.append(spatial_axis[bin_line_location] + spatial_floor)
        draw_lines.append(np.ones(len(spatial_axis[bin_line_location])) * b[1] + wavelength_start)
        draw_lines.append(spatial_axis[bin_line_location] + spatial_floor)

    title = "2D Spectrum with Boundaries of Localisation Bins"
    show_img(data_2d, axes_dict, header_parameters, draw_lines, title)

    return None


def get_draw_lines(axes_dict, boundaries):
    '''
    Creates the extraction limit lines drawn in diagnostic plots
    presented when MOTES is run.
    
    Args:
     -- axes_dict (dict)
          A dictionary containing axes and axis metadata for the
          current extraction.
     -- boundaries (numpy.ndarray)
          Array containing the extraction limits to be plotted.
    
    Return:
     -- draw_lines (list)
          List of arrays containing the resampled extraction limits for
          plotting.
    '''
	
    dispersion_length = axes_dict["dispersion_axis_length"]
    wave_start = axes_dict["wavelength_start"]
    spatial_floor = axes_dict["data_spatial_floor"]
    
    i=0
    if len(boundaries)>2:
        extraction_limit_x_axis = boundaries[0] + wave_start
        i+=1        
    else:
        extraction_limit_x_axis = np.array(range(dispersion_length)) + wave_start
    
    draw_lines = [
        extraction_limit_x_axis, boundaries[0+i] + spatial_floor,
        extraction_limit_x_axis, boundaries[1+i] + spatial_floor,
    ]
     
    return draw_lines


def get_png_name(each_bin, wavelength_start, file_name, sky_ext):
    """
    Constructs the name string for a png figure of a  Moffat profile fit
    plotted against the data it has been fitted to.
    
    Args:
     -- each_bin (list)
          List containing the bounds of the current bin.
     -- wavelength_start (int)
          Wavelength column where the extracted spectrum starts.
     -- file_name (str)
          Name of the current input file.
     -- sky_ext (str)
          String indicating whether the Moffat fitting is for sky
          subtraction or spectrum extraction.
    
    Return:
     -- png_name (str)
          Name of the png plot file.
    """
    first_bin_col = str(each_bin[0] + wavelength_start)
    first_bin_col = ("0" * (5 - len(first_bin_col))) + first_bin_col
    last_bin_col = str(each_bin[1] + wavelength_start - 1)
    last_bin_col = ("0" * (5 - len(last_bin_col))) + last_bin_col
    png_name = (
        "m" + file_name[:-5] + "_" + sky_ext + "_" 
        + first_bin_col + "_" + last_bin_col + "_moffat.png"
    )
    
    return png_name


def plot_fitted_spatial_profile(
    spatial_axis, bin_data, 
    hi_resolution_spatial_axis, 
    bin_moffat_parameters, 
    image_start, 
    header_parameters,
    fig_file_name
):
    """
    Plot the spatial profile of a collapsed spectrum or a collapsed bin
    therein, and plot the Moffat function fitted to the data on top.

    Args:
     -- spatial_axis (numpy.ndarray)
          The spatial, or x, axis of the profile.
     -- bin_data (numpy.ndarray)
          The binned data that has been fitted with a Moffat profile.
     -- hi_resolution_spatial_axis (numpy.ndarray)
          The supersampled spatial axis used only for plotting purposes.
     -- bin_moffat_parameters (list)
          Parameters defining the Moffat profiel that was fitted to
          bin_data.
     -- image_start (int)
          The limit of the spatial axis after the original 2D spectrum
          was cut down to the region defined in reg.txt
     -- header_parameters (dict)
          A dictionary of parameters pulled from the header of the
          current datafile
     -- fig_file_name (str)
          Name of the current file being processed.

    Returns: None
    """
    
    parameter_box = dict(boxstyle="round", alpha=0.8, facecolor="white")
    fwhm = common.moffat_fwhm(bin_moffat_parameters[2], bin_moffat_parameters[3])
    textstr = "\n".join(
        (
            r"A = " + str(round(bin_moffat_parameters[0], 2)),
            r"c = " + str(round(bin_moffat_parameters[1], 2) + image_start),
            r"$\alpha$ = " + str(round(bin_moffat_parameters[2], 2)),
            r"$\beta$ = " + str(round(bin_moffat_parameters[3], 2)),
            r"B = " + str(round(bin_moffat_parameters[4], 2)),
            r"m = " + str(round(bin_moffat_parameters[5], 2)),
            "FWHM = " + str(round(fwhm, 2)) + " pixels"
        )
    )
    
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 4.5)
    ax.plot(
        hi_resolution_spatial_axis + image_start,
        common.moffat(bin_moffat_parameters, hi_resolution_spatial_axis),
        color="r",
        linewidth=3,
        label="Fitted Moffat Profile",
    )
    ax.text(
        0.66,
        0.95,
        textstr,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=parameter_box
    )
    ax.plot(spatial_axis + image_start, bin_data, color="k", label="Spatial Profile")
    ax.grid(linestyle="dashed", color="gray")
    ax.legend()
    ax.set_title("Spectrum Spatial Profile and Fitted Moffat Profile")
    ax.set_xlabel("Spatial Axis, Pixels")
    ax.set_ylabel("Median Flux, " + header_parameters["flux_unit"])
    plt.savefig(fig_file_name)
    plt.close()

    return None


def plot_sky_fit(
        good_sky_axis,
        good_sky_pixels,
        good_sky_errs,
        column_axis,
        column_sky_model,
        column_sky_model_err,
        col,
        file_name,
        flux_unit
):
    """
    Plot and save a figure showing sky pixels and the polynomial fitted
    to them.
    
    Args:
     -- good_sky_axis (numpy.ndarray)
          Spatial axis of the sky pixels in the current column.
     -- good_sky_pixels (numpy.ndarray)
          Sky pixels in the current column.
     -- good_sky_errs (numpy.ndarray)
          Uncertainties of the sky pixels in the current column.
     -- column_axis (numpy.ndarray)
           Spatial axis of the entire current column.
     -- column_sky_model (numpy.ndarray)
           Sky model fitted to the sky pixels.
     -- column_sky_model_err (numpy.ndarray)
          Uncertainty of the sky model fitted to the sky pixels.
     -- col (numpy.int64)
          Number of the current column.
     -- file_name (str)
          Name of the file currently being processed.
     -- flux_unit (str)
          Flux unit of the data.
     
    Return: None
    """
    colstring = ("0" * (5 - len(str(col)))) + str(col)
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 4.5)
    plt.errorbar(
        good_sky_axis, 
        good_sky_pixels, 
        yerr=good_sky_errs, 
        color="k", 
        linewidth=2., 
        label="Sky Data"
    )
    plt.errorbar(
        column_axis, 
        column_sky_model, 
        yerr=column_sky_model_err, 
        color="r", 
        linewidth=2., 
        label="Polynomial Sky Fit"
    )
    ax.grid(linestyle="dashed", color="gray")
    ax.legend()
    ax.set_title("Spectrum Spatial Background and Fitted Sky Polynomial: Column " + str(col))
    ax.set_xlabel("Spatial Axis, Pixels")
    ax.set_ylabel("Flux, " + flux_unit)
    plt.savefig("m" + file_name[:-5] + "_sky_fit_" + colstring + ".png")
    plt.close()

    return None


def plot_spectrum(wavelength_axis, extracted, header_parameters):
    """
    Plots the extracted spectra with matplotlib.
	
    Args:
	 -- wavelength_axis (numpy.ndarray)
	      Numpy array containing the wavelength axis of the spectra.
	 -- extracted (list)
	      List of numpy arrays containing the extracted spectrum data
	      and the associated uncertainties.
	 -- header_parameters (dict)
	      A dictionary of parameters pulled from the header of the
          current datafile
    
    Return: None
    """
    
    # DIAGNOSTICS, EXTRACTED SPECTRUM
    a_label = "Aperture Spectrum"
    o_label = "Optimal Spectrum"
    plt.figure(figsize=(9, 6))
    plt.errorbar(wavelength_axis, extracted[2], extracted[3], color="k", marker=".", label=a_label)
    plt.errorbar(wavelength_axis, extracted[0], extracted[1], color="r", marker=".", label=o_label)
    plt.grid(alpha=0.5, linestyle="dotted")
    plt.title("Extracted 1D Spectrum")
    plt.ylabel("Flux, " + header_parameters["flux_unit"])
    plt.xlabel("Wavelength, " + header_parameters["wavelength_unit"])
    plt.legend()
    plt.show()


def show_img(data_2d, axes_dict, header_parameters, draw_lines, title):
    """
    Takes an input image and line data to be drawn on that image and
    creates a figure to be shown on screen.

    Args:
     -- data_2d (numpy.ndarray)
          The input image.
     -- axes_dict (dict)
          A dictionary containing the spatial and spectral axes of the
          input image.
     -- header_parameters (dict)
          A dictionary containing the header parameters of the input
          image.
     -- draw_lines (list)
          A list of line data to be drawn on the input image.
     -- title (str)
          The title of the figure.

    Returns: None
    """
    
    wavelength_axis = axes_dict["wavelength_axis"]
    wavelength_start = axes_dict["wavelength_start"]
    wavelength_end = wavelength_start + len(wavelength_axis)
    spatial_axis = axes_dict["spatial_axis"]
    spatial_floor = axes_dict["data_spatial_floor"]
    spatial_start = spatial_axis[0] + spatial_floor
    spatial_end = spatial_axis[-1] + spatial_floor

    power = int(np.floor(np.log10(np.abs(np.nanmean(data_2d)))))
    data_2d = copy.deepcopy(data_2d) / 10**power

    figwidth = 10.0
    fig = plt.figure(figsize=(figwidth, figwidth / 1.9))
    gs = gridspec.GridSpec(18, 33)
    ax = plt.subplot(gs[:, :32])
    color_axis = plt.subplot(gs[1:, 32])
    masked_data2d = np.ma.masked_where(data_2d == 0, data_2d)
    cmap = matplotlib.cm.inferno
    cmap.set_bad(color="red")
    
    v_max = np.ma.median(masked_data2d) + (0.5 * np.ma.std(masked_data2d))
    s = ax.imshow(
        masked_data2d,
        aspect="auto",
        vmin=0,
        vmax=v_max,
        origin="lower",
        cmap=cmap,
        extent=[wavelength_start, wavelength_end, spatial_start, spatial_end]
    )

    for i in range(int(len(draw_lines) / 2)):
        ax.plot(draw_lines[i * 2], draw_lines[(i * 2) + 1], color="white")
    cbar = fig.colorbar(s, cax=color_axis)
    cbar.ax.yaxis.set_offset_position("left")
    cbar.ax.set_ylabel("Pixel Flux, x10^" + str(power) + " " + header_parameters["flux_unit"])
    ax2 = ax.twiny()
    ax2.plot(wavelength_axis, data_2d[0, :], alpha=0)
    ax2.set_xlim(wavelength_axis[0], wavelength_axis[-1])
    ax2.set_xlabel("Wavelength, " + header_parameters["wavelength_unit"])
    ax.set_ylim(spatial_axis[0] + spatial_floor, spatial_axis[-1] + spatial_floor)
    ax.set_ylabel("Spatial Axis, Pixels")
    ax.set_xlim(wavelength_start, wavelength_end)
    ax.set_xlabel("Dispersion Axis, Pixels")
    plt.title(title, y=1.095)

    # Add interactive scaling bar to figures if the 2D spectrum isn't
    # flux calibrated. Again for some reason teeny tiny numbers cause
    # things to break.
    fig.subplots_adjust(bottom=0.2)
    ax_vmin = plt.axes([0.1, 0.05, 0.8, 0.03])
    ax_vmax = plt.axes([0.1, 0.01, 0.8, 0.03])
    slider_min = Slider(
        ax_vmin,
        "LowCut",
        0,
        np.ma.max(masked_data2d) - 1,
        valinit=0.0,
        valstep=0.001 * np.ma.max(masked_data2d),
    )
    slider_max = Slider(
        ax_vmax,
        "HighCut",
        1.0,
        np.ma.max(masked_data2d),
        valinit=np.ma.median(masked_data2d) + (3.0 * np.ma.std(masked_data2d)),
        valstep=0.001 * np.ma.max(masked_data2d),
    )

    def update(val):
        vmax = slider_max.val
        slider_min.valmax = vmax - 1
        vmin = slider_min.val
        slider_max.valmin = vmin + 1

        s.set_clim(vmin, vmax)
        cbar.ax.set_ylabel(
            "Pixel Flux, x10^" + str(power) + " " + header_parameters["flux_unit"]
        )
        fig.canvas.draw_idle()

    slider_min.on_changed(update)
    slider_max.on_changed(update)

    plt.show()

    return None

