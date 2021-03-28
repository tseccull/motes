# Change Log
All notable changes to this package will be documented here. 
This changelog follows the format described [here](https://keepachangelog.com/en/0.3.0/). [Semantic Versioning](https://semver.org/) is followed.

## 0.3.0_dev - 2021-03-28
Updates by D. Kiersz

... which probably need a second readthrough.

#### Progress at the Abyss

This major _development_ version is based on MOTES version 0.1.1 and it incorporates *some* changes introduced in 0.2.0. The development branch arose purely from X-Shooter spectrograms and more testing may be necessary with other instruments to ensure stability at release. There are also hope that this entry will summarise the GME 2.0 (aka GME Plus) report from January 2021, unifying our development information here. As mentioned, some changes seen in MOTES 0.2.0 are imported here (for example, the use of **np.nan** functions). Some other changes conflicted with this version that may need looking at. MOTES' 0.2.0 changes to cosmic ray detection and masking were **not** imported here. Additionally, this version may need some stylistic changes (i.e. cleanup) before release.

At last, this version incorporates an adapted version of an optimal extraction algorithm from Horne (1986) and changes to fitting routines. This version also introduces a handful of bug-fixes from the early days of the GME (for example, common.get_bins() corrections). A number of changes were fixed by the Tom in MOTES 0.1.0 while our development versions were done in solitude.

### Added
- A description of reg.txt was added to **motesparams.txt**
- Information on optimal X-Shooter dispersion (i.e. wavelength) limits and two new parameters (*-SUBTRACT_SKY_SINGLE_BIN* and *-ABBA*) were added to **motesparams.txt**. A full description will follow in the manual.
- Following functions have been added to common.py. Please refer to the code (and the upcoming manual!) for more information. These extra functions come from GME Plus' **xshoo_extra.py** and are now part of MOTES' **common.py**.
    * **filter_data()** is a simple method that takes in the data2D and errs2D and outputs frames where any NaN or Inf values turn to 0.0. This is useful for fitting data.
    * **qual_to_molecfit_compat()** is a method that ensures that MOTES doesn’t ignore the pixels that have a QUAL extension which are deemed to be interpolated or have a bad B-Spline fit as a result of ESO pipelines. This method is also used to convert spectrograms QUAL frame to a GME/MOTES complaint binary (i.e. 1 is good pixel, 0 is bad pixel).
    * **optimal_extraction()** adds in modified Horne's (1986) algorithm for extracting spectra optimally on bin-by-bin basis. A full description will be included in the manual but most of it should be in code.
    * **GME_Moffat1D_curvefit()** and **GME_Moffat1D_profile()** are two different profiles used by moffat_weighted_curvefit() and the first Moffat bin. The choice depends on the nature of the spectrogram.
    * **moffat_weighted_curvefit()** is a new non-linear weighted fitting using SciPy's **curve_fit()** and will be used for fitting profiles to the collapsed bins/spectra.
- Comments in-code for development purposes - most of which will need transferring onto a manual before release.
- An overhaul of the .yml Anaconda package dependancies to include all required modules. It was tested on a separate Linux machine with Anaconda so the package *should* work out of the box. The next step would be to add a bash script that will install the enviroment.

### Changed
- This changelog and **README.md**
- Minor changes to the code layout in accordance to AUTOPEP8 standards. Usually, this is done automatically in Microsoft VSCode (or Code-OSS, as used in many Linux distros).
- In the **harvester.harvest_xshoo()** add **common.qual_to_molecfit_compat(_imgqual_)** for converting framedict['qual'] frames to MOTES compliant binary and ensuring certain types of pixels are not treated as bad pixels.
- Force the TkAgg backend to enable saving any displayed figure within the package as an image. It also ensures that the diagnostic plots are working all the time. This is located in **motes.py**.
- During the 'First Bin' phase of MOTES, the spectrogram is collapsed via a median *and* the uncertanities of the results spatial profile are determined using the median absolute deviation (MAD). Across a whole 2D spectrum, this provides a robust measure of the average profile more resistant to outliers at the ends of each arm.
- In the NODDING reduced spectrograms of the X-Shooter, the data frame used to determine the datascale variable can include negatives. Correct for that by using the **np.abs()** and np.percentile function with the 86th percentile.
    * The use of np.abs() is similar to MOTES version 0.2.0, submitted during the development of this version.
    * (Conflict with 0.2.0) In version 0.2.0 **np.percentile()** is replaced with **np.nanmedian()**. However, specifically in NODDING spectrograms there is a possibility of a median being closer zero. To ensure that datascale is significant, using 86% of the _dataiscollapse_ yields more consistent results.
- During the first bin profile fitting (i.e. whole 2D spectrogram collapsed on the dispersion axis), use a *weighted* least-squares fit a Moffat distribution instead, as done in **common.moffat_weighted_curvefit()**. 
- (Conflict with 0.1.0) This relates to the **first use of common.extraction_limits()** right after obtaining the first profile, in order to establish what area of the spectrogram is used to calculate the threshold S/N. This which in turn is used to determine the bins by **common.get_bins()**. 
    * In old GME versions, these variables were subsequently reused in sky subtraction routines. These variables were overwritten later in code, which was a bug). 
    * In GME Plus, that was fixed by renaming variables. Additionally, the _lowext_ and _highext_ were determined with a constant 2.0 FWHM from the centre (instead of the original GME's 1.5 FWHM). In return, this would make the bins larger then expected. However, this change complied with the optimal extraction routines to utilise most of the objects flux for photometric accuracy and should be remembed when choosing this parameter.  
    * In MOTES, the constant FWHM have now been replaced with _params['-CR_FWHM_MULTIPLIER']_ for the user to alter. That means MOTES' 0.1.0 _lowext_ and _highext_ seemed to be unique to cosmic ray masking and removal since sky subtraction is now a module on its own.
    * To resolve any potential overwritting, _lowext_ and _highext_ are saved as *lowext_median* na *highext_median*. I also rename other instances of **common.extraction_limits()** products as appriopiate to avoid further issues issues.
- In **common.mask_cosmic_ways()** add the following two bugfixes:
    * If the limits _lext_ and _hext_ are outside the spectrogram (eg. individual ABBA frames reduced in STARE mode of the X-Shooter ESO pipeline), then add a correction that to force _lext_ and _hext_ to be within the frame.
    * As a result of the above, if there are only two sections of the 2D spectrum (an individual ABBA frame, where the object spectrum is not centered as a result of dithering), then correct the assembly of CR masks before returning any masks.
- **common.get_bins()** function consists of two loops which determine the bin span and spatial limits in dispersion and spatial direction respectively. Assume the following modifications are applicable for both:
    * Using **common.filter_data()**, filter out any NaN/Inf data points before determining bins and ensure the errors are positive.
    * Rename _datacol_ and _errcol_ to _data_slice_ and _errs_slice_ respectively.
    * After a single bin collapse, there may be scenarios where spatial limits defined in the ’First Bin’ _slow_ < 0 and _shigh_ > _bindatacol.shape[0]_ (again, possible in individual A and B STARE frames from an ABBA dithering pattern). Add an exception for that.
    * Add an exception added to handle anomalies in any results prior to estimating the S/N.
        * if ((_rssnoise_ == 0) or (np.isnan(_rssnoise_)) or (np.isnan(_signal_))): _snrestimate_ = 1
    * Add an exception to handle diving by zero or an invalid value when determining the CR/bad pixel replacement near the end of the loop. Otherwise, the scale variable may be a NaN which breaksthe  replacement of bad pixels.
    * (Conflict with 0.2.0) In response to above issues, MOTES version 0.2.0 uses np.nansum() and np.nanmax(). At first glance they should serve the same purpose to above changes. This should be verified later and we will change it as necessary.
- During the sky subtraction phase, change the following:
    * Rename _extractronlimits_ to *extractionlimits_sky*
    * Modify the fitting procedures similarly to the extraction phase explained below this section.
    * Rename _LowExt_, _HighExt_, _fwhm_ and _centre_ to *lowext_sky*, *highext_sky*, *fwhm_sky* and *centre_sky* respectively.
    * Add a commented out p-value diagnostic for judging the overall quality of fittings with a distribution of all bins. Needs more work...
    * _skyextractionlims_ is the result of **common.interpolate_extractionlims()**. This is then converted to standard resolution.
    * I added an option to determine and use a single extraction aperture for sky subtraction, similar to the classic extraction means. In rare circumstances, the sky lines were bright to the point of not allowing sky/object limits to be determined correctly. This is done via assigning *-SUBTRACT_SKY_SINGLE_BIN=1* in **motesparams.txt**
- During the extraction phase, modify and add the following:
    * Rename _extractionlims_ to _extractionlims_optimal_ for clarity.
    * Add a *p_value_extraction* array for statistical measures.
    * In the bin fitting loop, add **common.filter_data()** for reasons explained previously.
    * For fitting, use **np.nanmedian()** to collapse the bin in the dispersion direction and use MAD for the measure of spread (of each row in the dispersion direction). Additionally, ensure that the modification for GMOS remains.
    * Rename _LowExt_, _HighExt_, _fwhm_ and _centre_ to _lowext_optimal_, _highext_optimal_, _fwhm_optimal_ and _centre_optimal_ respectively.
    * Rename _binpars_ to _bin_save_ which serves for file saving.
    * Add a commented out p-value diagnostic for judging the overall quality of fittings with a distribution across all bins. Again, needs more work.
    * Rename results of **common.interpolate_extraction_lims()** from _finalextractionlims_ to *finalextractionlims_optimal*
    * Maintain supersampling of limits (and within their diagnostic plots) however do bare in mind these will be rounded to the upper/lower pixel in optimal extractions. This is reflected in the drawlines variable where np.ceil() and np.floor() are used and subsequently plotted.
    * Add **common.optimal_extraction()** with relevant inputs. This will need a through description in the manual. Some details are availiable in code.
- Upgrade the final plot.
- Save results of optimal extraction 'optdata1D' and 'opterrs1D'. We can add stddata1D and stderrs1D in the future.
- If sky subtraction is off make sure that skybinpars = 0.0 and skyextlims = 0.0 before attempting to save it.
- Minor edit to common.make_wav_axis() for clarity.
- In common.mask_cosmic_rays() and common.subtract_sky() use np.nan functions when calculating median etc. 

### Removed

- In **harvester.harvest_xshoo()**, remove the old qual2D conversion.
- I commented out old fitting procedures.
- Previous standard extraction code and its supersampling. 
- Removed _crmask_multiplier_=1. from **common.subtract_sky()** - nothing in the method used it. 

### Known Bugs

i.e. things that need ironing out in the future

- If *['-MASK_CR']* = 0 and *['-DIAG_PLOT_CRREPLACE']* = 1, code crashes as there is no error handling to contradicting parameters.
- Variable _chipgap_ seemgly required for GMOS is not present in this iteration of the optimal extraction routine.
- When attempting to supersample data before extractions, the numerical error between standard and optimal extracted flux becomes significant. On further inspection, this was caused by poor profile fitting given supersampling. For now, no supersampled data is extracted and only whole pixels are taken into account.
- I am also not sure if the description of **common.subtract_sky()** is valid anymore.

### From GME Plus

This is the record on what was changed in my version of the GME (GME Plus) independently from the main branch of MOTES 0.1.1 (and below). Any of these changes have some implementation in MOTES already thus left untouched.

- In the **harvester.py**, the user does already get information about both the spatial and the _dispersion_ pixel resolution which aligns with GME Plus.
- The early version of GME (prior to March 2020) used a **np.sum()** instead of a **np.nanmedian()** to collapse the whole spectrum and determine the overall spatial profile regardless of wherever the method used an weighted fitting or not. Both GME Plus and MOTES now use a median method instead. 
- Minor conflict on the _datascale_ definition between our versions should be now resolved as explained above.
- The derivation of _lowext_ and _highext_ when using **common.extraction_limits()** after the first bin differed in both MOTES and GME Plus. Renaming subsequent variables and sticking to MOTES' FWHM variable *params['-CR_FWHM_MULTIPLIER']* sufficed.
- In **common.get_bins()**, instances of **np.sum()** in GME where changed to **np.nansum()** in both GME Plus and MOTES. This is good and left alone.
- In the primodal GME (specifically **xshoo()**), a simple sum was derived before bin fitting. MOTES and GME Plus have changed that to a nanmedian (the latter had both options implemented with errors). For simplicity, the **np.nanmedian** is used along with the MAD.
- In both GME Plus and MOTES, _bgfwhm_, _bglowext_ and _behighext_ were removed since they served no purpose.
- Both GME Plus and MOTES solved the issues of crashing when loading onto Fits Viewer (FV) if sky subtraction did not take place and when the file was saved. This was a separate issue to not not assigning _skybinpars_ = 0.0 and _skyextlims_ = 0.0 before saving.
- In **common.mask_cosmic_rays()** and **common.subtract_sky()** both GME Plus and MOTES used nan functions.

## 0.1.1 - 2021-03-04
Updates by T. Seccull

### Changed
- Updated wording of python compatibility in README.md.

### Removed
- Unused lines related to marking the chip gaps in GMOS spectra with NaNs.


## 0.1.0 - 2021-03-03
Updates by T. Seccull

#### The Dark Ages
The majority of this package was written by Tom Seccull during his PhD at 
Queen's University Belfast from 2015-2019, under the tutelage of Wes Fraser (at 
HIA Victoria, BC, at time of writing). The process was messy, haphazard, and 
initially documented only within the fuzzy mind of the lead developer (needless 
to say better practices have been adopted since then). For much of 
its development MOTES was known as GME (the Grand Moff Extractor); the current 
name has been selected to better represent the function of the package and 
increase its accessibility to new users. Since March 2020 the primary source of 
documentation for this software has been Chapter 3 of Seccull's PhD Thesis 
titled ["Revealing Refractory Materials on Trans-Neptunian Objects and Centaurs 
via Reflectance Spectroscopy"](https://pure.qub.ac.uk/en/studentTheses/revealing-refractory-materials-on-trans-neptunian-objects-and-cen). The thesis text will be 
released from embargo in the summer of 2022, but better documentation for MOTES 
is expected to be provided before then. All changes listed here relay 
developments to MOTES following 2021-03-01.

### Added
- This changelog.
- if name == main has been introduced to motes.py so the MOTES functions may be 
called within another python script, or within a python session.
- Demo data for FORS2, GMOS, and X-Shooter, all bright targets.
- Scale bars on the figures now allow users to adjust the cut in the plotted 
images.

### Changed
- The entire package has been made modular and processes all spectra from 
different instruments in the same way. The process generally follows that define 
for X-Shooter spectra in GME.
- harvester.py has been created (modified from the version that handled FORS2 
spectra) and is capable of unpacking data and relevant FITS header parameters 
from any FITS formatted spectrogram. data_harvest() calls a more specialised 
harvester function depending on which instrument a spectrogram has come from. 
harvester.py has been designed to allow the development and inclusion of new 
specialised harvester functions for other instruments that can extract and 
convert the relevant data from other instruments into a format that MOTES can 
process.
- Functions readparfile() and readregions(), which respectively read 
motesparams.txt and reg.txt have been moved to their own module called 
startup.py
- Controls for how the data frames are initially sliced have been removed from 
motesparams.txt entirely; instead the user must define the number of rows to 
remove from the spatial axis of the spectrum, and the wavelength range of the 
spectrum for each input file in reg.txt
- motesparams.txt must now be present in the same directory as the input data 
files. This has been done to minimise the need to change the parameters for 
different setups between testing runs. It should also make life easier for users
who may be using MOTES with multiple datasets requiring different parameters.
- The main script that was forsextraction() or xshooextraction() has been 
incorporated into motes.py. The motes() function is now the main function of the 
package. Sky subtraction processes and cosmic ray handling have respectively 
been moved to their own functions within motes.py; this is to reflect the fact 
that they are optional and are not required for motes to run.
- Small rearrangements in the order of the output file header metadata have been
done.
- Data frames are now almost always stored in a dedicated dictionary referred to
as framedict or some similar name depending on which function it's called from.
- Similarly to data frames, axes and parameters related to them are stored in 
axesdict or similar.
- All references to GME changed to MOTES; this includes references within the 
Python scripts themselves and within comments.
- Version number has been dropped from v0.9.0 to v0.1.0, to better reflect 
reality (I'm also dumb and should have read up on semantic versioning before 
sticking a number on this, ha!).
- Output frames like skymod, skybins, skyextractionlims and crmask are only 
added to the frame dictionary if those sections of code are run. This may
fix a bug where spectrograms that are already skysubtracted are extracted
but error out during saving.
- Output files from MOTES are prepended with 'm' rather than having '_GME_'
included in the middle of the filename.
- MOTES no longer walks through an ESO-like directory tree to find input files.
All input FITS files are expected to be present in the current working 
directory, and ordered in the same way as their associated regions in reg.txt.
