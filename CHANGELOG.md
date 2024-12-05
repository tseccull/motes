# Change Log
All notable changes to this package will be documented here. 
This changelog follows the format described 
[here](https://keepachangelog.com/en/1.0.0/). [Semantic Versioning](https://semver.org/) is 
followed.


## 1.0.0 2024-12-05
Updates by T. Seccull

This update makes MOTES v1.0.0 Its core is done, or rather I am done with it. This
update is also the reason why I will never be hired to a professional software engineering position. 
The source code has been almost entirely gutted and replaced. I am also removing CI/CD 
(sorry Dominik) as frankly this package will never meet the bar for linting and formatting that 
is required in industry (not without a prohibitively time consuming redesign). It works, and is now
good enough, even if it isn't the prettiest or easiest codebase to manage. From here, the only 
major updates should be ones adding/updating modules for I/O for instruments. Docs also need to be 
added. Various other changes and fixes have been implemented. Most of which I have pulled out of 
the swamps of my recent memory and outlined below. No doubt, I have forgotten some of the more minor 
fixes or changes I've made over the past several months. Once documentation and minor testing are 
completed it will be possible to tentatively release this as v1.0. I live in hope.

### Added
- Python's logging module is now used to log all messages both to files and to the terminal.
  By default, warnings and critical messages are printed in the terminal, and info level
  messages are logged to a file. Verbose mode can be called to print all messages above info
  level to the terminal as well as to the log file. A log file is produced for each new
  extraction in the chain when MOTES is run. A central `motes.log` file also contains
  general information about the run.
- Diagnostic plotting for the sky subtraction has now been implemented. When activated
  a plot of each sky fit in each wavelength column and its uncertainties are plotted to
  a figure that is saved in the working directory.

### Changed
- The source code has overall been cleaned, broken into more functions, and partitioned into
  modules roughly by function. Any changes that need to be made in future should be much more
  straightforward. A few instances of repeated code have also been removed, shrinking the codebase
  a little bit. The only instrument module currently fully compatible with this format is `gmosio.py`.
  Modules for other instruments still need to be updated.
- MOTES can now be called from the command line. `motesparams.txt` is no longer used as an input
  file. All parameters in `motesparams.txt` are now either called from the command line using
  argparse, or are supplied in a Namespace object to the `motes()` function in `motes.py`. Docs
  will be written to outline these options to users.
- Changes to `gmosio.py` now mean that MOTES only accepts GMOS data if it has been reduced with
  DRAGONS. Gemini IRAF is no longer supported by MOTES. As DRAGONS has been the default data
  reduction method for GMOS longslit data for a few years now, this should not be a massive issue.
- `reg.txt` has been relocated to `\inputs` alongside the input files. The order in which input
  files are processed is no longer done in the alphabetised order of files in the inputs directory.
  Each line of `reg.txt` contains five comma-separated variables: (1) the number of spatial rows to
  slice of the top of the 2D data frame, (2) the number of spatial rows to slice off the bottom
  of the 2D data frame, (3) the lower limit of the desired wavelength range, (4) the upper limit
  of the desired wavelength range (5) the name of the file containing the data to be processed.
  With this new layout, files can be processed more than once within a single run of MOTES,
  providing a rudimentary MOS functionality.
- `motes.py` has had a number of changes performed, primarily to remove repeating sections of code.
- `motesio.py` now contains all I/O functions that are not directly specific to an instrument. This
  includes `read_regions()`.
- `common.py` is greatly reduced and only contains functions that are called by functions in more
  than one other module.
- `diagnostics.py` contains all functions related to plotting and printing of diagnostic data that
  isn't logging.
- `extraction.py` contains all functions that are directly involved with spectrum extraction.
- `logs.py` contains a function defining the configuration of the `motes` logger.
- `sky.py` contains all functions directly involved with sky subtraction.
- `tracing.py` contains all functions related to localising the spectrum and sky regions of the 
  2D data.
- Sky subtraction has been overhauled to use `numpy.polyfit()` to perform the fitting of polynomial
  functions to the background that can be subtracted to remove the sky contribution. Polynomials of
  any order can be used, and setting the order to 0 leads to a simple median subtraction being performed.
  For each spatial column polynomial fits to the background are bootstrapped within the uncertainties of
  the data to produce an estimate of the median background and uncertainties on the sky subtraction.
  These uncertainties are propagated into the sky subtracted error frame and is accounted for further
  down the line during the spectrum extraction.
- Like the sky subtraction diagnostics plots, diagnostic plots of the Moffat profile fits to each
  localisation bin are saved as individual figures in the working directory of MOTES. It is no longer
  possible to view the moffat fit while MOTES is running.

### Removed
- `.githubworkflows` has been removed from the repo, along with its contents. New workflows for basic
  linting may be added again at a later date, but right now they are a hindrance to completion of this project.
- The github site directory has been removed. Standard readthedocs format will be used for documentation.

## 0.4.7 2024-04-XX
Updates by T. Seccull

This is a large update that provides various bug fixes, removes unused code and variables, and adjusts 
almost all variables to conform to snake_case. Many changes are based on the analysis of linting 
software like flake8 and pylint. Linting isn't yet fully complete, but this update hopefully represents a 
solid start to that work. Large portions of the codebase can still be simplified and organized with a 
cleaner structure. MOTES is still fully functional, and has been tested following the changes in this 
update. I did make a few changes to the motes.yml file, but reverted them as Dominik is working on a 
better version for a future update; its best I leave the CI stuff to someone who knows what they're doing.

### Fixed
- Bad pixel columns containing no good data (this doesn't include flagged columns in detector chip gaps)
  were causing divide by zero runtime errors in `optimal_extraction()` in cases where the `MEDIAN` sky
  subtraction mode was used. The extraction would still complete successfully with the bad data removed
  from the spectrum, but the script's behaviour was not consistent between all sky subtraction modes.
  To fix this I added an if statement to catch columns of all-zero data and extract the optimal spectrum
  with zero counts and zero errors before continuing on to the next good column. This if statement
  bypasses the actual optimal extraction that would throw the runtime warning if the bad data was fed to
  it.
- Using the `LINEAR`, `POLY2`, or `POLY3` sky subtraction modes would cause the sky model frame to be
  saved with a transposed (i.e. x and y axis both flipped) orientation. I added an extra `.T` to the
  `sky_model` frame in the frame dictionary before the end of `subtract_sky()` for these modes to fix
  this. 

### Changed
- Some instances of boolean logic statements (like some cases of `np.where()` for example) have been
  reformatted as recommended by flake8.
- Minor formatting changes have been done to satisfy flake8 and pylint.
- 200+ keywords, variable names, and funtion names across the codebase have been carefully adjusted to be
  more descriptive and follow the snake_case convention.
- The file name format for the output files for MOTES have been changed. From now on the file name will be
  the same as the name of the input file, but prepended with an `m`.
- `extraction_limits()` function name changed to `set_extraction_limits()` to avoid conflicts with internal
  variable names. The maintain alphabetization of the funtions in common.py, `set_extraction_limits()` has
  been moved further down in the file.
- 

### Removed
- Trailing spaces at the ends of some lines.
- Previously deprecated and commented out sections of the get_bins() function in common.py have finally been
  removed. This includes all sections of code triggered by `-REPLACE_CRBP` being set to True in
  motesparams.txt. As `-REPLACE_CRBP` is no longer used, it can be removed from motesparams.txt and should
  not be included in that file alongside any demo data that gets added to the repo.
- A few stdout print statements were being duplicated when the code was run because similar statements were
  being printed both before the call to a function and at the start of the function itself. In each case
  one of these duplicates has been removed.
- References to `-INTERP_KIND` have been removed, as the extraction limits are now hard coded to be
  be interpolated linearly between the centres of each bin fitted with a Moffat profile. The `-INTERP_KIND`
  keyword appeared only in the save_fits() function, and has no bearing on the function of MOTES; it should be
  removed from `motesparams.txt`.
- A few unused variables have been removed.

## 0.4.6 2024-02-13
Updates by T. Seccull & D. Kiersz

Removed disabled github actions and unneccessary CI/CD components from the repo. A robust linting 
pipeline has been added.

### Added
- `.github/workflows/motes.yml` is a linting pipeline that activates on push and pull requests. It runs flake8,
  pylint, black, and bandit.

### Removed
- `.github/workflows/motes-cicd.yml` has been replaced with `.github/workflows/motes.yml`
- `.github/workflows/labeler.yml` has been removed.
- `.github/workflows/delete-old-runs.yml` has been removed.

## 0.4.5 2024-02-09
Updates by T. Seccull

Running tests on MOTES in its current main branch form to make sure basic functionality is good. 

### Changed
- Tidied up formatting for script comments, function descriptions, the CHANGELOG, and text printing 
  to the screen for consistency and ease of reading in a text editor.
- Converted some `print()` statements in `common.py` to `sys.stdout.write()` for consistency with 
  the rest of MOTES.
 

## 0.4.4 2023-02-16
Updates by D.Kiersz

We are not ready for full CI/CD yet, but this small update is a step in the right direction as it 
introduces a CI pipeline. Deployment to [PyPI](https://pypi.org/) is the end goal, we are unable to 
register new projects at this time.

### Added
- A CI pipeline has been added to the repository. This will run the linting/vurnebility checks 
  against the code on every push to the main branch, but not enforce them for now. This should help 
  improve the quality of the code down the line.
- Added a workflow to remove old workflows.
- Badges to README.md.

### Changed

- Improvement to .gitignore by using a Python template from 
  [gitignore.io](https://www.toptal.com/developers/gitignore).
- Formatted with `black` and `isort` to conform to PEP8.
- Job names in pipelines are more systematic.

## 0.4.3 2023-04-21
Updates by T. Seccull

### Fixed
* A bug was fixed in the `optimal_extraction()` function in `common.py` related to the enforcement 
  of positivity required for the spatial profile when performing optimal extraction as described by 
  Horne (1986). Previously this was misinterpreted as setting the both the spatial profile and its 
  associated data to zero where data pixels were negative. Instead, the enforcement of positivity 
  applies only the the spatial profile. Because, however, the spatial profile in MOTES is a model 
  Moffat PSF, it has no negative values to begin with and positivity therfore does not need to be 
  enforced. The effect of this change will be minimal for spectra of bright targets that have few 
  negative pixels in their spatial profile. This update prevents overestimation of counts in faint 
  regions of the spectrum with higher numbers of negative pixels. 

## 0.4.2 2023-03-05
Updates by D.Kiersz

### Changed
* `inputs` folder to store fits files, change output file format name to include 1D.
* Conversion to docstrings adhering to Google's style guide: 
  https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings, added 
  to all functions. This is the first iteration of in-line documention that will need refining.
* Refactored loops in startup.py.
* `Black` formatter applied across all python files.
* Remove redundant param_dict from data_harvest()
* Remove redundant imports from common.py
* Remove redundant `interpkind` var from interpolate_extraction_lims and related functions
* Remove redudnant variable declaration `wav_min = scihead["CRVAL1"]` from harvester.py
* Import order of modules to conform to PEP8.
* Instances of exit() to sys.exit() to conform to PEP8.
* Updated link in CHANGELOG.md

## 0.4.1 2023-02-15
Updated by D. Kiersz

### Added
* Installation environment for Anaconda and dependencies with instructions. This is the first step 
  for contenarising the application.

### Changed
* Updated README.md, including a spelling mistake.
* In Removed the following line 
  `sys.stdout.write('     Extraction of ' + direc + ' could not be completed.\n')`.
* In `startup.py` use a list-comprehension instead of a loop.

## 0.4.0 2022-07-27
Updates by T. Seccull

### Added
- An optimal extraction routine developed by Dominik Kiersz has now been fully integrated into 
  MOTES. Some changes have been made to the original design by TS to increase the efficiency of the 
  process, and ensure that is compatible with data from multiple instruments. It still operates on 
  the original fundamental principles developed by DAK, however. MOTES now produces both aperture, 
  and optimally extracted spectra by default.

### Changed
- The comments related to many functions have been updated to more clearly describe their inputs, 
  outputs, and purpose. This work is ongoing along with removing deprecated sections of code.
- The GMOS harvester function has been updated to handle the detector chip gaps more efficiently.
- There have been multiple changes to how the extraction limits are defined to account for the new 
  optimal extraction routine. The main one being that the data is no longer supersampled prior to 
  being extracted in the aperture method.  

## 0.3.0 2022-04-06
Updates by T. Seccull

### Added
- Harvester function has been added for the FLOYDS spectrographs on the Las Cumbres Observatory 2 m 
  telescopes.

### Changed
- Astroscrappy has been run on the GMOS demo data to remove cosmic rays.
- Demo data has been removed from the github repository for now as the GMOS frames are too large. 
  Return of the FORS2 and X-Shooter demo data frames to the repository should be fine in the near 
  future.

### Deprecated
- All cosmic ray handling functions and lines have been commented out pending removal. This is 
  because there are much better software tools, (e.g. Astroscrappy) that can be used for removing 
  and replacing cosmic rays. Testing of the software on the demo data shows that it is fully 
  functional without the cosmic ray handling sections. These regions will need to be tidied up.


## 0.2.1 2021-03-26

### Changed
- Colormap for plots updated to inferno. Using a perceptually uniform sequential colormap really 
  makes a difference in clearly defining structure in the spectrum.

## 0.2.0 - 2021-03-19
Updates by T. Seccull

### Changed
- Simplified the way the GMOS harvester handles incoming qual frame and converts it so good 
  pixels=1 and bad pixels=0, while ensuring the chip gaps aren't flagged as bad pixels.
- Overhauled cosmic ray handling so CRs and badpixels can be replaced by the input quality frame 
  alone without the need for motes to use its, frankly piss-poor, CR detection routine. If software 
  like astroscrappy is used to create a cosmic ray mask, it needs to be combined with the 
  spectrogram's input quality mask for motes to recognise the detected CRs as bad.
- In `motesparams.txt`, the `-MASK_CR` keyword is now `-IDENT_CR` and is specifically used to refer 
  to the motes CR identification routine independently of the `-REPLACE_CR` keyword. Setting either 
  `-IDENT_CR` or `-REPLACE_CR` to True (1) results in an updated `cr_handlin()` function in 
  `motes.py` being called.

### Fixed
- Fixed incorrect header parameter dictionary call in `save_fits()` when recording the wavelength 
  unit to the header of the output file.
- When determining datascale in `motes()`, the base 10 log of the absolute value of the median 
  collapsed profile is now calculated. Previously this would fail if the median value turned out to 
  be negative.
- Fixed outdated method of calling `common.get_bins_output()` in the case where CR masking is not 
  done. This has now been updated in line with other calls to this function from `motes()`.
- Fixed a bug where cosmic ray replacement would completely fail on spatial pixel columns with no 
  valid data (e.g. all nans or zeros), and would wrongly propagate bad pixel flags to adjacent 
  columns with good data. `common.get_bins()` now simply ignores columns with all bad data, as they 
  are pretty much unrecoverable. 

## 0.1.1 - 2021-03-04
Updates by T. Seccull

### Changed
- Updated wording of python compatibility in `README.md`.

### Removed
- Unused lines related to marking the chip gaps in GMOS spectra with NaNs.


## 0.1.0 - 2021-03-03
Updates by T. Seccull

#### The Dark Ages
The majority of this package was written by Tom Seccull during his PhD at Queen's University 
Belfast from 2015-2019, under the tutelage of Wes Fraser (at HIA Victoria, BC, at time of writing). 
The process was messy, haphazard, and initially documented only within the fuzzy mind of the lead 
developer (needless to say better practices have been adopted since then). For much of its 
development MOTES was known as GME (the Grand Moff Extractor); the current name has been selected 
to better represent the function of the package and increase its accessibility to new users. Since 
March 2020 the primary source of documentation for this software has been Chapter 3 of Seccull's 
PhD Thesis titled ["Revealing Refractory Materials on Trans-Neptunian Objects and Centaurs via 
Reflectance Spectroscopy"](https://pure.qub.ac.uk/en/studentTheses/revealing-refractory-materials-on-trans-neptunian-objects-and-cen).
The thesis text will be released from embargo in the summer of 2022, but better documentation for 
MOTES is expected to be provided before then. All changes listed here relay developments to MOTES 
following 2021-03-01.

### Added
- This changelog.
- `if name == main:` has been introduced to `motes.py` so the MOTES functions may be called within 
  another python script, or within a python session.
- Demo data for FORS2, GMOS, and X-Shooter, all bright targets.
- Scale bars on the figures now allow users to adjust the cut in the plotted images.

### Changed
- The entire package has been made modular and processes all spectra from different instruments in 
  the same way. The process generally follows that defined for X-Shooter spectra in GME.
- `harvester.py` has been created (modified from the version that handled FORS2 spectra) and is 
  capable of unpacking data and relevant FITS header parameters from any FITS formatted 
  spectrogram. `data_harvest()` calls a more specialised harvester function depending on which 
  instrument a spectrogram has come from. `harvester.py` has been designed to allow the development 
  and inclusion of new specialised harvester functions for other instruments that can extract and 
  convert the relevant data from other instruments into a format that MOTES can process.
- Functions `readparfile()` and `readregions()`, which respectively read `motesparams.txt` and 
  `reg.txt` have been moved to their own module called `startup.py`.
- Controls for how the data frames are initially sliced have been removed from `motesparams.txt` 
  entirely; instead the user must define the number of rows to remove from the spatial axis of the 
  spectrum, and the wavelength range of the spectrum for each input file in `reg.txt`
- `motesparams.txt` must now be present in the same directory as the input data files. This has been 
  done to minimise the need to change the parameters for different setups between testing runs. It 
  should also make life easier for users who may be using MOTES with multiple datasets requiring 
  different parameters.
- The main script that was `forsextraction() or `xshooextraction()` has been incorporated into 
  `motes.py`. The `motes()` function is now the main function of the package. Sky subtraction 
  processes and cosmic ray handling have respectively been moved to their own functions within 
  `motes.py`; this is to reflect the fact that they are optional and are not required for motes to 
  run.
- Small rearrangements in the order of the output file header metadata have been done.
- Data frames are now almost always stored in a dedicated dictionary referred to as framedict or 
  some similar name depending on which function it's called from.
- Similarly to data frames, axes and parameters related to them are stored in axesdict or similar.
- All references to GME changed to MOTES; this includes references within the Python scripts 
  themselves and within comments.
- Version number has been dropped from v0.9.0 to v0.1.0, to better reflect reality (I'm also dumb 
  and should have read up on semantic versioning before sticking a number on this, ha!).
- Output frames like skymod, skybins, skyextractionlims and crmask are only added to the frame 
  dictionary if those sections of code are run. This may fix a bug where spectrograms that are 
  already skysubtracted are extracted but error out during saving.
- Output files from MOTES are prepended with 'm' rather than having '_GME_' included in the middle 
  of the filename.
- MOTES no longer walks through an ESO-like directory tree to find input files. All input FITS 
  files are expected to be present in the current working directory, and ordered in the same way as 
  their associated regions in `reg.txt`.

