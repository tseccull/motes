# Change Log
All notable changes to this package will be documented here. 
This changelog follows the format described [here](https://keepachangelog.com/en/1.0.0/). [Semantic Versioning](https://semver.org/) is followed.

## 0.4.3 2023-04-21
Updates by T. Seccull

### Changed
* A bug was fixed in the `optimal_extraction()` function in common.py related to the enforcement of positivity required for the spatial profile when performing optimal extraction as described by Horne (1986). Previously this was misinterpreted as setting the both the spatial profile and its associated data to zero where data pixels were negative. Instead, the enforcement of positivity applies only the the spatial profile. Because, however, the spatial profile in MOTES is a model Moffat PSF, it has no negative values to begin with and positivity therfore does not need to be enforced. The effect of this change will be minimal for spectra of bright targets that have few negative pixels in their spatial profile. This update prevents overestimation of counts in faint regions of the spectrum with higher numbers of negative pixels. 

## 0.4.2 2023-03-05
Updates by D.Kiersz

### Changed
* `inputs` folder to store fits files, change output file format name to include 1D.
* Conversion to docstrings adhering to Google's style guide: https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings, added to all functions. This is the first iteration of in-line documention that will need refining.
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
* Installation environment for Anaconda and dependencies with instructions. This is the first step for contenarising the application.

### Changed
* Updated README.md, including a spelling mistake.
* In Removed the following line `sys.stdout.write('     Extraction of ' + direc + ' could not be completed.\n')`.
* In `startup.py` use a list-comprehension instead of a loop.

## 0.4.0 2022-07-27
Updates by T. Seccull

### Added
- An optimal extraction routine developed by Dominik Kiersz has 
now been fully integrated into MOTES. Some changes have been made to the
original design by TS to increase the efficiency of the process, and ensure that
is compatible with data from multiple instruments. It still operates on the 
original fundamental principles developed by DAK, however. MOTES now produces
both aperture, and optimally extracted spectra by default.

### Changed
- The comments related to many functions have been updated to more clearly 
describe their inputs, outputs, and purpose. This work is ongoing along with 
removing deprecated sections of code.
- The GMOS harvester function has been updated to handle the detector chip gaps 
more efficiently.
- There have been multiple changes to how the extraction limits are defined to 
account for the new optimal extraction routine. The main one being that the data 
is no longer supersampled prior to being extracted in the aperture method.  

## 0.3.0 2022-04-06
Updates by T. Seccull

### Added
- Harvester function has been added for the FLOYDS spectrographs on the Las Cumbres 
Observatory 2 m telescopes.

### Changed
- Astroscrappy has been run on the GMOS demo data to remove cosmic rays.
- Demo data has been removed from the github repository for now as the GMOS frames 
are too large. Return of the FORS2 and X-Shooter demo data frames to the repository 
should be fine in the near future.

### Deprecated
- All cosmic ray handling functions and lines have been commented 
out pending removal. This is because there are much better software tools, (e.g. 
Astroscrappy) that can be used for removing and replacing cosmic rays. Testing of 
the software on the demo data shows that it is fully functional without the cosmic 
ray handling sections. These regions will need to be tidied up, however.


## 0.2.1 2021-03-26

### Changed
- Colormap for plots updated to inferno. Using a perceptually uniform sequential
colormap really makes a difference in clearly defining structure in the 
spectrum.

## 0.2.0 - 2021-03-19
Updates by T. Seccull

### Changed
- Simplified the way the GMOS harvester handles incoming qual frame and converts
it so good pixels=1 and bad pixels=0, while ensuring the chip gaps aren't
flagged as bad pixels.
- Overhauled cosmic ray handling so CRs and badpixels can be replaced by the
input quality frame alone without the need for motes to use its, frankly 
piss-poor, CR detection routine. If software like astroscrappy is used to 
create a cosmic ray mask, it needs to be combined with the spectrogram's
input quality mask for motes to recognise the detected CRs as bad.
- In motesparams.txt, the -MASK_CR keyword is now -IDENT_CR and is specifically
used to refer to the motes CR identification routine independently of the
 -REPLACE_CR keyword. Setting either -IDENT_CR or -REPLACE_CR to True (1)
results in an updated cr_handlin() function in motes.py being called.

### Fixed
- Fixed incorrect header parameter dictionary call in save_fits() when 
recording the wavelength unit to the header of the output file.
- When determining datascale in motes(), the base 10 log of the absolute value 
of the median collapsed profile is now calculated. Previously this would fail 
if the median value turned out to be negative.
- Fixed outdated method of calling common.get_bins_output() in the case where 
CR masking is not done. This has now been updated in line with other calls to 
this function from motes().
- Fixed a bug where cosmic ray replacement would completely fail on spatial
pixel columns with no valid data (e.g. all nans or zeros), and would wrongly 
propagate bad pixel flags to adjacent columns with good data. common.get_bins()
now simply ignores columns with all bad data, as they are pretty much 
unrecoverable. 

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
different instruments in the same way. The process generally follows that 
defined for X-Shooter spectra in GME.
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

