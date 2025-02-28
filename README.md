# MOTES

![Licence](https://img.shields.io/github/license/tseccull/motes)

Version 1.0.0

Last Updated: 2025-02-28 - T. Seccull

Developers: Tom Seccull & Dominik Kiersz

MOTES stands for the Modular and Optimized Tracer and Extractor of Spectra. 
This is a python package for tracing and extracting spectra of faint 
astronomical point sources from 2D spectrograms. MOTES can optionally also
sky subtract 2D spectra prior to extraction. 

Citation: We will have a Zenodo DOI soon, and eventually a publication.

# Project Status
Currently at v1.0.0. Extraction of longslit spectra observed by GMOS-N or 
GMOS-S and reduced by DRAGONS is fully supported. Other instruments with 
harvester functions (X-Shooter, FORS2, FLOYDS) still need to be grandfathered
into this latest version of the codebase. Documentation is being written. 
Watch this space! We aim to submit an accompanying article to one of the AAS 
journals, likely PASP.

# Supported Instruments/Data
[GMOS-N](https://www.gemini.edu/instrumentation/gmos) longslit spectra reduced with [DRAGONS](https://doi.org/10.3847/2515-5172/ad0044). Spectra reduced with IRAF are no longer supported.
[GMOS-S](https://www.gemini.edu/instrumentation/gmos) longslit spectra reduced with [DRAGONS](https://doi.org/10.3847/2515-5172/ad0044). Spectra reduced with IRAF are no longer supported.

# Installation
Python 3 is recommended for running MOTES. We intend to make it installable 
through pip. Watch this space.

# Usage
readthedocs documentation will be written for MOTES shortly. We also expect
to provide demo data and tutorials to help newcomers learn about MOTES. 

# Contributing
Procedures for users to supply new functionality, in particular new harvesters 
for different instruments, will be instated alongside documentation.

# License
[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

# References
MOTES relies on a number of Python packages that deserve recognition. If you
use MOTES in your work, please also be sure to cite the following:
[Astropy - Astropy Collaboration et al. 2013, A&A, 588, A33](https://doi.org/10.1051/0004-6361/201322068)
[Matplotlib - Hunter 2007, CSE, 9, 90](https://doi.org/10.1109/MCSE.2007.55)
[NumPy - Harris et al. 2020, Nature, 585, 357](https://doi.org/10.1038/s41586-020-2649-2)
[SciPy - Virtanen et al. 2020, NatMe, 17, 261](https://doi.org/10.1038/s41592-019-0686-2)

The optimal extraction method used by MOTES was first described by [Horne 1986, PASP, 98, 609](https://doi.org/10.1086/131801). 
Please also cite this article if you use MOTES to process your spectra.
