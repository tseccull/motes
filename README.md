# MOTES: a Modular, Optimised Tracer and Extractor of Spectra 

Version: 0.2.0_dev

Last Updated: 2021-03-04 - T. Seccull

Developers: Tom Seccull (Lead) & Dominik Kiersz

MOTES is a Python package used for tracing and extracting spectra of faint 
astronomical point sources from 2D spectrograms. MOTES is also capable of sky 
subtraction and cosmic ray masking and subsequent rejection. 

Citation: Please cite .... paper to be written .... if you use MOTES in your 
research. Need a Zenodo DOI?

# Project Status
Currently in final stages of development toward v1.0 and release. Updates and 
inclusion of optimal extraction capability are to be included in the near 
future. Documentation is still to be written, as well as an accompanying article 
that will likely be submitted to one of the AAS journals.

This development branch is setup to transfer relevant methods and fixes from Dom's GME (locally called GME Plus) to Tom's MOTES master branch. Primarly to transfer supersampled optimal extraction,improvements to fitting methods where possible and optimilisations where possible. EFOSC2 harvester needs to be developed into the new standard with the test data provided by Dom (from the P106 observing run).

# Installation
Python 3 is recommended for running MOTES. It is compatible with Python 2.7, but is also noticably slower.

For dependencies, installing an Anaconda enviroment from a .yml file will suffice in most cases. I've uploaded one as an example and the following terminal command can install the enviroment locally.

conda env create --name motes --file=motes.yml

The above needs testing to ensure out-of-the-box compatibility. The rest of this section needs to be worked out.

# Usage
PDF Manual currently in construction on Overleaf.

https://www.overleaf.com/project/604614a9482e712b1ea2fbbd

# Contributions
Procedures for users to supply new functionality, in particular new harvesters 
for different instruments.

# License
[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
