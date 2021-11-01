# MOTES: a Modular, Optimised Tracer and Extractor of Spectra 

**Version:** 0.3.1_dev

**Last Updated:** 2021-11-01 - D. Kiersz

**Developers:** Tom Seccull (Lead) & Dominik Kiersz

MOTES is a Python package used for tracing and extracting spectra of faint 
astronomical point sources from 2D spectrograms. MOTES is also capable of sky 
subtraction and cosmic ray masking and subsequent rejection. 

**Citation:** Please cite *(this paper)* if you use MOTES in your 
research. 

## Project Status
Currently in final stages of development toward v1.0 and release. Documentation is still to be written, as well as an accompanying article 
that will likely be submitted to one of the AAS journals (MNRAS maybe?)

EFOSC2 harvester needs to be developed into the new standard with the test data provided by Dominik (from the P106 observing run).

## Installation
Python 3 is recommended for running MOTES. It is compatible with Python 2.7, but is also noticably slower.

For dependencies, installing an Anaconda enviroment from a `.yml` file will suffice in most cases. I've uploaded one as an example and the following terminal command can install the enviroment locally.

conda env create --name motes --file=PATH/to/motes.yml

## Usage
PDF Manual currently in construction on [Overleaf](https://www.overleaf.com/project/604614a9482e712b1ea2fbbd). 

## Contributions
Procedures for users to supply new functionality, in particular new harvesters 
for different instruments.

## License
[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
