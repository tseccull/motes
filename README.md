# MOTES: a Modular, Optimised Tracer and Extractor of Spectra 

[![Motes CI](https://github.com/tseccull/motes/actions/workflows/motes.yml/badge.svg)](https://github.com/tseccull/motes/actions/workflows/motes-cicd.yml)
![Licence](https://img.shields.io/github/license/tseccull/motes)

Version 0.4.7

Last Updated: 2024-02-25 - D. Kiersz

Developers: Tom Seccull & Dominik Kiersz

MOTES is a python package used for tracing and extracting spectra of faint 
astronomical point sources from 2D spectrograms. MOTES is also capable of sky 
subtraction and , to a limited degree, cosmic ray rejection/masking. 

Citation: Please cite .... paper to be written .... if you use MOTES in your 
research. Need a Zenodo DOI?

# Project Status
Currently in final stages of development toward v1.0 and release. Dominik's 
optimal extraction routine is now formally included in MOTES. Documentation is 
still to be written, as well as an accompanying article that will likely be 
submitted to one of the AAS journals.

# Installation
Python 3 is recommended for running MOTES. It is compatible with Python 2.7, but is also noticably slower.

## Anaconda

Installing an Anaconda enviroment from a `.yml` and running MOTES from the environment file will suffice in most cases. 

You can install the Anaconda pseudo-distribution for Python from [here](https://www.anaconda.com/).

Provided you are in the root directory of this repository, execute the first command to create the environment and install dependencies automatically. The second command will activate the environment.

```shell
conda env create --name motes --file=./conda/motes.yml
conda activate motes
```

MOTES can then be executed from `motes.py` in the root directory.

# Usage

PDF Manual currently in construction on [Overleaf](https://www.overleaf.com/project/604614a9482e712b1ea2fbbd). 

# Contributions
Procedures for users to supply new functionality, in particular new harvesters 
for different instruments.

# License
[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
