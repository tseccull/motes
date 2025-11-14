# MOTES: a Modular, Optimised Tracer and Extractor of Spectra 

![GitHub License](https://img.shields.io/github/license/tseccull/motes)

Version 1.0.0

Last Updated: 2024-12-05 - T. Seccull

Developers: Tom Seccull & Dominik Kiersz

MOTES is a python package used for tracing and extracting spectra of faint 
astronomical point sources from 2D spectrograms. MOTES is also capable of sky 
subtraction and , to a limited degree, cosmic ray rejection/masking. 

Citation: Please cite .... paper to be written .... if you use MOTES in your 
research. Need a Zenodo DO

# Project Status
Currently in final stages of development toward v1.0 and release. Dominik's 
optimal extraction routine is now formally included in MOTES. Documentation is 
still to be written, as well as an accompanying article that will likely be 
submitted to one of the AAS journals.

# Installation
Python 3 is recommended for running MOTES. We intend to make it installable through pip and/or conda. watch this space.

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

readthedocs documentation will be written for MOTES shortly. 

# Contributions
Procedures for users to supply new functionality, in particular new harvesters 
for different instruments.

# License
[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
