# Integrating an Instrument I/O Module

An instrument I/O module is a python script containing functions that
are critical for MOTES to correctly read and write data files produced 
by a particular astronomical spectrograph and/or data reduction
pipeline. Astronomical spectroscopic data frames from different
sources do not often share a common data format and metadata
structure. This means that `motesio.py`, the the primary module
used by MOTES to read and write data, must call secondary 
instrument-specific functions designed to handle datasets from each
unique source. These instrument-specific functions live in an
`instrumentio.py` module.

If you're a Python veteran and have some familiarity with astronomical
data, some of what follows may be familiar to you. It's hoped that
starting from the basics, however, will be helpful to anyone approaching
this kind of data processing/reduction work for the first time.

## Creating an Instrument I/O Module.

The following instructions will guide you through the process of 
creating your own instrument I/O module. The `instrumentio.py`
[template](instrumentio.py) also contains comments which will hopefully be useful to you.
For inspiration, you can also take a look at the instrument I/O 
modules that are already provided with MOTES in `/src/motes/io/`.
If you experience difficulties and would like some guidance, feel free
to open an [issue](https://github.com/tseccull/motes/issues) in the
MOTES GitHub repo.

To start with, you should copy the `instrumentio.py` template that's
provided with these instructions. By convention, the names of these
modules are the instrument's name followed by `io.py`. For example,
the I/O module for the GMOS spectrographs is `gmosio.py`. It isn't
mandatory to follow this convention when building an I/O module for
yourself, but if you intend to share it and request that it be added
to the base version of MOTES, you will be asked to follow it to ensure
consistency. Don't forget to alter the docstring at the top of your
copied `instrumentio.py` template and also the docstrings of each 
function. Sections of the docstring enclosed in square brackets should
be updated with relevant details about you and
the module you're creating.

### MOTES Style Conventions
A template instrument I/O module, `instrumentio.py`, has been provided
in the same directory as these instructions. In most cases, editing a
copy of `instrumentio.py` will be the easiest way to design an I/O
module for an instrument not already served by MOTES. When editing
the template, be mindful to match its style and formatting, especially
if you intend to share your module and request its addition to the
base version of MOTES.

- MOTES uses four spaces as indents, not tabs.
- Code line lengths should be kept below 80 characters wherever
possible, but line lengths up to 100 characters are permitted if the
alternative is difficult to read.
- Comment and docstring line lengths should be no longer than 72
characters.
- `snake_case` is the standard for MOTES, not `camelCase`.
- Keep variable names descriptive. Do not abbreviate them so much
that a newcomer cannot quickly figure out what they are.

## The Harvester Function
Instrument I/O modules contain two functions, the first and more
important of these is the harvester function. The harvester function's
job is to take the list of Header Data Units (HDUs) read in from the 
input .fits file by the `data_harvest()` function in `motesio.py`, 
retrieve required data frames and metadata from that HDU list, and 
reformat them into what MOTES expects as input.

### The List of Header Data Units
The harvester function receives only one argument from `data_harvest()`,
`input_fits_hdu`. This is the full list of HDUs extracted from the input
file. In many cases, this HDU List will contain all the required data
frames and header metadata to run MOTES. In some cases, there will be
extra HDUs/frames/headers that are not needed, and in other cases some
required data will either need to be calculated or generated in the
harvester function. As you read on you'll see what this means. 

### Accessing the HDU List
Astropy has extensive [documentation](https://docs.astropy.org/en/stable/io/fits/index.html) 
about how it reads/writes both header and image data from/to FITS files.
Here we'll just cover enough to ensure you can write a harvester
function capable of accessing the correct data frames and header data 
needed to run MOTES.

It's possible to get astropy to print much of the information needed to
access the different parts of an HDU list, including the names and 
indices of each HDU and the number of cards in each HDU's header. For
example, running the following commands in python on a file containing
the partially reduced GMOS spectrum of a star...

```python
import astropy.io.fits as fits
with fits.open("N20231106S0163_HD_54351_2D.fits") as input_file:
	print(input_file.info())
```

...results in the following output in the terminal...

```
Filename: N20231106S0163_HD_54351_2D.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU     188   ()      
  1  SCI           1 ImageHDU        52   (3138, 512)   float32   
  2  VAR           1 ImageHDU        52   (3138, 512)   float32   
  3  DQ            1 ImageHDU        54   (3138, 512)   int16 (rescales to uint16)   
  4  HISTORY       1 BinTableHDU     17   13R x 4C   [128A, 192A, 28A, 28A]   
  5  MDF           1 BinTableHDU     53   3R x 12C   [J, E, E, 20A, J, E, E, E, E, E, E, E]   
  6  PROVENANCE    1 BinTableHDU     17   4R x 4C   [28A, 128A, 128A, 128A]   
None
```

The name of an 
instrument's harvester function generally has the form 
`harvest_instrument()`.

## The Save Function

## Editing motesio.py and Finishing Up

  
