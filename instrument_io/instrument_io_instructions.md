# Integrating an Instrument I/O Module

An instrument I/O module is a python script containing functions that
are critical for MOTES to correctly read and write data files produced 
by a particular astronomical spectrograph and/or data reduction
pipeline. Astronomical spectroscopic data frames from different
sources do not often share a common data format and metadata
structure. This means that `motesio.py`, the the primary module
used by MOTES to read and write data, must call secondary 
instrument-specific functions designed to handle datasets from each
unique source. These instrument-specific functions live in the
`instrumentio.py` module.

## Creating an Instrument I/O Module.

The following instructions will guide you through the process of 
creating your own instrument I/O module. The `instrumentio.py`
template also contains comments which will hopefully be useful to you.
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
consistency. Don't forget to alter the docstring at the top of your copied 
`instrumentio.py` template. Sections of the docstring enclosed in
square brackets should be adjusted to match

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

## The Save Function

## Editing motesio.py and Finishing Up

  
