###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import os
import sys

###############################################################################
# FUNCTIONS ///////////////////////////////////////////////////////////////// #
###############################################################################
# ////////////////////////////////////////////////////////////////////////////#
# Import parameters from motesparams.txt parameter file into a dictionary.
# INPUTS:  None
# OUTPUTS: pars - dictionary containing parameters read in from motesparams.txt
def read_parfile():

    sys.stdout.write(' >>> Reading in parameters from motesparams.txt. ')
    sys.stdout.flush()

    # Read in MOTES parameter file line by line and filter out the empty lines.
    with open('motesparams.txt') as parfile:
        parlines = parfile.read().splitlines()
        parlines = filter(None, parlines)

    # Flatten the 2D list of parameters and keywords into a 1D list where each
    # parameter's value follows its associated keyword.
    lumpy_param_list = [x.split('=') for x in parlines if x[0] == '-']
    flat_param_list = [y for x in lumpy_param_list for y in x]

    # Convert all numerical values in the parameter list to floats.
    paramlist = []
    for i in flat_param_list:
        try:
            paramlist.append(float(i))
        except ValueError:
            paramlist.append(i)

    # Assign parameters and their associated keywords to a dictionary.
    pars = dict(zip(paramlist[::2], paramlist[1::2]))

    sys.stdout.write('DONE.\n')

    return pars


# /////////////////////////////////////////////////////////////////////////// #
# Search for, and read in, reg.txt file in root working directory.
# Reads an input line from reg.txt and returns a list of integers defining the 
# boundaries of the region of the 2D data that contains the spectrum to be 
# extracted.
# INPUTS:  None
# OUTPUTS: intregion - a list, which for each file contains a list of integers that define the boundaries of the region of the 2D spectum that will be used for the extraction. For each file, the 
#                      first two integers are the number of pixel rows to remove from each end of the spatial axis of the spectrum, and the last two are the upper and lower wavelength bounds of the 
#                      region on the dispersion axis.
def read_regions():

    # Search for reg.txt and read in the list that it contains.
    if os.path.isfile('reg.txt'):
        sys.stdout.write(' >>> reg.txt file found. Reading reg.txt. ')
        sys.stdout.flush()

        with open('reg.txt') as reg:
            reg = reg.read().splitlines()
            intregion = [[int(lim) for lim in x.split(',')] for x in reg]

        sys.stdout.write('DONE.\n')

    # Complain and quit MOTES if reg.txt isn't found.
    else:
        sys.stdout.write(' >>> reg.txt file not found in root working directory.\n')
        sys.stdout.write('     Root working directory is:\n')
        sys.stdout.write('     ' + os.getcwd() + '\n')
        sys.stdout.write('     Please provide a reg.txt file to tell MOTES the spatial extent of the 2D spectra.\n')
        sys.stdout.write('     See the docs for further info.\n')
        sys.stdout.write('     Terminating MOTES.\n\n')
        exit()

    return intregion
