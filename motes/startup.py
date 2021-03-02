###############################################################################
# IMPORT MODULES /////////////////////////////////////////////////////////////#
###############################################################################
import os
import sys

###############################################################################
# FUNCTIONS ///////////////////////////////////////////////////////////////// #
###############################################################################
# ////////////////////////////////////////////////////////////////////////////#
# IMPORT PARAMETERS FROM PARAMETER FILE INTO A DICTIONARY
def read_parfile():

    sys.stdout.write(' >>> Reading in parameters from motesparams.txt. ')
    sys.stdout.flush()

    # Read in MOTES parameter file line by line and filter out the empty lines.
    with open(sys.path[0] + '/motesparams.txt') as parfile:
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
# SEARCH FOR, AND READ IN, reg.txt FILE IN ROOT WORKING DIRECTORY

def read_regions():

    # Search for reg.txt and read in the list it contains.
    if os.path.isfile('reg.txt'):
        sys.stdout.write(' >>> reg.txt file found. Reading reg.txt. ')
        sys.stdout.flush()

        with open('reg.txt') as reg:
            reg = reg.read().splitlines()
            intregion = []
            for x in reg:
                lims = x.split(',')
                intregion.append([int(lims[0]), int(lims[1]), int(lims[2]), int(lims[3])])

        sys.stdout.write('DONE.\n')

    # Complain and quit if reg.txt isn't found.
    else:
        sys.stdout.write(' >>> reg.txt file not found in root working directory.\n')
        sys.stdout.write('     Root working directory is:\n')
        sys.stdout.write('     ' + os.getcwd() + '\n')
        sys.stdout.write('     Please provide a reg.txt file to tell MOTES the spatial extent of the 2D spectra.\n')
        sys.stdout.write('     See the docs for further info.\n')
        sys.stdout.write('     Terminating MOTES.\n\n')
        exit()

    return intregion
