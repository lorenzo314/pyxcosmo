# IPython log file

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
import os.path
from fitsutils import mrdfits
from fitsutils import mwrfits

from icosmo_cosmo import xgen
# can be replaced by np.logspace

def make_mesindexcrtab(original, cr_min, cr_max, n_cr, verbose=False, save=True):

    """
    Compute the location of a vector of cr in a given table of cr
    The vector of cr is log linear, with 
    INPUT : 
         orignal :  file name of the input cr table
         cr_min : minimum of the vector of cr
         cr_max : maxium of  the vector of cr
         n_cr : number of points of the vector of cr
         outfile : the name of the output file
         save : save a file, default True
         verbose: default false
         
         
 OUTPUT : this function returns a list containing 
            1) the locations
            2) the original table
    """
    if (verbose):print original, cr_min, cr_max, n_cr, verbose, save

    #  compute the step log linear of the vector of cr
    dlog_cr = (math.log(cr_max) - math.log(cr_min))/(n_cr-1.)
    # compute the origin
    cr_start = cr_min*math.exp(-dlog_cr/2.)
    cr_stop = cr_max*math.exp(dlog_cr/2.)
    if (verbose): print 'step=', dlog_cr, ' start=', cr_start,' stop=', cr_stop 


    input_table = mrdfits(original,1)
    crtab = input_table['CR_PERL_TIMESD2']
    locations = (np.log(crtab)-math.log(cr_start))/dlog_cr

    vec_cr_bounds = xgen(cr_start, cr_stop, npoints=n_cr+1, logplot=True)
    # np.logspace(math.log10(cr_start), math.log10(cr_stop), num=n_cr+1, endpoint=True)


    result = {'INDEX_CR_PERL_TIMESD2':locations,
              'ORIGINAL_CR_PERL_TIMESD2':crtab, 'VEC_CR_BOUNDS':vec_cr_bounds}

    if (save): mwrfits(result, original+'.index')

    return result
