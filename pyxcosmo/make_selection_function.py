#---------------------- Dedicated Function: make_selection_function
#  for python equivalent to readcol.pro see:
#  http://www.astrobetter.com/wiki/tiki-index.php?page=Python+Switchers+Guide
#---------------------- -----------------------
import numpy as np
from fitsutils import mrdfits
from fitsutils import mwrfits

def make_selection_function(template, background, out=None, verbose=False):
    """
    Create a table of probability of detection function of the count rate and the core radius, 
    for a given background level.

    Algorithm
    ---------
    1) It reads an ascii file containing the array of background value for which the table has already been computed,  
    with the corresponding 'middlefix'.
    2) Check that the given background is in the array of precomputed background.
    3) Read a the two tables of probabilty versus count rate and core radius for the background below and above the given background.
    4) Interpolate between these two tables. 
  
    Arguments
    ---------
    template: string, the template of the filename of the ascii file list of backgrounds.
    background: scalar, number : the value of the background.
 
    Keyword argument
    ---------
   verbose: bool, flag, if set print some messages.
   out: if it exists, name of the output file

    Returns
    -------
    A dictionary of arrays, with:
         the probability (array 2d, n1, n2)
         the count rate (array 1d, n1)
         the core radius (array 1d, n2)
         the background (scalar)


    Update:
    ------
    3 Oct 2014 created by Rene Gastaud.
    2 April 2015 add write fits keyword out

    """
    if (verbose): print ' function make_selection_function'
    # Check the known values
    mdtype={'names': ('value', 'name'),'formats': ('f4', 'S4')}
    back = np.loadtxt(template+'_knownback.txt', delimiter=' ', dtype=mdtype)

    if (background > back['value'].max()):
        print 'make_selection_function error background too big', background, back['value'].max()
        return -1

    if (background < back['value'].min()):
        print 'make_selection_function error background too small', background, back['value'].min()
        return -2


    i1 = np.where(back['value'] <=  background)[0].max()
    i2 = np.where(back['value'] >=  background)[0].min()

    name = template+'_'+back['name'][i1]+'.fits'
    name
    tab = mrdfits(name,1)
    result1 = tab['PROBA']


    name = template+'_'+back['name'][i2]+'.fits'
    name
    tab = mrdfits(name,1)
    result2 = tab['PROBA']

    if(i1==i2):
        result=result1
    else:
        dx = back['value'][i2]-back['value'][i1]
        x1 = (back['value'][i2]-background)/dx
        x2 = (background-back['value'][i1])/dx
        result = result1*x1 + result2*x2

    n_rcore = tab['RCORE'].size
    n_countr = tab['COUNTR'].size/n_rcore
    countr = tab['COUNTR'][0:n_countr]

    # trick for mwrfits
    background = np.array(background)
    str_prob = {'BACKGROUND':background,'RCORE':tab['RCORE'], 'COUNTR':countr, 'PROBA':result}

    if not(out is None): mwrfits(str_prob, out)

    return str_prob
