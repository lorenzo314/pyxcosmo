# IPython log file

import pytest
import os
from pyxcosmo import growth_a_ode
from fitsutils import mrdfits

def test_growth_a_ode(verbose=False):

    threshold = 1e-7
    if (verbose):print 'test_growth_a_ode', threshold
    # read the input
    mydir = os.environ['TESTDATADIR']
    fitsname = mydir+'growth_a_ode_in_out.fits'
    cosmo    = mrdfits(fitsname,'cosmo')
    if (verbose): print 'cosmo', cosmo
    aa    = mrdfits(fitsname,'aa')
    if (verbose): print 'aa',  aa.shape, aa.min(), aa.max()
    reference = mrdfits(fitsname,'python_out')
    if (verbose):print 'reference keys', reference.keys()
    # call the function
    output_py = growth_a_ode(cosmo, aa)
    # check the output
    diff = output_py-reference['SORTIE_ODE_PY_3']
    if (verbose):print 'diff', diff.min(), diff.max()
    assert(diff.min() > (-threshold) )
    assert(diff.max() <  threshold)
    return

