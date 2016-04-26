#!/usr/bin/env python

#  16 March 2015  Rene Gastaud use is None
#---------------------- Dedicated Function: mk_mt
#---------------------- -----------------------

import numpy as np

from mk_lt import  broken_power_law_invert, broken_power_law

def mk_ml(cosmo, evol, ML_struct, mass=None, lum=None, redshift=0, verbose=False):
    """
    compute the mass function of the luminosity or the
    luminosity function of the mass using X-ray scaling
    relations M = Mstar*(L/Lstar)^pow * (1+z)^zevol * E(z)^hevol
           
    INPUT:
    cosmo: cosmological parameter structure 
    evol: evolution structure produced by mk_evol.pro
    ml_struct :

    OPTIONAL INPUT: 
    mass: 
    lum: luminosity (not temporary)
    redshift (default:z=0)

    OUTPUT: mass or luminosity: vector

    called functions: 
    """

    if ((mass is None) and (lum is None)): 
        print '%MK_ML: Mass or Luminosity required as input keyword' 
        return 0 

    if (not(mass is None) and not(lum is None)): 
        print '%MK_ML: Mass and Luminosity keywords cannot be set simultaneously' 
        return 0 

    Mstar = ML_struct['MSTAR']
    Lpow = ML_struct['POW']
    broken1 = ML_struct['BROKEN1']
    broken2 = ML_struct['BROKEN2']
    broken3 = ML_struct['BROKEN3']
    broken4 = ML_struct['BROKEN4']
    if (verbose): print Mstar, Lpow, broken1, broken2, broken3, broken4
    powers  = [broken1, broken2, broken3, broken4, Lpow]
    lboundaries = np.array([1., 2., 3., 4.])

    zevol = ML_struct['ZEVOL']
    hevol = ML_struct['HEVOL']
    if (verbose): print zevol, hevol

    hz = np.interp(redshift, evol['Z'], evol['HC'])/cosmo['H0']

    if (mass==None):
        if (verbose): print '%MK_ML: Luminosity to Mass'
        lum = np.asarray(lum)
        mass = broken_power_law(lum, powers, lboundaries, Mstar, zevol, hevol, hz, redshift, verbose=verbose)
        return mass 

    if (lum==None):
        if (verbose): print '%MK_ML: Mass to luminosity' 
        mboundaries = broken_power_law(lboundaries, powers, lboundaries, Mstar, zevol, hevol, hz, redshift)
        luminosity = broken_power_law_invert(mass, powers, mboundaries, lboundaries, Mstar, zevol, hevol, hz, redshift, verbose=verbose)
        return luminosity

    return -1
