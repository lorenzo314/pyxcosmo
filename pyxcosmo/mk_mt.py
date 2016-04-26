#!/usr/bin/env python

#  16 March 2015  Rene Gastaud use is None
#---------------------- Dedicated Function: mk_mt
#---------------------- -----------------------

import numpy as np

from mk_lt import  broken_power_law_invert, broken_power_law

def mk_mt(cosmo, evol, MT_struct, mass=None, temp=None, redshift=0, verbose=False):
    """
    compute the mass function of the temperature, or the
    temperature function of the mass using X-ray scaling
    relations M = Mstar*(T/4keV)^pow * (1+z)^zevol * E(z)^hevol
           
    INPUT: 
    cosmo: cosmological parameter structure 
    evol: evolution structure produced by mk_evol.pro mt_struct :


    OPTIONAL INPUT: 
    mass : 
    temp: temperature (not temporary)
    redshift (default:z=0)
 

    OUTPUT: mass or temperature: vector

    called functions: 
    """
   
    if ((mass is None) and (temp is None)): 
        print '%MK_MT: Mass or Temperature required as input keyword' 
        return 0 

    #if ( (mass.size >0) and (temp.size >0)): 
    if (not(mass is None) and not(temp is None)): 
        print '%MK_MT: Mass and Temp keywords cannot be set simultaneously' 
        return 0 

    Mstar = MT_struct['MSTAR']
    Tpow = MT_struct['POW']
    broken1 = MT_struct['BROKEN1']
    broken2 = MT_struct['BROKEN2']
    broken3 = MT_struct['BROKEN3']
    broken4 = MT_struct['BROKEN4']
    if (verbose): print Mstar, Tpow, broken1, broken2, broken3, broken4
    powers  = [broken1, broken2, broken3, broken4, Tpow]
    tboundaries = np.array([1., 2., 3., 4.])

    zevol = MT_struct['ZEVOL']
    hevol = MT_struct['HEVOL']
    if (verbose): print zevol, hevol

    hz = np.interp(redshift, evol['Z'], evol['HC'])/cosmo['H0']

    if (mass==None): 
        if (verbose): print '%MK_MT:Temperature to Mass'
        temp = np.asarray(temp)
        mass = broken_power_law(temp, powers, tboundaries, Mstar, zevol, hevol, hz, redshift, verbose=verbose)
        return mass 

    if (temp==None): 
        if (verbose): print '%MK_MT: Mass to Temperature' 
        mboundaries = broken_power_law(tboundaries, powers, tboundaries, Mstar, zevol, hevol, hz, redshift)
        temperature = broken_power_law_invert(mass, powers, mboundaries, tboundaries, Mstar, zevol, hevol, hz, redshift, verbose=verbose)
        return temperature

    return -1
