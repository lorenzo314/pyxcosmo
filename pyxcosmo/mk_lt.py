#!/usr/bin/env python
#  12 March 2015  R Gastaud : get_scalar_rg
#   now woks with broken a scalar or a numpy array

#---------------------- Dedicated Function: mk_lt
#---------------------- -----------------------

import numpy as np

def get_scalar_rg(x):
    if isinstance( x, np.ndarray):
        y = x[0]
    else:
         y = x
    return x
 

def broken_power_law(x, powers, xboundaries,  ystar, zevol, hevol, hz, z, verbose=False):
    #nn = len(x)  does not work on numpy array of one element
    nn = x.size
    ni = len(powers)
    nj = len(xboundaries)
    if ((ni-1) != nj):
        print "dimension powers boundaries not consistent"
        return -1
    #
    Tpow = np.tile(powers[ni-1], nn)
    norm = np.ones(nn)
    i = ni-1
    #print 'loop',i,   powers[i]
    #
    i = ni-2  # broken4
    ii = np.where(x < xboundaries[i])
    ii = ii[0]
    if (len(ii) > 0): Tpow[ii] = powers[i]
    if (verbose): print 'loop',i,  xboundaries[i],  powers[i]
    #
    i = ni-3
    mynorm=1.
    #print 'mynorm', mynorm
    for i in xrange(ni-3, -1, -1): 
        #print 'loop',i,  xboundaries[i],'/',  boundaries[ni-2], powers[i+1],'-',  powers[i]
        ii = np.where(x <  xboundaries[i])
        ii = ii[0]
        if (len(ii) > 0): 
            Tpow[ii] =  powers[i]
            mynorm =( xboundaries[i]/ xboundaries[ni-2])**( powers[i+1] - powers[i])*mynorm
            if (verbose): print 'mynorm', mynorm
            norm[ii] = mynorm
        #
    y = norm * ystar * (x/4.)**Tpow * (1+z)**zevol * hz**hevol
    # norm * Lstar * (Temp/4.)^Tpow * (1+z)^zevol * hz^hevol * tz^tevol
    return y



def  broken_power_law_invert(y, powers, yboundaries, xboundaries,  ystar, zevol, hevol, hz, z, verbose=False):
    #nn = len(y) does not work on numpy array of one element
    nn = y.size
    ni = len(powers)
    nj = len(xboundaries)
    if ((ni-1) != nj):
        print "dimension powers and xboundaries are not consistent"
        return -1
    nk = len(yboundaries)
    if ( nk != nj):
        print "dimension xboundaries and yboundaries are not consistent"
        return -1
 
    #
    if (verbose): print ' debug ystar, zevol, hevol, hz, z', ystar, zevol, hevol, hz, z
    Tpow = np.tile(powers[ni-1], nn)
    norm = np.ones(nn)
    i = ni-1
    #
    i = ni-2  # broken4
    ii = np.where(y < yboundaries[i])
    ii = ii[0]
    if (len(ii) > 0): Tpow[ii] = powers[i]
    #
    i = ni-3
    mynorm=1.
    #print 'mynorm', mynorm
    for i in xrange(ni-3, -1, -1): 
        if (verbose): print 'loop',i,'xboundaries',  xboundaries[i],'/',  xboundaries[ni-2],'powers', powers[i+1],'-',  powers[i]
        ii = np.where(y <  yboundaries[i])
        ii = ii[0]
        if (len(ii) > 0): 
            Tpow[ii] =  powers[i]
            mynorm =( xboundaries[i]/ xboundaries[ni-2])**( powers[i+1] - powers[i])*mynorm
            if (verbose): print 'mynorm', mynorm
            norm[ii] = mynorm
        #
    # 4.*(lum/Lstar/norm * (1+z)^(-zevol)*hz^(-hevol))^(1./Tpow)*tz^(-tevol) 
    x = 4.*(y/ystar/norm * (1+z)**(-zevol) * hz**(-hevol))**(1./Tpow)
    return  x

def mk_lt(cosmo, evol, LT_struct, lum=[], temp=[], redshift=0, verbose=False):


    """compute the luminosity function of the temperature, or the
    temperature function of the luminosity using X-ray scaling relations 
    luminosity = Lstar*(T/4keV)^pow * (1+z)^zevol * (H(z)/H0)^hevol

    INPUT:
             cosmo: cosmological parameter structure 
             evol: evolution structure produced by mk_evol.pro mt_struct :
             lt_struct :

 OPTIONAL INPUT: 
     luminosity: 
     temp: temperature (not temporary)
     redshift (default:z=0)
 

 OUTPUT: 
      luminosity or temperature: vector

 Called Functions: 
      broken_power_law, broken_power_law_invert

Remark:
   evol not wanted here!

  """
   
    if ((len(lum)==0) and (len(temp)==0)): 
        print '%MK_LT:  Luminosity or Temperature required as input keyword' 
        return 0 
    # 
    if ((len(lum)>0) and (len(temp)>0)): 
        print '%MK_LT:  Luminosity and Temperatrure keywords cannot be set simultaneously' 
        return 0 
    #

    if not(evol):
         mk_evol, cosmo, evol

    Lstar   = get_scalar_rg(LT_struct['LSTAR'])
    Tpow    = get_scalar_rg(LT_struct['POW'])
    # if read from fits, need  [0], if computed by check_param, don't need
    broken1 = get_scalar_rg( LT_struct['BROKEN1'])
    broken2 = get_scalar_rg( LT_struct['BROKEN2'])
    broken3 = get_scalar_rg( LT_struct['BROKEN3'])
    broken4 = get_scalar_rg( LT_struct['BROKEN4'])
    powers  = [broken1, broken2, broken3, broken4, Tpow]
    tboundaries = np.array([1., 2., 3., 4.])

    zevol = get_scalar_rg(LT_struct['ZEVOL'])
    hevol = get_scalar_rg(LT_struct['HEVOL'])
    tevol = get_scalar_rg(LT_struct['TEVOL'])
    
    hz = np.interp(redshift, evol['Z'], evol['HC'])/cosmo['H0']
    tz = (evol['AGE'] - np.interp(redshift, evol['Z'], evol['T'])) / evol['AGE']
    if (verbose): print 'hz=', hz, 'hevol=', hevol, 'tz=', tz, 'tevol=', tevol
    # 
    if ((len(lum)==0) and (len(temp)>0)): 
        if (verbose): print '%MK_MT:Temperature to Luminosity'
        luminosity = broken_power_law(temp, powers, tboundaries, Lstar, zevol, hevol, hz, redshift, verbose=verbose)
        # if tevol == 0 next line does nothing
        # anything**0 = 1
        luminosity = luminosity*tz**tevol

        return luminosity 
    # 
    if ((len(lum)>0) and(len(temp)==0)): 
        if (verbose): print '%MK_MT:  Luminosity to Temperature' 
        lboundaries = broken_power_law(tboundaries, powers, tboundaries, Lstar, zevol, hevol, hz, redshift, verbose=verbose)
        temperature = broken_power_law_invert(lum, powers, lboundaries, tboundaries, Lstar, zevol, hevol, hz, redshift, verbose=verbose)
        # if tevol == 0 next line does nothing
        # anything**0 = 1
        temperature =  temperature/tz**tevol

        return temperature

    return -1




##########
