#!/usr/bin/env python
#
#---------------------- -----------------------

import numpy as np

"""
THERE ARE TWO CAUSES OF DISCREPANCY BETWEEN US AND ANREA

1: IN THE NORMALIZATION EXPRESSION 10**norm_mt ANDREA USES 0.45961403
   WHEREAS ASpix ROUNDS IT TO 0.46; IT MAKES SOME DIFFERENCE

2: ANDREA'S COMPUTES MSTAR AS 1d14*10.^param.norm_mt/0.7*cosmo.h
   WHERE 0.7 IS h FROM ARNAUD 05 AND cosmo.h = 0.72 IS OUR h

   THE CORRECT EXPRESSION HOWEVER IS 1d14*10.^param.norm_mt * 0.7 / cosmo.h
   WHICH IS THE ONE THAT IS USED HERE

SO, USE AT LEAST 0.45961403 INSTEAD OF 0.46

"""

def get_Mstar_Arnaud05(h):
    """
    Gives the normalization of the MT relation from Arnaud et al.05
    for delta=200, pertaining to the kT>3.5 sample and with a pivot at 4 keV
    rescaled to an arbitrary value of h = H0 / (100 km/sec/Mpc)
    """

    # norm_mt = 0.46
    # Andrea uses the value below; since this is an exponent it matters
    norm_mt = 0.45961403

    h_Arnaud = 0.7

    Mstar = 1e14 * 10**norm_mt * h_Arnaud / h

    return Mstar

def mass_to_temperature(mass, redshift, cosmo, evol, param = None, verbose=False):


    """compute the temperature function of the mass,
        using X-ray scaling relations :
    temperature = Lstar*(mass/Mstar)^pow * (1+z)^zevol * (H(z)/H0)^hevol

    INPUT:
            mass
            redshift
            cosmo: cosmological parameter structure
            evol: evolution structure produced by mk_evol.pro mt_struct :

 OUTPUT: 
      temperature vector

 BEWARE:
      Default values of constant are taken from ASpix

 SEE ALSO:
      check_param and mk_mt
  """

    ###############################################################
    # CAREFUL CAREFUL CAREFUL
    #
    # THE NORMALIZATION OF THE RELATION 10*0.46=2.88 COMES FROM 
    # THE MT RELATION FOR THE kT > 3.5 CLUSTER SAMPLE OF ARNAUD 05 AT M200
    #
    # THE AMPLITUDE OF THAT RELATION IS 5.74 IN UNITS OF 10^14 MSUN
    # (NOT M / h !!!!!!)
    #
    # NOTE ALSO THAT h(Arnaud) = 0.7 AND THE PIVOT OF THE RELATION IS AT 5 KeV
    # 
    # FOR US h = 0.72 AND THE PIVOT IS AT 4 KeV
    #
    # THEREFORE THE NORMALIZATION IN UNITS OF 10^14 MSUN / h(ARNAUD)
    # WITH A PIVOT AT 4 KEv IS GIVEN BY
    # 5.74 x 0.7 x (4/5)^1.49 WHERE 1.49 IS THE EXPONENT OF THE POWER LAW
    # 5.74 x 0.7 x (4/5)^1.49 = 2.88 = 10^0.46, THE NORMALIZATION GIVEN IN ASpix
    # THEREFORE THE NORMALIZATION GIVEN IN ASpix IS FOR h = 0.7; WE MUST CONVERT
    # TO OUR h
    ####################################################################

    if cosmo is None:
        h = 0.72
    else:
        h = cosmo['H0'] / 100.0

    if param is None:
        Mstar = get_Mstar_Arnaud05(h)
        alpha_mt = 1.49
        zevol = 0.0
        hevol = -1.0
    else:
        if 'MSTAR' in param.keys():
            Mstar = param['MSTAR']
        else:
            Mstar = get_Mstar_Arnaud05(h)

        alpha_mt = param['POW_MT']
        zevol = param['ZEVOL_MT']
        hevol = param['HEVOL_MT']

    # Try to redo Andrea; 
    # First undo the Mstar = 1e14 * 10**norm_ml * h_Arnaud / h above
    Mstar *= (0.72 / 0.7)
    # Then fo the mstar:1d14*10.^param.norm_mt/0.7*cosmo.h in Andrea's code
    Mstar *= (0.72 / 0.7)
    # print "Our Mstar from mass_to_temperature", Mstar
    # print "mass_to_temperature BUGGY, JUST TO CHECK"

    Ez = np.interp(redshift, evol['Z'], evol['HC']) / cosmo['H0']

    temperature = 4.0 * (mass / Mstar * (1 + redshift)**(-zevol) * Ez**(-hevol))** (1.0 / alpha_mt)

    return temperature
