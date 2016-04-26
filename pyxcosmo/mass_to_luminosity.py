#!/usr/bin/env python
#
#---------------------- -----------------------

import numpy as np
from mass_to_temperature import get_Mstar_Arnaud05

"""
The ASpix values for ML come from combining these two relations

M = 10e14 * 10**0.46 * (T/4keV) ** 1.49
L = 10e44 * 10**0.40 * (T/4keV) ** 2.89

1.49/2.89 = 0.52 as in ASpix

They give

M = 10e14 * 10**0.46 * (L/1e44) ** 0.52 * ((1/10**0.4) ** 0.52)

((1/10**0.4) ** 0.52) * (10 ** 0.46) = 1.79

and 10 ** 0.25 = 1.78 (the ASpix normalization)

It is still necessary to properly take into account the factors of h
which differ from Arnaud 05 and us
"""

def get_Mstar_Arnaud05(h):
    # From ASpix
    norm_ml = 0.25

    h_Arnaud = 0.7

    Mstar = 1e14 * 10**norm_ml * h_Arnaud / h

    return Mstar

def mass_to_luminosity(mass, redshift, cosmo, evol, param = None, verbose=False):


    """compute the luminosity function of the mass,
        using X-ray scaling relations :
    luminosity = Lstar*(mass/Mstar)^pow * (1+z)^zevol * (H(z)/H0)^hevol

    INPUT:
            mass
            redshift
            cosmo: cosmological parameter structure
            evol: evolution structure produced by mk_evol.pro mt_struct :

 OUTPUT: 
      luminosity vector

 BEWARE:
      Default values of constant are taken from ASpix

 SEE ALSO:
      check_param and mk_lt
  """
    if cosmo is None:
        h = 0.72
    else:
        h = cosmo['H0'] / 100.0

    if param is None:
        Mstar = get_Mstar_Arnaud05(h)
        Lstar=1e44
        alpha_ml = 0.51557092
        zevol = 0.0
        hevol= -1.51557092
    else:
        if 'LSTAR' in param.keys():
            Lstar = param['LSTAR']
        else:
            Lstar = 1e44

        if 'MSTAR' in param.keys():
            Mstar = param['MSTAR']
        else:
            Mstar = get_Mstar_Arnaud05(h)

        alpha_ml = param['POW_ML']
        zevol = param['ZEVOL_ML']
        hevol = param['HEVOL_ML']

    Ez = np.interp(redshift, evol['Z'], evol['HC']) / cosmo['H0']

    # Try to redo Andrea; 
    # First undo the Mstar = 1e14 * 10**norm_ml * h_Arnaud / h above
    Mstar *= (0.72 / 0.7)
    # Then fo the mstar:1d14*10.^param.norm_mt/0.7*cosmo.h in Andrea's code
    Mstar *= (0.72 / 0.7)
    # print "Our Mstar from mass_to_luminosity", Mstar
    # print "mass_to_luminosity BUGGY, JUST TO CHECK"

    luminosity = Lstar * (mass / Mstar * (1+redshift)**(-zevol) * Ez**(-hevol)) ** (1.0 / alpha_ml)

    return luminosity
