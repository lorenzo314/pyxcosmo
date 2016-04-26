#!/usr/bin/env python
# RG 20 October 2015 version 4 : scatter > 0 processed by compute_dndzdt
# RG 09 October 2015 version 3 : add gtype  , see mk_nm, mk_evol, cmp_dndzdldt, combine_dndzdldt
#   4 May 2015    RG add test on n_t version 2
#  13 March 2015  R Gastaud
#  13 March 2015  R Gastaud

import numpy as np
import matplotlib.pyplot as plt
import math
import pdb

from icosmo_cosmo import mk_evol
from icosmo_cosmo import convert_nm
from mk_mt import mk_mt
from deriv_idl import deriv_idl
from scipy.interpolate import interp1d
from mk_lt import mk_lt
from cmp_aperture_luminosity_fast import cmp_aperture_luminosity_fast

__version__ = 5

from int_tabulated import int_tabulated

# RG 20 October 2015 version 4 : for scatter > 0
def compute_dndzdt(nm, ScTab, temperature, temp_out):
    ### compute
    ## ratio
    n_t = temp_out.size
    nn = temperature.size
    
    ratio = np.reshape(temp_out, [n_t,1])/np.reshape(temperature, [1,nn])
    
    result = np.exp(-0.5*(np.log(ratio)/ScTab)**2)
    
    ##  norm_ttab
    mynorm = np.reshape(nm['DNDMDOMDZ'], [nn,1])
    norm_ttab = mynorm/temp_out
    norm_ttab =  norm_ttab.T
    #print 'transpose'
    
    ##  dndzdmdt
    dndzdmdt = norm_ttab/math.sqrt(2*math.pi)/ScTab*result
    
    #  dndzdt sum over mass m
    dndzdt = np.zeros(n_t)
    for i in range(0, n_t):
        dndzdt[i]=int_tabulated(np.log(nm['M']), nm['M']*dndzdmdt[i,:])
    
    return dndzdt

def cmp_dndzdldt(cosmo, nm, mt_struct, lt_struct, fluxmes_struct, model_struct, \
    tmin=0.3, tmax=30.0, l_min=1.0E40, l_max=8.0E48, gtype='ODE', verbose=False):
    """
    Compute dndzdldt
    History:
    4 May 2015    RG add test on n_t version 2
    13 March 2015  R Gastaud
    RG 09 October 2015 version 3 : add gtype  , see mk_nm, mk_evol, cmp_dndzdldt, combine_dndzdldt
    RG 18 October 2015 version 4 : mt_struct['SCATTER'] >= 0
    """

    z = nm['Z']
    evol = mk_evol(cosmo,ran_z=[0,0.01]+z,n_z=2, gtype=gtype)

    if (model_struct.has_key('MLIMZ_M')):
        #print 'not yet implemented'
        #return float('NaN')
        raise NameError('mk_evol model_struct.MLIMZ_M not yet implemented')

    nmconv200c=convert_nm(cosmo, evol, nm, overdensity=200, rhoc=True)
    nmconv500c = convert_nm(cosmo, evol, nm, overdensity=500, rhoc=True)

    mytype = ''.join(mt_struct['MTYPE'])
    if(verbose): print 'mytype:', mytype, ' ', mt_struct['MTYPE']
    if (mt_struct['MTYPE'] == 'm200'):
        temperature = mk_mt(cosmo, evol, mt_struct, mass=nmconv200c['M'],redshift=z)
    elif (mt_struct['MTYPE'] == 'm500'):
        temperature = mk_mt(cosmo, evol, mt_struct, mass=nmconv500c['M'],redshift=z)
    else:
        raise NameError('mk_evol unknown mytype='+mytype)

    #### now work on temperature 
    n_t = model_struct['N_T']
    if isinstance(n_t, np.ndarray): n_t = n_t[0] # 4 May 2015
    temp_out = np.logspace(math.log(tmin), math.log(tmax), num=n_t, base=math.exp(1.))

    if (mt_struct['SCATTER'] > 0):
        if verbose: print 'scatter', mt_struct['SCATTER']
        dndzdt = compute_dndzdt(nm, mt_struct['SCATTER'], temperature, temp_out)
    else:
        myderiv = deriv_idl(np.log(temperature), np.log(nm['M']))
        dndzdt = nm['DNDMDOMDZ'] * nm['M']/temperature * myderiv
        interpolator = interp1d(temperature, dndzdt)
        if (verbose):
            print 'temperature in minmax', temperature.min(), temperature.max()
            print 'temperature out minmax', temp_out.min(), temp_out.max()
        dndzdt = interpolator(temp_out)

    #### now luminosity
    lum = mk_lt(cosmo, evol, lt_struct, temp=temp_out, redshift=z)
    # remove mynorm in mk_lt

    n_l = model_struct['N_L']  # 4 May 2015
    if isinstance(n_l, np.ndarray): n_l = n_l[0]
    lum_out = np.logspace(math.log(l_min), math.log(l_max), num=n_l, base=math.exp(1.))

    sc_tab = lt_struct['SCATTER'] + lt_struct['SCATTER_EVOL'] * math.log(1+z)

    # n_l, n_t
    norm     = dndzdt.reshape(1, n_t)
    Ttab     = temp_out.reshape(1, n_t)
    LTab     = lum_out.reshape(n_l, 1)
    LmeanTab = lum.reshape(1, n_t)

    # RG 08/01/16: Scatter in luminosity
    rapport = np.log(LTab/LmeanTab)/sc_tab
    if verbose: rapport.shape ,  n_l, n_t
    # (180, 150)
    # Laplace Gauss law
    dndzdldt = norm/math.sqrt(2.*math.pi)/sc_tab * np.exp(-0.5*rapport**2)

    #### thresholding temperature too low
    if(model_struct['TLIM'] > 0):
        ii = np.where(temp_out < model_struct['TLIM'])
        if (ii[0].size > 0):
            dndzdldt[:, ii[0]]=0

    ###### finite aperture
    interpolator = interp1d(evol['Z'], evol['HC'])
    Ez = interpolator(z)/cosmo['H0']

    m200c = mk_mt(cosmo,evol,mt_struct,temp=temp_out,redshift=z)
    # remove mynorm in mk_mt

    # interpolator = interp1d(nmconv200c['M'], nmconv500c['M'])
    # m500c = interpolator(m200c)
    # Replaced by a simple linear relation from A Valotti 20 Nov 2015 BEWARE
    m500c = m200c*0.714 # following Arya

    rhoc_0_e4 =  2.7889537E7
    rhoc_z = rhoc_0_e4*(Ez*cosmo['H0'])**2. 

    r500c = ( m500c*(.75/math.pi/ 500./rhoc_z / cosmo['H'] ))**(1/3.)

    # Parametrize the xc parameter of the AB model see Pratt & Arnaud (2002) 
    xcparameter = model_struct['XC_0']*(1+z)**(model_struct['XC_H'])*Ez**(model_struct['XC_H'])

    # 4 May 2015
    #print 'fluxmes_struct keys', fluxmes_struct.keys()
    aperture = fluxmes_struct['APERTURE']
    if isinstance(fluxmes_struct['APERTURE'], np.ndarray):
        n_apertures = fluxmes_struct['APERTURE'].size
        if(n_apertures > 1):
            aperture_0 = fluxmes_struct['APERTURE'][0]
        else:
            aperture_0 = fluxmes_struct['APERTURE']
    else:
        aperture_0 = fluxmes_struct['APERTURE']
        n_apertures = 1

    if (aperture_0 > 0):
        # 30.
        interpolator = interp1d(evol['Z'], evol['DA'])
        Da = interpolator(z)

        # conversion in arcseconds
        r500angular = (r500c/Da)*(180/math.pi)*3600

        # BEWARE BUG here when aperture is a vector
        phys_aperture = aperture/r500angular 
        fasttable = fluxmes_struct['PROFTABLE']

        if (verbose): print 'fasttable', fasttable

        fluxratio = cmp_aperture_luminosity_fast(fasttable,phys_aperture, r500angular, xc=xcparameter)

        ###  not useful, fluxratio does not depend upon luminosity (why ??)
        fluxratiotab = np.tile(fluxratio, n_l).reshape(n_l, n_t)
    else:
        fluxratiotab = np.ones([n_l, n_t]) # bug corrected 20 oct 2015

    # The rcore model (already included in cmp_aperture_luminosity)
    rcoretab = xcparameter*r500c*1E3

    res = { 'DNDZDLDT' : dndzdldt,'TEMP' : temp_out, 'L' : lum_out,'Z' : z,'FLUXRATIOTAB' : fluxratiotab, 'RCORETAB' : rcoretab }

    return res
