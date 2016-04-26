#!/usr/bin/env python
#  15 March 2015  use is None
#  14 March 2015  kpiv and deltar2 version 0.3

import os.path
from fitsutils import mrdfits
from fitsutils import mwrfits

from scipy import constants
import math

import pdb
version=0.3

def rd_cosmo(h=None, omega_c=None, omegac_h2=None, omega_b=None,omegab_h2=None,
             omega_l=None, omega_r=None, w_l=None, w1_l=None, age=None, 
             n_pk= None, pk_feature=None, sigma8=None, deltar2=None, 
             kpiv=None, gamma=None, tau=None, verbose=False):

    """
    Purpose
    -------------
    Read cosmological parameters. The conventions of 
    Peacock (1999) or Kaiser (1996) are adopted

    Inputs
    ----------
    h   : Hubble constant in units 100 km/s/Mpc
    omega_c   : Mean Cold Matter Density, also called omega_dm
    omega_r   : Mean Radiation Density
    omega_m   : Mean Matter Density
    omega_b   : Mean baryon density
    omega_l   : Mean dark energy/lambda density
    w_l       : Dark energy equation of state
    w1_l      : varying w parametrised as: w(a)=w_l+w1_l(1-a)
    n         : Slope of primordial fluctuation power spectrum
    sigma8    : Normalization of matter power spectrum at z=0 
    tau       : Optical Depth

    Returns
    -------
    cosmo: dictionnary with all information

    """


    #construct cosmo structure

    geom={'H':0.7 ,'H0':0. ,'RH':0. ,'OMEGA_DM':0.3, 'OMEGA_B':0.,
          'OMEGA_M':0., 'OMEGA_L':0.7,'OMEGA_R':0., 
          'OMEGA':0., 'OMEGA_K':0.,'W_L':-1., 'W1_L':0.,
          'K':0., 'SQRTK':0.,'R0':0.}

    #Hubble constant
    if (h!=None): geom['H']=h
 
    #  Baryon Density
    if (omega_b != None) and  (omegab_h2 != None) :
        print "%RD_COSMO: only omega_b or omegab_h2 should be provided"
        print 'OMEGA_B', omega_b
        return None

    if not(omega_r is None):  
        geom['OMEGA_R'] = omega_r 
        if (verbose): print 'omega_r updated', omega_r

    if not(omega_b is None):  
        geom['OMEGA_B'] = omega_b 
        if (verbose): print 'omega_b updated', omega_b


    if not(omegab_h2 is None):  
        geom['OMEGAB_H2'] = omegab_h2 /geom['H']**2 
        if (verbose): print 'omegab_h2 updated', geom['OMEGAB_H2']


       #  Cold Darl Matter Density
    if not(omega_c is None) and  (omegac_h2 != None) :
        print "%RD_COSMO: only omega_c or omegac_h2 should be provided"
        print 'OMEGA_C', omega_c
        return None

    if not(omega_c is None):  
        geom['OMEGA_DM'] = omega_c  
        if (verbose): print 'omega_c = DM updated', geom['OMEGA_DM']

    if not(omegac_h2 is None):  
        geom['OMEGAC_H2'] = omegac_h2 /geom['H']**2 
        if (verbose): print 'omegac_h2 updated', geom['OMEGAC_H2']

    ##  Dark Energy
    if not(omega_l is None):  
        geom['OMEGA_L'] = omega_l  
        if (verbose): print 'omega_l updated', geom['OMEGA_L']
    if not(w_l is None):  
        geom['W_L'] = w_l  
        if (verbose): print 'w_l updated', geom['W_L']
    if not(w_l is None):  
        geom['W1_L'] = w1_l  
        if (verbose): print 'w1_l updated', geom['W1_L']

    ## derived geometry parameters

    geom['H0'] = 100.*geom['H']
    geom['RH'] = constants.c/1E3/geom['H0'] #Hubble radius [Mpc]

    # omega_m total matter density (z=0)
    #db.set_trace()
    geom['OMEGA_M'] = geom['OMEGA_DM'] + geom['OMEGA_B']

    # omega Total mass density  (z=0)
    geom['OMEGA'] = geom['OMEGA_M'] + geom['OMEGA_L'] + geom['OMEGA_R']
 
    # Global curvature (K)
    geom['OMEGA_K'] = 1. - geom['OMEGA']
   

    ## useful factors
    if (geom['OMEGA'] > 1): # close universe
        if (verbose): print 'omega > 1'
        geom['K'] = 1 
        geom['SQRTK'] = math.sqrt(geom['OMEGA']-1.) 

    if (geom['OMEGA'] == 1):  # flat universe
        if (verbose): print 'omega = 1'
        geom['K'] = 0 
        geom['SQRTK'] = 1.

    if (geom['OMEGA'] < 1): # open universe
        if (verbose): print 'omega < 1'
        geom['K'] = -1 
        geom['SQRTK'] = math.sqrt(1.-geom['OMEGA']) 

    # ;scale radius at z=0 [Mpc]
    #geom.r0=geom.rh/geom.sqrtk
    geom['R0'] = geom['RH']/geom['SQRTK'] 
                                
    ## CDM Power Spectrum
    PowerSpec = {'N':n_pk, 'TAU':tau, 'NORM': 1, 'SIGMA8':sigma8, 'GAMMA':gamma,
                'DELTAR2':-1., 'KPIV':-1. }
    # Primordial Power Spectrum
    if (pk_feature != None):PowerSpec['PK_FEATURE'] = pk_feature

    #  Power spectrum normalisation
    if not(sigma8 is None) and  (deltar2 != None):
        print "%RD_COSMO: only sigma8 or deltar2 should be provided"
        print 'SIGMA8', sigma8
        return None

    if not(deltar2 is None) :
        if (verbose):print "deltar2"
        PowerSpec['DELTAR2'] = deltar2
        PowerSpec['NORM'] = 0
        PowerSpec['SIGMA8'] = -1
        PowerSpec['KPIV'] =  0.002
        if (kpiv != None): PowerSpec['KPIV']

    # Power Spectrum alteration
    if (gamma is None):
        if (verbose): print 'gamma missing'
        # ; spectrum shape factor Gamma^*see refs in Ma 1998, ApJ 508, L5
        gamma = geom['OMEGA_M']*geom['H']
        if (geom['OMEGA_B'] !=0 ):
            gamma *= math.exp(-geom['OMEGA_B']*(1.+math.sqrt(2*geom['H'])/geom['OMEGA_M']))  
        PowerSpec['GAMMA'] = gamma


        #construct cosmo structure
        cosmo = dict(geom, **PowerSpec) # dic0.update(dic1)
 

    if not(age is None) :
        if (verbose):print "age"
        ageh0 = t_a(cosmo,1.e-8)  #age of the universe
        ageh0 = ageh0/cosmo['H0']*3.085678E19/3.15581498E7/1E9
        cosmo['AGE'] = ageh0

    if (verbose): print 'rd_cosmo c est fini'

    return cosmo
