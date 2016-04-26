import numpy as np
import math
import pdb

from fitsutils import mrdfits
from fitsutils import mwrfits

def calc_likelihood_dndmdz(observable, model, background=1e-12, verbose=False, n_m=None):
    
    """
    PURPOSE:
    -------------
    compute the likelihood
    
    INPUTS:
    -------------
    (observable) : a dictionnary containing:
    * the density dn/(dmdOmegadz), a grid of dimension n_m*n_z
    * a vector of m  , dimension n_m
    * a vector of z  , dimension n_z
    
    (model)
    
    
    (background): a scalar
    
    
    OUTPUT:
    -------------
    Z: likelihood , a scalar
    
    HISTORY:
    -------------
    Rene Gastaud, 13 January 2016
    """
    # Convolve the model with the measurement error

    n_z = model['Z'].size
    if n_m is None: n_m = model['M'].size
    #n_m = n_m/n_z
    mass = model['M'][0:n_m]
    if (verbose): print n_z, n_m

    # LF 29/01/16
    area_sq_deg = 100.0
    # WRONG: TO CONVERT 100 SQ DEGREES TO STERADIANTS WE MUST DO
    area = area_sq_deg / (180.0 / math.pi)**2
    # OR EQUIVALENTLY
    # area = math.radians(10)**2 
    # area=math.radians(100)**2 GIVES AN AREA 100 TIMES LARGER

    dlogm = (np.log(mass.max())-np.log(mass.min()))/(n_m-1.)
    dlogz = (np.log(model['Z'].max())-np.log(model['Z'].min()))/(n_z-1.)

    tab_mass_z = mass.reshape(1, n_m)*model['Z'].reshape(n_z,1)

    model_dndmdomdz = model['DNDMDOMDZ'].reshape([n_z, n_m])
    model_zm = model_dndmdomdz*dlogm*dlogz * area * tab_mass_z

    obs_dndmdomdz = observable['DNDMDOMDZ'].reshape([n_z, n_m])
    
    mass2 = observable['M'][0:n_m]

    observable_zm = obs_dndmdomdz*dlogm*dlogz*area*tab_mass_z

    print "xxxx"
    print model_zm.shape
    print observable_zm.shape
    print "xxxx"

    res = (model_zm + background -  (observable_zm + background)*np.log(model_zm + background)).sum()

    return -res
