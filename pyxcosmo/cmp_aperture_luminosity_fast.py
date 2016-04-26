#!/usr/bin/env python
#  13 March 2015  R Gastaud

import numpy as np
import matplotlib.pyplot as plt
import __builtin__
from scipy import interpolate
from math import log

from scipy.interpolate import interp1d
from icosmo_cosmo import extrap1d
from scipy import ndimage
from fitsutils.fitstables import mrdfits
import os.path
__version__= 3

def cmp_aperture_luminosity_fast(table, outaperture, r500angular, xc=.303, inlum=1., verbose=False, order=1):
    """
    Compute the luminosity function of outaperture, r500angular, xc 
    This is done by interpolating a table
    Interpolation on xc , r500angular
    Extrapolation on outaperture

    INPUT : 
         table:  file name of the input table, luminosity versus outaperture, r500angular, xc 
         outaperture: vector, out aperture radius in which the luminosity is computed, in r500angular unit
         r500angular: vector, size of the galaxy, value in arcsec of r500 on the detector
         xc: the "xc" parameter in the AB model (default 0.303 Piffaretti 2010). 
                Useful to parametrize evolution or non-self-similar relation.
         inlum: input luminosity
         order: order of the interpolation, default 1, linear
         verbose: default false
         
         
    OUTPUT: interpolated luminosity

    Remark: input_aperture is hardcoded to 5 r500angular
    Bibliography: 
    Arnaud, M., Pratt, G.W., Piffaretti, R., Boehringer, H., Croston, J.H., Pointecouteau, E., 2010, A & A 517, 92,
    The universal galaxy cluster pressure profile from a representative sample of nearby systems (REXCESS) 
    and the Y_SZ-M_500 relation.
    ------
    01 May   2015 RG  interpolation logarithmic of outaperture
    30 April 2015 add parameter order, default 1 and use verbose for print
    13 March 2015 raise an exception when file is not found
    20 jan 2015 created by Rene Gastaud.
    """

    # verifications
    n_out = outaperture.size
    n_xc = xc.size
    n_ang = r500angular.size

    if(n_out==1): n_out = max([n_xc,n_ang])

    if (verbose):
        if(n_xc !=  n_out): print 'cmp_aperture_luminosity_fast warning n_xc !=  n_out', n_xc, n_out
        if(n_out != n_ang): print 'cmp_aperture_luminosity_fast  warning  n_ang !=  n_out', n_ang, n_out
        if(n_xc != n_ang): print 'cmp_aperture_luminosity_fast  warning n_ang !=  n_out', n_xc, n_out

    if (verbose): print 'n_out=', n_out, ' order=', order
    if(n_xc == 1): xc = np.repeat(xc, n_out)

    # I/Read table
    fstatus = os.path.isfile(table)
    print "Table", table
    if (verbose): print "read cmp_aperture_luminosity_fast file"+table+" ", fstatus
    if (fstatus):
        tab = mrdfits(table)
    else:
        raise Exception("cmp_aperture_luminosity_fast file"+table+" not found")

    m_ang = tab['R500ANGULAR'].size
    m_out = tab['OUTAPERTURE'].size
    m_xc = tab['XC'].size

    index_r500angular = np.arange(m_ang)
    index_xc = np.arange(m_xc)
    index_outaperture = np.arange(m_out)

    interpol_r = interp1d(tab['R500ANGULAR'], index_r500angular)
    ## interpolation does not work
    #r_interpolated = interpol_r(r500angular)
    extrapol_r = extrap1d(interpol_r)
    r_extrapolated = extrapol_r(r500angular)
    if (r500angular.max() > tab['R500ANGULAR'].max()):
        print 'WARNING cmp_aperture_luminosity_fast extrapolation on r500angular ', r500angular.max(), '>', tab['R500ANGULAR'].max()


    interpol_xc = interp1d(tab['XC'], index_xc)
    xc_interpolated = interpol_xc(xc)
    extrapol_xc =  extrap1d(interpol_xc)
    xc_extrapolated = extrapol_xc(xc)

    #interpol_outaperture = interp1d(log(tab['OUTAPERTURE']), index_outaperture)
    #extrapol_outaperture = extrap1d(interpol_outaperture)
    #outaperture_extrapolated = extrapol_outaperture(log(outaperture))

    #interpol_outaperture = interp1d(tab['OUTAPERTURE'], index_outaperture)
    interpol_outaperture = interp1d(np.log(tab['OUTAPERTURE']), index_outaperture)
    extrapol_outaperture = extrap1d(interpol_outaperture)
    if (outaperture.max() > tab['OUTAPERTURE'].max()):
        print 'WARNING cmp_aperture_luminosity_fast extrapolation on outaperture ', outaperture.max(),'>', tab['OUTAPERTURE'].max()
    outaperture_extrapolated = extrapol_outaperture(np.log(outaperture))

    #z = mrdfits('z.fits',0)
    #y = mrdfits('y.fits',0)
    #x = mrdfits('r_interpolated.fits',0)
    #coords = (z,y,x)
    #print 'xc_extrapolated ', xc_extrapolated.size

    coords = (outaperture_extrapolated, xc_extrapolated, r_extrapolated)

    outlum = inlum*ndimage.map_coordinates(tab['APERTURE_LUMINOSITY'], coords, order=order)
 
    return outlum
