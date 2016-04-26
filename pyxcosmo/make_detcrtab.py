#---------------------- Dedicated Function: make_detcrtab
#  BEWARE python first axis is the second!
#  BEWARE crtab['NH'] must be sorted

#---------------------- -----------------------

from fitsutils import mrdfits
from compare_idlpy import compare_idlpy
import numpy as np
import matplotlib.pyplot as plt
import __builtin__
from scipy import interpolate
from math import log

__version__ = 3

def make_detcrtab(tabdir,middlefix,filtab=None,nh=4e20,Ttab2D=None,Ztab2D=None, verbose=False):
    """
    Create a table of probability of detection count rate function of temperature and red shift z, 
    for a given NH level, if the filters are given.
    Without filters, it creates a table the hrtab versus temperate and red shift

    Algorithm
    ---------
    1) It reads an fits file containing the array of cr for which the table has already been computed,  
    with the corresponding 'middlefix'.
    2) Check that the given NH is in the array of precomputed NH.
    3) Select the two tables of CR for the NH below and above the given NH.
    4) Interpolate between these two tables.
    5) If asked, interpolate for the new array of temperature and new array of redshift.
    Beware, transposition needed in 5!
  
    Arguments
    ---------
    tabdir: string, the directory of the fits file of the cube.
    middlefix: string, the middle fix of the fits file of the cube.
 
    Keyword argument
    ---------
    filtab: string array, the name of the filters, if not set create hrtab
    nh: scalar, number : the value of the NH.
    Ttab2D: 2D array (n1,n2) , replicate Ttab2D[:,0]
    Ztab2D: 2D array (n1,n2) , replicate Ztab2D[0,:]
    verbose: bool, flag, if set print some messages.

    Returns
    -------
    A dictionary of arrays, with:
         the count rate (array 2d, n1, n2)  or  (array 2d, n2, n1)
         the temperature (array 1d, n1)
         the red shift (array 1d, n2)
         the NH (scalar)


    Update:
    ------
    11 May 2015 version 3 if filtab=None, it is make_meshr
    9 May 2015 version 2 temperature and z can be one dimensional (vector)
    3 oct 2014 created by Rene Gastaud.

    """
    if (verbose):
        if (filtab==None):
            print ' function make_detcrab use as make_meshr'
        else:
            print ' function make_detcrab'
            
    if(filtab != None):
        dict = {'Thin1': 1, 'Medium': 2, 'Thick': 3}
        m1f = dict[filtab[0]]
        m2f = dict[filtab[1]]
        pnf = dict[filtab[2]]
        filt_id = (m1f + 3*m2f + 9*pnf)
        filt_id = '0'+str(filt_id)
        ####  beware to be replaced only filter thin1 implemented now
        filt_id = '00'
        filename = tabdir+'/CRtabs_courant/alldet_b'+middlefix+'_FiltID'+filt_id+'_crtab.fits'
    else:
        #  make meshr and not detcrtab ??
        filename = tabdir+'/CRtabs_courant/alldet_'+middlefix+'_crtab.fits'
        # make_meshrtab

    if (verbose): print 'read ', filename

    crtab = mrdfits(filename,1)

    if (nh > crtab['NH'].max()):
        print "warning nh too big"
        ii = crtab['NH'].argmax()
        tab =  crtab['CR_PERL_TIMESD2'][ii,:, :]

    elif (nh < crtab['NH'].min()):
        print "warning nh too small"
        ii = crtab['NH'].argmin()
        tab =  crtab['CR_PERL_TIMESD2'][ii,:, :]

    else:
        i1 = np.where(crtab['NH'] <=  nh)[0].max()
        i2 = np.where(crtab['NH'] >=  nh)[0].min()
        if (verbose): print 'i1', i1, 'i2', i2
        tab1 = crtab['CR_PERL_TIMESD2'][i1,:, :]
        tab2 = crtab['CR_PERL_TIMESD2'][i2,:, :]
        if (i1 != i2):
           dx = log(crtab['NH'][i2]) - log(crtab['NH'][i1])
           x1  = (log(crtab['NH'][i2])-log(nh))/dx
           x2 = (log(nh)-log(crtab['NH'][i1]))/dx
           if (verbose): print 'x1', x1, 'x2', x2
        else:
            x1 = 1
            x2 = 0    
        tab = x1*tab1 + x2*tab2

    ### warning 
    ##   
    ## crtab['CR_PERL_TIMESD2']   7, 300, 39
    ## crtab['T']  39, crtab['Z']  300

    ##  reference
    ## reference['CR_PERL_TIMESD2']   150, 64
    ## reference['T']  150, referenc['Z']  64

    ###  TRANSPOSiTION
    tab_transposed = tab.transpose()

    if (Ztab2D != None) and  (Ttab2D != None):
        #  beware interpolate does not check the consistency at the creation
        #  interpolator = interpolate.interp2d(x, y, z)
        #  x, y with x the second axis, and y the first axis!
        if (verbose): print 'interpolate the table'
        interpolator = interpolate.interp2d(crtab['Z'], crtab['T'],  tab_transposed)
        if (Ztab2D.ndim == 2):
            z =  Ztab2D[0,:] #  (64,)
        else:
            z =  Ztab2D

        if(Ttab2D.ndim ==2):
            temperature = Ttab2D[:,0]  # ((150,)
        else:
            temperature = Ttab2D

        # attention python invert the axis
        cr_perL_timesD2 = interpolator(z, temperature)
    else:
        cr_perL_timesD2 = tab 
        temperature = crtab['T']
        z = crtab['Z']

    detcrtab = {'NH':nh,'T':temperature, 'Z':z, 'CR_PERL_TIMESD2':cr_perL_timesD2}

    return detcrtab
