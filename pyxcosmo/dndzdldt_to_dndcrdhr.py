###  copy of dndzdldt_to_dndcrdhr_3.py in /Users/gastaud/icosmo/codePython/brouillon13/delivery
###  correct selftype in seltype rg Lun 26 jan 2015 12:20
###  replace fast_histogram with call_histogram_integer  rg  3 May 2015 version 2
# 7 May 2015 : model_struct['PHYS_RC'] may be a scalar, version 3
# RG 11 May 2015 version 4 get nx, ny

import numpy as np
import __builtin__
import math
import time
import scipy.ndimage as ndimage

from get_scalar_argument import get_int_argument
from get_scalar_argument import get_scalar_argument

__version__ = 4

from call_histogram_integer import call_histogram_integer

def compute_virtual_index(grid, values, flag=True):
    """
        Compute the virtual indexes.

        INPUTS
        grid position of each sample

        OUTPUTS
        indexes 

        EXAMPLE
        index = compute_virtual_index(grid, values)

	HISTORY:
        RG 11 May 2015 version 4 get nx, ny
          count rate idl 400 and python only 40, the only relevant ones!
        R Gastaud fecit,  24 oct 2014

    """
    x = grid
    nx = x.size
    steps = x[1:]-x[0:-1]
    diff = steps.max()-steps.min()
    if (diff > 5e-7):
        print 'error in compute_virtual_index'
        return -1

    step=steps.mean()
    origin = x[0]
      
    virtual_i = (values-origin)/step
    virtual_i =  virtual_i.clip(min=0, max=nx-1)

    return  virtual_i


def compute_probdet(proba, x, y, xnew, ynew):

    """
   Interpolate detection probability 
      func_count_rate nx
      func_rcore  ny 

       count_rate ( nL, nT, nZ)
       rcore (nT, nZ)
    """
    #start_time = time.time()

    virtual_iy = compute_virtual_index(y, ynew)
    virtual_ix = compute_virtual_index(x, xnew)

    myshape = virtual_ix.shape # nL, nT, nZ
    nL=(virtual_ix.shape)[0]  # nL=180

    nn = virtual_iy.size # # nT, nZ
    virtual_iy = virtual_iy.reshape(nn)
    # nn is nT*nZ = 150*64

    virtual_ixx = virtual_ix[0,:,:]
    virtual_ixx = virtual_ixx.reshape(nn)
    coords = np.array([virtual_iy, virtual_ixx])    
    probdet = np.zeros((nL, nn))

    for i in range(0, nL):
        virtual_ixx = virtual_ix[i,:,:]
        virtual_ixx = virtual_ixx.reshape(nn)
        #coords = np.array([virtual_iy, virtual_ixx])
        # do not create coords in a loop =>
        #              avoid memory problem, coords is an object
        coords[1,:] = virtual_ixx
        probdet[i,:] = ndimage.map_coordinates(proba, coords, prefilter=True, order=1)
        # , output=probdet[i,:,:]

    probdet = probdet.reshape(myshape)
    #elapsed_time = time.time() - start_time
    #print 'interpolation elapsed time', elapsed_time

    return probdet

########################

#from fast_histogram import fast_histogram
# 
def compute_conversion(probdet, dndzdldt, model_struct, aperture, index_mes_cr, mes_hr ):

    #start_time = time.time()

    n_apertures = aperture.size
 
    integrand = dndzdldt['TTAB_TIMES_DNDZDLDT']*probdet

    d_linz = get_scalar_argument(dndzdldt['DLINZ'])
    d_logl = get_scalar_argument(dndzdldt['DLOGL'])
    d_logtemp = get_scalar_argument(dndzdldt['DLOGTEMP'])

    in_voxel = d_linz*d_logtemp*d_logl

    n_cr = get_int_argument(model_struct['N_CR'])
    n_z =  get_int_argument(model_struct['N_Z'])
    n_l =  get_int_argument(model_struct['N_L'])
    n_m =  get_int_argument(model_struct['N_M'])
    n_t =  get_int_argument(model_struct['N_T'])
    n_hr = get_int_argument(model_struct['N_CR_2'])

    #print 'n_l=', n_l, ' n_t=', n_t, ' n_z=', n_z
    ### compute binlocate1  CR ##
    # index_mes_cr comes from FLUXMES_STRUCT.indexcrtab.CR_PERL_TIMESD2
    #print 'before logfluxratio_sur_dlogcr.shape=', dndzdldt['LOGFLUXRATIOTAB_SUR_DLOGCR'].shape
    #print ' logfluxratio_sur_dlogcr=', dndzdldt['LOGDIMTAB_SUR_DLOGCR'].shape, 'index_mes_cr=', index_mes_cr.shape

    logfluxratio_sur_dlogcr = dndzdldt['LOGFLUXRATIOTAB_SUR_DLOGCR']
    if (logfluxratio_sur_dlogcr.ndim == 4):
        logfluxratio_sur_dlogcr = logfluxratio_sur_dlogcr[:,:,0,:]
        logfluxratio_sur_dlogcr = logfluxratio_sur_dlogcr.reshape(n_l, n_t, n_z)

    logdim_sur_dlogcr = dndzdldt['LOGDIMTAB_SUR_DLOGCR']
    if (logdim_sur_dlogcr.ndim == 2):
        logdim_sur_dlogcr = logdim_sur_dlogcr.reshape(n_l, 1, n_z)

    index_mes_cr = index_mes_cr.reshape(1, n_t, n_z)
    #print 'after logfluxratio_sur_dlogcr.shape=', logfluxratio_sur_dlogcr.shape
    #print ' logfluxratio_sur_dlogcr=', logfluxratio_sur_dlogcr.shape, 'index_mes_cr=', index_mes_cr.shape

    binlocate1 = (np.floor(index_mes_cr + logdim_sur_dlogcr + logfluxratio_sur_dlogcr)).astype(np.int64)
    #print ' binlocate1.shape=', binlocate1.shape

    binlocate1 = binlocate1.clip(min=-1, max=n_cr)
    binlocate1.min(), binlocate1.max()
    binlocate1 += 1


    ### compute binlocate2  HR ##
    vec_hr_bounds = dndzdldt['VEC_HR_BOUNDS']
    ratio = vec_hr_bounds[1:]/vec_hr_bounds[0:-1]
    step_hr = math.log(np.median(ratio)) 
    origin_hr = math.log(min(vec_hr_bounds))

    #  mes_hr comes from FLUXMES_STRUCT_CRTAB_2
    log_mes_hr = np.log(mes_hr)
    binlocate2 = (np.floor((log_mes_hr-origin_hr)/step_hr)).astype(np.int64)

    binlocate2 = binlocate2.clip(min=-1, max=n_hr)
    binlocate2 +=1

    ### call histogram_integer  
    sizex = binlocate1.max()+1
    #print 'sizex', sizex, 'n_cr+2', (n_cr+2)

    sizey = binlocate2.max()+1
    #print 'sizey', sizey, 'n_hr+2', (n_hr+2)


    #print 'shape before', binlocate1.shape, binlocate2.shape
    #print 'dtype ', binlocate1.dtype, binlocate2.dtype
    nn = binlocate1.size
    binlocate1 = binlocate1.reshape(nn) 
    binlocate2_3d = np.tile(binlocate2, (n_l, 1, 1))
    binlocate2_3d = binlocate2_3d.reshape(nn)
    #print 'shape after', binlocate1.shape, binlocate2_3d.shape

    ri = call_histogram_integer(binlocate1, binlocate2_3d,  (n_cr+2),  (n_hr+2), nn, location=model_struct['LIBDIR'])
    # checking ri, it worded...
    #from fitsutils import mrdfits
    #ref_ri = mrdfits('../fast_histogram/ri.fits',0)
    #ref_ri.dtype, ref_ri.size, ri.dtype, ri.size
    #diff = ri-ref_ri
    #print ' ri difference ', diff.min(), diff.max()

    integrand = integrand.reshape(integrand.size)

    ###  counting ...  ##
    ### 
    dndcrdhr = np.zeros([n_apertures, n_hr,n_cr])
    #dndCRdHR = fltarr(nCR,nHR,n_apertures)

    for j in range(0, n_hr):
        #print j
        for i in range(0, n_cr):
            #  n_cr+2 = sizex
            ij = (i+1) + (n_cr+2)*(j+1)
            if(ri[ij+1] > ri[ij]):
                #print i,j
                # beware last is not included
                select = integrand[ri[ri[ij]:ri[ij+1]]]
                dndcrdhr[0,j,i] = select.sum()
    #
    # normalise
    dndcrdhr *= in_voxel
    #  (64, 64, 1)   and (64, 64)  does not broadcast???
    #for i in range(0, n_apertures):
        #dndcrdhr[:,:,i] /= dndzdldt['BINDCR_TIMES_BINDHR']

    dndcrdhr /= dndzdldt['BINDCR_TIMES_BINDHR']
    #elapsed_time = time.time() - start_time
    #print 'conversion elapsed time', elapsed_time

    return dndcrdhr

######################################

def  dndzdldt_to_dndcrdhr(dndzdldt, cosmo_h, seltype, crtab, sfunctab, model_struct, aperture, index_mes_cr, mes_hr, verbose=False):

    """
  Convert the dndzdldt data to dndcrdhr data.
  This is done in two steps
    1) Interpolation of the pre-computed table cr, hdr
    2) Conversion 

  dndzdldt:
       dndzdldt['RCORETAB']    (150, 64)  (n_t, n_Z)
       dndzdldt['DIMTAB']      (180, 150, 64) (n_l, n_t, n_z) 
         only  dndzdldt['DIMTAB'][:,0,:] (n_l, n_z)  usefull, dim means attenuated
       dndzdldt['DATAB']  (150, 64)  (n_t, n_z)
      in conversion:
       dndzdldt['TTAB_TIMES_DNDZDLDT']  (180, 150, 64) (n_l, n_t, n_z) 
       dndzdldt['DLINZ']    scalar, float, logarithmic increment of Z, 0.0277777
       dndzdldt['DLOGTEMP'] scalar, float, logarithmic increment of Temperature, 0.0309
       dndzdldt['DLOGL']    scalar, float, logarithmic increment of Luminosity,  0.1145
       dndzdldt['LOGDIMTAB_SUR_DLOGCR']       (180, 150, 64) (n_l, n_t, n_z)
       dndzdldt['LOGFLUXRATIOTAB_SUR_DLOGCR']  (180, 150, 64) (n_l, n_t, n_z)
       dndzdldt['VEC_HR_BOUNDS']    vector n_z+1
       dndzdldt['BINDCR_TIMES_BINDHR']   (64, 64)  (n_z, n_z)


    cosmo_h: scalar, float

    seltype: one string selection function ?

    crtab:
       crtab['CR_PERL_TIMESD2']  count rate (150, 64)  (n_t, n_z)
 
    sfunctab
       sfunctab['RCORE']    (10,)   ny=10
       sfunctab['COUNTR']   (400,) = only usefull [0:40], nx=40
       sfunctab['PROBA']    (400,)=  (10, 40)  (ny, nx)
     
    model_struct
        n_z =  int(model_struct['N_Z'][0])  64
        n_l =  int(model_struct['N_L'][0]) 180
        n_t =  int(model_struct['N_T'][0]) 150
        (180, 150, 64)   = (n_l, n_t, n_z)
        n_m =  int(model_struct['N_M'][0])   200
        n_cr = int(model_struct['N_CR'][0])   64
        n_hr = int(model_struct['N_CR_2'][0]) 64
        location=model_struct['LIBDIR']
        model_struct['PHYS_RC'][0] 

    aperture: vector of float, only the first value, aperture[0] , is used, and the size of the vector

     index_mes_cr: table of float  (150, 64) (n_t, n_z)
     mes_hr: table of float        (150, 64) (n_t, n_z)

XAMIN_struct
	XAMIN_struct.seltype  string
	XAMIN_struct.crtab    structure
	XAMIN_struct.sfunctab structure

	HISTORY:
	7 May 2015 : model_struct['PHYS_RC'] may be a scalar, version 3
 
    """

    start_time = time.time()

    # model_struct['PHYS_RC'][0]
    #  is an array of dimension one, value -1
    #   or a scalar 7 May 2015 version 3
    phys_rcore_0 = get_scalar_argument(model_struct['PHYS_RC'])
    if (phys_rcore_0 > 0):
        phys_rcore = model_struct['PHYS_RC']
        physical_rcore = phys_rc / cosmo_h
    else:
        physical_rcore = dndzdldt['RCORETAB']


    if (seltype != 'selfunc'): 
        print 'no interpolation'
        print ' je ne sais pas quoi faire'
        raise Exception(" seltype is not selfunc in dndzdldt_to_dndcrdh")

    #pass

    count_rate = crtab['CR_PERL_TIMESD2']
    # (150, 64)  (nT, nZ)
    (n_t, n_z) = count_rate.shape
    count_rate  = count_rate.reshape(1,n_t, n_z)  #7_May
    #  remark dim means attenuation
    #  does not depend upon Temperature
    # dndzdldt['DIMTAB'].shape
    # (180, 150, 64)   = (nL, nT, nZ)
    #  we can use dndzdldt['DIMTAB'][:,0,:]
    #  (nL, nZ)
    dim_tab = dndzdldt['DIMTAB']
    if (dim_tab.ndim == 2):
        (n_l, n_z) = dim_tab.shape
        #print ' dim tab ', n_l, n_z
        dim_tab = dim_tab.reshape(n_l, 1, n_z)

    #print ' count_rate.shape before=', count_rate.shape, ' dim_tab ', dim_tab.shape
    count_rate = count_rate*dim_tab
    #print ' count_rate.shape after=', count_rate.shape
    # (180, 150, 64)   = (nL, nT, nZ)

    n_l = get_int_argument(model_struct['N_L'])

    if 'DATAB' in dndzdldt:
        if(verbose): print "use DATAB"
        da_tab = dndzdldt['DATAB']
    else:
        if(verbose): print "use DA"
        da_tab = dndzdldt['DA']

    #print 'da_tab.shape', da_tab.shape, 'dndzdldt[RCORETAB]',dndzdldt['RCORETAB'].shape
    rcore = dndzdldt['RCORETAB']/da_tab*(180*3.6/math.pi)
    # (150, 64)  (nT, nZ)

    ## we want to interpolate on count_rate and rcore 

    func_rcore = sfunctab['RCORE']
    # 10
    ny = func_rcore.size

    func_proba = sfunctab['PROBA']
    #  1 dimension, nx*ny
    nn = func_proba.size
    nx = nn/ny 
    func_proba = func_proba.reshape(ny,nx)
    #  proba in python (10, 40)  ny, nx


    func_count_rate = sfunctab['COUNTR']
    #  1 dimension, nx*ny  in IDL, but in python only nx
    # nx = func_count_rate.size/ny this is true for idl only
    func_count_rate = func_count_rate[0:nx]  # this always works 

    #print 'func shapes', func_count_rate.shape, func_rcore.shape, func_proba.shape
    # ((40,), (10,), (10, 40))


    #probdet =interpol_2dfast(sfunctab.proba,sfunctab[0].countr,sfunctab.rcore,cr,rcore)

    # (180, 150, 64)
    #print 'output shape', count_rate.shape, rcore.shape

########################################
    probdet = compute_probdet(func_proba, func_count_rate, func_rcore, count_rate, rcore)


    elapsed_time1 = time.time() - start_time
    #print 'dndzdldt interpolation time ', elapsed_time1

    dndcrdhr = compute_conversion(probdet, dndzdldt, model_struct, aperture, index_mes_cr, mes_hr )


    elapsed_time2 = time.time() - start_time
    if (verbose): print 'total elapsed time for dndzdldt_to_dndcrdhr', elapsed_time2

    #print 'conversion elapsed time ', elapsed_time2-elapsed_time1


    return dndcrdhr
