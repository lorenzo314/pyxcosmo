###


import numpy as np
import time
import __builtin__
 
from icosmo_cosmo import mk_evol
from xgen_bound import xgen_bound
from icosmo_cosmo import xgen
from compute_step import compute_step

from get_scalar_argument import get_int_argument
from get_scalar_argument import get_scalar_argument

__version__ = 5

# LF 08/12/15: look at combine_dndzdldt.pro in combine_dndzdldt.pro
# RG: init_cmp_dndcr does the computation for a few redshifts; this one interpolates between missing redshifts 
def combine_dndzdldt (data_list, cosmo, model_struct, fluxmes_struct, gtype='ODE', verbose=False):
    """
        Resample dndzdldt, flux_ration, rcore,  on a new grid of redshift, given in the data_list.
        If needed, interpolate data between two different redshifts.


        INPUTS
        grid position of each sample

        OUTPUT

        dndzdldt:
       dndzdldt['RCORETAB']    (150, 64)  (n_t, n_Z)
       dndzdldt['DIMTAB']      (180, 150, 64) (n_l, n_t, n_z) 
         only  dndzdldt['DIMTAB'][:,0,:] (n_l, n_z)  usefull, dim means attenuated
       dndzdldt['DATAB']  (150, 64)  (n_t, n_z)
      
       dndzdldt['TTAB_TIMES_DNDZDLDT']  (180, 150, 64) (n_l, n_t, n_z) 
       dndzdldt['DLINZ']    scalar, float, logarithmic increment of Z, 0.0277777
       dndzdldt['DLOGTEMP'] scalar, float, logarithmic increment of Temperature, 0.0309
       dndzdldt['DLOGL']    scalar, float, logarithmic increment of Luminosity,  0.1145
       dndzdldt['LOGDIMTAB_SUR_DLOGCR']       (180, 150, 64) (n_l, n_t, n_z)
       dndzdldt['LOGFLUXRATIOTAB_SUR_DLOGCR']  (180, 150, 64) (n_l, n_t, n_z)
       dndzdldt['VEC_HR_BOUNDS']    vector n_z+1
       dndzdldt['BINDCR_TIMES_BINDHR']   (64, 64)  (n_z, n_z)
        

        EXAMPLE
        index = compute_virtual_index(grid, values)

        HISTORY:
        RG 09 October 2015 version 5 : add gtype, see init_cmp_dndcr, mk_evol
        RG 8 May add bind_cr_times_bind_hr in the output
        RG 7 May the get_scalar functions are in a different file version 4
        RG 4 May 2015 input can be scalar or arry of 1 version 3
        R Gastaud fecit,  1 December 2014

    """
    start_time = time.time()

    # zvector = all the red shifts
    zvector = data_list.pop(0)  # 64
    commb = data_list.pop(0)  # 33  # selection fo z vector
    nn = len(data_list)

    in_z = np.zeros(nn)
    for i in range(0, nn): 
        in_z[i] = data_list[i]['Z']

    # check consistency of data list
    if max(abs(np.sort(in_z)-in_z)) > 0. : 
        print 'combine_dndzdld in_z not sorted'
        # raise an exception ??
        #return -1

    n_z = get_int_argument(model_struct['N_Z']) # 64
    n_l = get_int_argument(model_struct['N_L']) #180
    n_t = get_int_argument(model_struct['N_T']) #150
    #(180, 150, 64)   = (n_l, n_t, n_z)
    if (verbose): print 'n_l=',n_l, 'n_t=', n_t, 'n_z=', n_z

    if np.isscalar(fluxmes_struct['APERTURE']):
        n_ap = 1
    else:
        if np.atleast_1d(fluxmes_struct['APERTURE']):
            n_ap = fluxmes_struct['APERTURE'].size
        else:
            n_ap = len(fluxmes_struct['APERTURE'])  # list ??

    full_dndzdldt = np.zeros([n_l, n_t, n_z])

    full_fluxratio =  np.zeros([n_l, n_t, n_ap, n_z])
    #  (180, 150, 1, 64)

    full_rcore = np.zeros([n_t, n_z])
    # full_rcore.shape   (150, 64)

    red = 0
    for z in zvector:
       index = np.argmin(abs(in_z-z))
       diff = z-in_z[index]
       if (abs(diff) < 1e-6):
           if (verbose):print 'no interpolation for ', red, index
           data = data_list[index]
           full_dndzdldt[:,:,red] = data['DNDZDLDT']
           # n_ap = 1 ==> 0 hardcoded beware !!!!
           full_fluxratio[:,:,0,red] = data['FLUXRATIOTAB']
           full_rcore[:,red] = data['RCORETAB']
           # check that the temperature and luminosity vectors are the same
           #data['TEMP']
           #data['L']
       else:
           # interpolation
           #print 'interpolation for ', red
           if (diff > 0):
               redinf=index
               redsup = index+1
           else:
               redinf=index-1
               redsup = index
           t = (z-in_z[redinf])/(in_z[redsup]-in_z[redinf])
           #
           data1 = data_list[redinf]
           data2 = data_list[redsup]
           if (verbose):print 'interpolation red, i, t', red, i, t, redinf, redsup
           if (verbose):print ' '

           full_dndzdldt[:,:,red] = (1-t)*data1['DNDZDLDT'] + t*data2['DNDZDLDT']
            # n_ap = 1 ==> 0 hardcoded beware !!!!
           full_fluxratio[:,:,0,red] = (1-t)*data1['FLUXRATIOTAB'] + t*data2['FLUXRATIOTAB']
           full_rcore[:,red] = (1-t)*data1['RCORETAB'] + t*data2['RCORETAB']
       red=red+1
       #


    ###########  mk_ evol
    evol = mk_evol(cosmo, z=zvector, gtype=gtype)
    #evol = mrdfits('input_data/evol.fits')

    # evol['DL'].shape # (64,)  nz
    # evol['DA'].shape # (64,)  nz

    Dimtab = data['L'].reshape(n_l,1)/1.E44/(np.square(evol['DL'])).reshape(1,n_z)
    # Dimtab.shape  (180, 64)   (n_l, n_z)

    full_Ttab_times_dndzdldt   = data['TEMP'].reshape(1,n_t,1)*full_dndzdldt

    dlog_temperature = compute_step(data['TEMP'], n_t, log_flag=True)
    dlog_luminosity = compute_step(data['L'], n_l, log_flag=True)
    dlin_redshift = compute_step(zvector, n_z)

    n_cr = get_int_argument(model_struct['N_CR'])
    n_hr = get_int_argument(model_struct['N_CR_2'])
   
    min_cr = get_scalar_argument(model_struct['CR_MIN'])
    max_cr = get_scalar_argument(model_struct['CR_MAX'])
    min_hr = get_scalar_argument(model_struct['CR_MIN_2'])
    max_hr = get_scalar_argument(model_struct['CR_MAX_2'])

    vec_cr = xgen(min_cr, max_cr, n_cr, True)
    dlog_cr =  compute_step(vec_cr, n_cr, log_flag=True)

    vec_hr = xgen(min_hr, max_hr, n_hr, True)
    dlog_hr =  compute_step(vec_hr, n_hr, log_flag=True)

    vec_cr_bounds = xgen_bound(min_cr, max_cr, n_cr, log_flag=True)
    vec_hr_bounds = xgen_bound(min_hr, max_hr, n_hr, log_flag=True)

    bind_cr = vec_cr * dlog_cr
    # bindcr = rebin(reform(bindcr,nCR,1,1),nCR,nHR,n_ap)
    bind_hr = vec_hr * dlog_hr  
    #bindhr = rebin(reform(bindhr,1,nHR,1),nCR,nHR,n_ap)
    bind_cr_times_bind_hr = bind_cr * bind_hr

    if (fluxmes_struct['OBSTYPE'] == 'crhr'):
        #full_fluxratio.shape (180, 150, 1, 64)
        log_flux_ratiotab_sur_dlogcr = np.log(full_fluxratio)/dlog_cr
        log_dim_sur_dlogCR = np.log(Dimtab)/dlog_cr


    dndzdldt = { 'L':data['L'], 'Z':zvector, 'TEMP': data['TEMP'],
                 'NL':n_l,'NZ':n_z,'NT':n_t,
                 'DLOGL':dlog_luminosity, 'DLINZ':dlin_redshift, 'DLOGTEMP':dlog_temperature,
                 'DNDZDLDT':full_dndzdldt, 'TTAB_TIMES_DNDZDLDT':full_Ttab_times_dndzdldt,
                 'LOGFLUXRATIOTAB_SUR_DLOGCR':log_flux_ratiotab_sur_dlogcr, 'LOGDIMTAB_SUR_DLOGCR':log_dim_sur_dlogCR,
                 'VEC_CR':vec_cr, 'VEC_HR':vec_hr, 'VEC_CR_BOUNDS':vec_cr_bounds, 'VEC_HR_BOUNDS':vec_hr_bounds,
                 'DLOG_CR':dlog_cr, 'DLOG_HR':dlog_hr,
                 'RCORETAB':full_rcore, 'DA':evol['DA'], 'DIMTAB':Dimtab,
                 'BINDCR_TIMES_BINDHR':bind_cr_times_bind_hr}

    # FLUXRATIOTAB is not used  by dndzdldt_to_dndcrdhr
    #  DATAB   (150,64) = (n_t, n_z)  is replaced by DA     (64) = n_z
    #  TTAB2D  (150,64) = (n_t, n_z)   is replaced by TEMP (150) =n_t
    #  ZTAB2D  (150,64) = (n_t, n_z)   is replaced by Z     (64) = n_z
    #
    elapsed_time = time.time() - start_time
    if (verbose): print 'combine_dndzdldt elapsed time', elapsed_time

    return dndzdldt
