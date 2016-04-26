#!/usr/bin/env python
# 24 february 2015 R Gastaud  brouillon 26

import numpy as np

from icosmo_cosmo import mk_evol
from icosmo_cosmo import mk_nm # brouillon 5

from check_param import check_param # brouillon24
from cmp_dndzdldt import  cmp_dndzdldt # brouillon15

from fitsutils import mwrfits
import pdb

__version__ = 3

# LF 08/12/15: look at init_cmp_dndcr.pro in /dsm/sappccosmos0/work/lfacciol/faites_par_nicolas
# RG: output is dn / (dl dt dz)
# init_cmp_dndcr (given by Nicolas) seems misleading, maybe means it initiates the computation of dncdr
def init_cmp_dndcr(cosmo, MT_struct, LT_struct, FLUXMES_struct, MODEL_struct, verbose=False, directory=None, gtype='ODE'):
    """
        Initialise dndcrdhr from dndzdldt, mainly calling
        cmp_dndzdldt

        Inputs
        ----------
        cosmo:
        fluxmes_struct:
        model_struct:


        Returns
        -------
        list of dictionaries

        HISTORY:
         RG 09 October 2015 version 3 : add gtype  , see mk_nm, mk_evol, cmp_dndzdldt, combine_dndzdldt
         RG 04 May 2015 version 2 : remove unused argument clean
    """

    # Check structures and initialize fundamental parameters if necessary
    
    save = False
    if not (directory == None): save=True
    # diri='intermediates/'
    result = []
    #print 'lt_struc keys xxx ', LT_struct.keys(),  LT_struct['STRUCT_NAME']

    MT_struct = check_param(MT_struct, cosmo, 'mt', verbose=verbose)
    LT_struct = check_param(LT_struct, cosmo, 'lt', verbose=verbose)
    MODEL_struct = check_param(MODEL_struct, cosmo, 'model', verbose=verbose)
    FLUXMES_struct = check_param(FLUXMES_struct, cosmo, 'fluxmes', verbose=verbose)
    print 'init_cmp_dndcr xxx FLUXMES_struct.keys()', FLUXMES_struct.keys()

    # /1/ Set up mass function
    z_ran = MODEL_struct['Z_RAN']
    n_z = MODEL_struct['N_Z']
    m_ran = MODEL_struct['M_RAN']
    n_m = MODEL_struct['N_M']
    if (verbose): print 'z_ran=',z_ran, 'n_z=',n_z, 'm_ran=', m_ran,'n_m=',n_m

    # mk_nm_array replaced by mk_evol+mk_nm
    evol = mk_evol(cosmo, ran_z=z_ran, n_z=n_z, gtype=gtype)

    # fit_tk and ctype are integrers, not logicals
    nmz = mk_nm(cosmo,evol, z=evol['Z'], m_ran=m_ran, n_m=n_m, profile=True, fit_nm=4, tinkerdelta=200, fit_tk=1, ctype=0, tinkerrhoc=True, gtype=gtype)

    result.append(evol['Z'])

    # now select only part of the data
    peigne =  MODEL_struct['Z_PEIGNE']
    index = np.arange(0, n_z, peigne)

    # add the last point if it was not selected
    if (index[-1] < (n_z-1)):
        index = np.hstack([index, n_z-1])
    result.append(index)

    # send_to_cmp_dndzdldt : loop on z
    # need a new dictionary for nmz[i]
    new_dict = dict()
    for i in index:
        for key in nmz.keys():    
            new_dict[key]=nmz[key][i]
        new_dict["M"]=nmz["M"]

        for key in new_dict.keys(): 
            if(new_dict[key].size == 1): 
                new_dict[key]=new_dict[key].reshape(1) 

        # M and Z are the independant variables, with # dimensions (e.g. M (60,) Z (64,)
        if (verbose): print 'i=', i, new_dict['Z'], nmz['Z'][i]
        #print 'init_cmp_dndcr fluxmes_struct keys', FLUXMES_struct.keys()
        #print 'init_cmp_dndcr fluxmes_struct aperture ', FLUXMES_struct['APERTURE']
        #pdb.set_trace()
        interm = cmp_dndzdldt(cosmo, new_dict, MT_struct, LT_struct, FLUXMES_struct, MODEL_struct, gtype=gtype, verbose=True)
        # def    cmp_dndzdldt( cosmo, nm,       mt_struct, lt_struct, fluxmes_struct, model_struct, verbose=False):
 
        #new_dict.keys()
        #print new_dict
        #for key in nmz.keys():print key, nmz[key].shape
        #for key in new_dict.keys(): print key, new_dict[key].shape
        if(save):
            print 'save for index', i
            suffix = "{0:02d}".format(i)
            mwrfits(new_dict, directory+'new_dict'+suffix+'.fits' )
            mwrfits(interm, directory+'interm'+suffix+'.fits' )
       
        result.append(interm)
    
    if(save):
        mwrfits(nmz, directory+'nmz.fits' )
        #mwrfits(cosmo, diri+'cosmo.fits' )
        #mwrfits(MT_struct, diri+'MT_struct.fits' )
        #mwrfits(LT_struct, diri+'LT_struct.fits' )
        #mwrfits(FLUXMES_struct, diri+'FLUXMES_struct.fits' )
        #mwrfits(MODEL_struct, diri+'MODE_struct.fits' )
 
    return result 

