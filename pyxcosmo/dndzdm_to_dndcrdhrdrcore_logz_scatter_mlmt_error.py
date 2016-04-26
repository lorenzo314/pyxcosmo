# This is version 10 of dndzdm_to_dndcrdhrdrcore_logz_scatter_mlmt
# 26/02/2016

import numpy as np

from scipy import interpolate
from scipy.interpolate import interp2d

import math

from pyxcosmo.icosmo_cosmo import convert_nm
from pyxcosmo.interpol2d import interpol2d

from pyxcosmo.extract_dict import extract_dict

from mass_to_luminosity import mass_to_luminosity
from mass_to_temperature import mass_to_temperature

import pdb
from fitsutils import mwrfits
import time

from pyxcosmo.check_param import check_param
from pyxcosmo.mk_mt import mk_mt
from pyxcosmo.mk_lt import mk_lt
from compute_virtual_index import compute_virtual_index
from scipy.ndimage.interpolation import map_coordinates

from convol_3d_fix_error import convol_3d_fix_error

def compute_selfun(sfunctab, this_rc, this_cr):
    func_count_rate = sfunctab['COUNTR']
    func_rcore = sfunctab['RCORE']
    func_proba = sfunctab['PROBA']

    this_cr = np.array(this_cr)
    this_rc = np.array(this_rc)

    ncr = this_cr.size
    nrc = this_rc.size

    this_cr = this_cr.reshape(ncr)
    this_rc = this_rc.reshape(nrc)

    if ncr == 1 and nrc > 1: this_cr = np.tile(this_cr[0], nrc)

    flag = True

    virtual_iy = compute_virtual_index(np.log(func_rcore), np.log(this_rc), flag=flag)
    if virtual_iy.size == 1: virtual_iy = virtual_iy[0]

    virtual_ix = compute_virtual_index(np.log(func_count_rate), np.log(this_cr), flag=flag)
    if virtual_ix.size == 1: virtual_ix = virtual_ix[0]

    coords = np.array([virtual_iy, virtual_ix])
    coords2 = coords.reshape(2, nrc)

    probdet = map_coordinates(func_proba, coords2, prefilter=True, order=1)

    return probdet

def our_remove(index, x1, x2, x3, x4):
     x1 = x1[index]
     x2 = x2[index]
     x3 = x3[index]
     x4 = x4[index]
     return x1, x2, x3, x4

"""
    Compute count rate and hardness ratio from mass and redshift
	INPUT :
        nmz: dn / (dmdz)
    			nmz.z, nmz.m
		param
            param['XC_0'], 
            param['SCATT_RC']  this is variance sigma**2  on core radius
            param['OMEGA_M'], param['SIGMA8'] not needed for computation

        model_struct: dictionary containing tags:
			cr: must be sorted in ascending order and the ratio of consecutive elements must be constant
			hr: must be sorted in ascending order and the ratio of consecutive elements must be constant
			rc: must be sorted in ascending order and the ratio of consecutive elements must be constant

	        crtab and hrtab have the same Temperature and Redshift vectors

        KEYED INPUT:
                scatter: dictionary containing a grid of scatter for luminosity and temperature
                verbose
                dndcrdhrdrcore
"""

def dndzdm_to_dndcrdhrdrcore_logz_scatter_mlmt_error(nmz, param, param1_name, param2_name, evol, cosmo, crtab, hrtab, sfunctab, \
    model_struct, area, threshold=0.005, scatter=None, verbose=False, dndcrdhrdrcore=None, andrea_scaling_laws=None, error=None):

    __version = 1

    start = time.time()
 
    if verbose and threshold == None: print "Threshold is not defined, I will assume 0.005 c/s"

    area_sqrad = area/math.degrees(1.)/math.degrees(1.)
    if verbose is not None: print 'area ', area_sqrad, 'steradians, threshold=', threshold, 'count/second '

    n_z = nmz['Z'].size
    dlogz = (np.log(nmz['Z'].max()) - np.log(nmz['Z'].min()))/(n_z-1)
    if verbose: print " n_z", n_z, dlogz

    n_m = nmz['M'].size/n_z
    diff = nmz['M'][0:n_m] - nmz['M'][n_m:2*n_m]
    flag = False
    if (diff.min() < 0)or (diff.max() > 0):
        n_m = nmz['M'].size
        flag = True

    if scatter is not None:
        n_t = scatter['VEC_T'].size
        n_l = scatter['VEC_L'].size

    # 03/02/16: RG Do not need to do n_m = nmz['M'].size//n_z
    # if it comes from python
    # if n_m is rebinned (comse from IDL) we divide otherwise no
    #n_m = nmz['M'].size//n_z  #  // for integer division

    dlogm = (np.log(nmz['M'].max()) - np.log(nmz['M'].min())) / (n_m - 1)
    if verbose: print 'n_m', n_m, 'dlogm=', dlogm, 'Python flag=', flag

    hr = model_struct['HR']
    cr = model_struct['CR']
    rc = model_struct['RC']

    nhr = hr.size
    ncr = cr.size
    nrc = rc.size

    dlogcr = (math.log(cr[-1])-math.log(cr[0])) / (ncr-1.0)
    dloghr = (math.log(hr[-1])-math.log(hr[0])) / (nhr-1.0)
    dlogrc = (math.log(rc[-1])-math.log(rc[0])) / (nrc-1.0)
    
    if (verbose):
        print 'ncr', ncr, 'nhr', nhr, 'nrc', nrc, 'dlogcr=', dlogcr,'dloghr=', dloghr, 'dlogrc=', dlogrc

    tab_transposed_cr = crtab['CR_PERL_TIMESD2'].T
    interpolator_cr = interpolate.interp2d(hrtab['Z'], hrtab['T'],  tab_transposed_cr)

    tab_transposed_hr = hrtab['CR_PERL_TIMESD2'].T
    interpolator_hr = interpolate.interp2d(hrtab['Z'], hrtab['T'],  tab_transposed_hr)

    func_count_rate = sfunctab['COUNTR']
    func_rcore = sfunctab['RCORE']
    func_proba = sfunctab['PROBA']

    da = evol['DA']  ## size n_m e.g. 200
    dl = evol['DL']  ## size n_m e.g. 200
    # Ez = evol.hc / cosmo.h0  code IDL
    Ez = evol['HC'] / cosmo['H0'] ## size n_m e.g. 200
    # In Nicolas code Ez was computed inside the loop on z

    # initialisation of the result
    if dndcrdhrdrcore is None:
        dndcrdhrdrcore = np.zeros([ncr, nhr, nrc])
    else:
        dndcrdhrdrcore.fill(0)
        #dndcrdhrdrcore[:] = 0.0   equivalent

    dndcrdhr = np.zeros([ncr, nhr])
    # z dependency

    # See below
    evol2 = evol

    # LF: THIS IS THE CRITICAL DENSITY TODAY DIVIDED BY 100^2 IN UNITS OF
    # SOLAR MASSES PER MEGAPARSEC CUBE
    # rhoc_0_e4 * 10000 IS THE CRITICAL DENSITY TODAY FOR H0 = 100 km / s / Mpc
    # THE CRITICAL DENSITY TODAY IS 2.7889537 x 10^11 h^2 Msun / Mpc^-3
    rhoc_0_e4 =  2.7889537E7

    probdet1 = np.zeros([ncr, nrc])

    for i in np.arange(0, n_z):
    
        dndcrdhr.fill(0)
        nm = extract_dict(nmz,i)
        z0 = nm['Z']
 
        # It was in Nicolas' code but it is not clear; Andrea says it is not necessary
        # evol2 = mk_evol(cosmo,ran_z=np.array([0,0.01])+z0,n_z=2)

        #  In Nicolas code Ez was computed here
        # LF: THIS IS THE CRITICAL DENSITY AT z
        # IT INCLUDES h (cosmo['H0']=0.72) SO THE UNITS ARE Msun / Mpc^-3, NOT h^2 Msun / Mpc^-3
        rhoc_z = rhoc_0_e4 * (Ez[i]*cosmo['H0']) ** 2.0

        nmconv200c = convert_nm(cosmo, evol2, nm, overdensity=200, rhoc=True)

        for  j in np.arange(0, n_m):
            # compute current point density ??
            delta_n = nm['DNDMDOMDZ'][j] * nm['M'][j] * dlogm * dlogz * z0 * area_sqrad

            if (delta_n <= 0):
                #if verbose: print 'delta_n zero', i, j, delta_n
                continue

            # compute current point temperature
            if andrea_scaling_laws is None:
                temperature_scalar = mass_to_temperature(nmconv200c['M'][j], z0, cosmo, evol2, param=param)
            else:
                # From IDL code: mstar:1d14*10.^param.norm_mt/0.7*cosmo.h IT IS WRONG !!!!!!!!!!!!!!!!!

                #if verbose: print "Applying Andrea's scaling laws"

                m_star = (1e14*10 ** param['NORM_MT']) / 0.7 * cosmo['H']
                mt_struct1 = {'POW': param['POW_MT'], 'MSTAR': m_star, 'HEVOL': param['HEVOL_MT'], \
                    'ZEVOL': param['ZEVOL_MT'],'SCATTER': param['SCATT_MT'], 'STRUCT_NAME':'mt'}
                mt_struct2 = check_param(mt_struct1, cosmo, 'mt', verbose = False) 

                ######################################
                # RG: bug in check_param ??
                # LF: put mt_struct2['MSTAR'] = m_star
                mt_struct2['MSTAR'] = m_star 
                ######################################

                temperature_scalar = mk_mt(cosmo, evol2, mt_struct2, mass=nmconv200c['M'][j],redshift=z0)

            # If scatter temperature is a vector of size n_t
            if scatter is not None:
                temperature = temperature_scalar * np.exp(scatter['VEC_T'])
                z_vec = np.tile(z0, n_t)
            else:
                temperature = temperature_scalar
                z_vec = np.array(z0)

            # compute current point count rate from luminosity
            if andrea_scaling_laws is None:
                luminosity_scalar = mass_to_luminosity(nmconv200c['M'][j], z0, cosmo, evol2, param=param)
            else:
                # From IDL code: lstar:1d44*10.^param.norm_lt

                #if verbose: print "Applying Andrea's scaling laws"

                l_star = 1e44*10 ** param['NORM_LT']
                lt_struct1 = {'POW': param['POW_LT'], 'LSTAR': l_star, 'HEVOL': param['HEVOL_LT'], \
                    'ZEVOL': param['ZEVOL_LT'], 'SCATTER': param['SCATT_LT'], 'SCATTER_EVOL': param['SCATTER_LT_EVOL'], 'STRUCT_NAME':'lt'}
                lt_struct2 = check_param(lt_struct1, cosmo, 'lt', verbose = False)

                luminosity_scalar = mk_lt(cosmo, evol2, lt_struct2, temp=temperature, redshift=z0)

            # Adding scattering on luminosity
            if scatter is not None:
                # Luminosity is now scattered from scalar to 1D vector of n_phys_l elements
                luminosity = luminosity_scalar * np.exp(scatter['VEC_L'])
                # Luminosity is now rebinned as a matrix of n_phys_t, n_phys_l elements
                luminosity = np.reshape(luminosity, (1, n_l))
            else:
                luminosity = luminosity_scalar
                # If no scatter maybe is necessary
  
            #count_rate = interpolator_cr(z0, temperature)
            #  beware this is not the count rate yet, have to multiply by luminosity
            this_cr = interpol2d(tab_transposed_cr, crtab['T'], crtab['Z'], temperature, z_vec)
            if scatter is not None:
                # this_hr, this_rc, this_cr are now rebinned as matrces of n_t, n_l elements
                this_cr = np.reshape(this_cr, (n_t, 1))
                    
            dimming = luminosity / 1.E44 / (np.square(dl[i]))
            this_cr = this_cr * dimming
            # if no scatter, this_cr is a scalar, else a matrix (n_t, n_l)

            #AV 14/01/2016: apply threshold   nselect=0 or continue?
            #if (this_cr < threshold):  continue
            
            # compute current point hardness ratio hr  does not depend upon Luminosity
            #  if no scatter, this_hr is a scalar, else vector 1d n_t
            this_hr=interpol2d(tab_transposed_hr, hrtab['T'], hrtab['Z'], temperature, z0 )
            if scatter is not None:
                this_hr = np.reshape(this_hr, (n_t, 1))
            
            # compute current point mass
            # Conversion from m200c to m500c maybe from Arnaud's 2010
            m500c = nmconv200c['M'][j] * 0.714

            # compute current point radius core rc
            # r500c = ( m500c*(.75/math.pi/ 500./rhoc_z / cosmo['H'] ))**(1/3.)
            #############################################################
            # BETTER WRITTEN
            # ALL MASSES SO FAR IN UNITS Msun / h
            # IT IS NOW TIME TO EXPLICITELY INSERT 1 / h TO TAKE OUT THIS DEPENDENCE
            # FOR rhoc_z ANY DEPENDENCE ON h HAS ALREADY BEEN TAKEN OUT SO r500c IS IN Mpc, NOT h^-1 Mpc
            r500c = ((3.0 * m500c / cosmo['H']) / (4.0 * math.pi * 500.0 * rhoc_z)) ** (1.0 / 3.0)
            #############################################################

            #  conversion from radian to arcseconds
            this_rc = r500c / da[i] * math.degrees(1) * 3600.0 * param['XC_0']
            #  this_rc is a scalar
            if scatter is not None:
                this_rc = np.reshape(this_rc, (1, 1))

            ########   Compute RC, HR  ############
            if scatter is not None:
                # MAYBE IT'S TRANSPOSE ??????
                # TRY WITH TRANSPOSE
                ###############################################
                # LF 05/04/16: I THINK IT SHOULD BE TRANSPOSE
                ###############################################
                nselect = delta_n * scatter['DUMMY'].T
                ###############################################

                #nselect = delta_n * scatter['DUMMY']
            else:
                nselect = np.array([delta_n])

            ii = np.round((np.log(this_cr)-np.log(cr[0])) / dlogcr)

            jj = np.round((np.log(this_hr)-np.log(hr[0])) / dloghr)
            if scatter is not None:
                jjj=np.tile(jj, [1, n_l])
            else:
                jjj = jj

            kk = np.round((np.log(this_rc)-np.log(rc[0])) / dlogrc)
            if scatter is not None: 
                kkk = np.tile(kk, [n_t, n_l])
            else:
                kkk = kk

            # AV 14/01/2016: apply threshold to this_cr
            index = np.where(this_cr >= threshold)
            if len(index[0]) == 0: continue
            ii, jjj, kkk, nselect = our_remove(index, ii, jjj, kkk, nselect)

            index = np.where(nselect > 0)
            if len(index[0]) == 0: continue
            ii, jjj, kkk, nselect = our_remove(index, ii, jjj, kkk, nselect)

            index = np.where((ii >= 0) & (ii < ncr))
            if len(index[0]) == 0: continue
            ii, jjj, kkk, nselect = our_remove(index, ii, jjj, kkk, nselect)

            index = np.where((jjj >= 0) & (jjj < nhr))
            if len(index[0]) == 0: continue
            ii, jjj, kkk, nselect = our_remove(index, ii, jjj, kkk, nselect)

            index = np.where((kkk >= 0) & (kkk < nrc))
            if len(index[0]) == 0: continue
            ii, jjj, kkk, nselect = our_remove(index, ii, jjj, kkk, nselect)

            ok = nselect.size
            if scatter is not None and param['SCATT_RC'] != 0:
                scatter_rc = param['SCATT_RC'] ## BEWARE  this is variance sigma**2
                profil = 1./math.sqrt(2*math.pi*scatter_rc)*np.exp(-0.5*((np.log(rc)-np.log(this_rc))**2/scatter_rc))*dlogrc
                profil=profil.flatten()
                diff = np.abs(profil.sum()-1.)
                #if (diff > 0.1):
                #    print 'Warning: profile is not one', profil.sum(), i, j, this_rc

                ##  normalisation
                profil = profil / profil.sum()

                for v in range(ok):
                    #probdet = compute_selfun(sfunctab, rc, this_cr[v])
                    nselect2 = profil * nselect[v] # * probdet
                    dndcrdhrdrcore[ii[v], jjj[v], :] += nselect2
            else:
                for v in range(ok):
                    probdet = compute_selfun(sfunctab, this_rc, this_cr[v])
                    dndcrdhrdrcore[ii[v], jjj[v], kkk[v]] += nselect[v] * probdet

    # Now apply probdet
    if scatter is not None:
       probdet_interpolator = interp2d(np.log(func_count_rate), np.log(func_rcore), func_proba, kind='linear', bounds_error=False)
       probdet = probdet_interpolator(np.log(cr), np.log(rc))
       probdet2 = (probdet.T).reshape([ncr, 1, nrc])
       dndcrdhrdrcore2 = dndcrdhrdrcore * probdet2
    else:
       dndcrdhrdrcore2 = dndcrdhrdrcore

    if error is not None:
        temporary = {'DNDCRDHRDRCORE':dndcrdhrdrcore2, 'CR':cr, 'HR':hr, 'RCORE':rc}
        intermediate_result = convol_3d_fix_error(temporary, error=error)
        dndcrdhrdrcore3 = intermediate_result['DNDCRDHRDRCORE']
    else:
        dndcrdhrdrcore3 =  dndcrdhrdrcore2

    result = {'DNDCRDHRDRCORE':dndcrdhrdrcore3, 'CR':cr, 'HR':hr, 'RCORE':rc, \
                param1_name:param[param1_name], param2_name:param[param2_name],\
                'Z':nmz['Z'], 'M':nm['M'], \
                'DNDCRDHRDRC_NO_ERR':dndcrdhrdrcore2, 'DNDCRDHRDRC_NO_ERR_NO_SEL':dndcrdhrdrcore}

    end = time.time()

    print "dndzdm_to_dndcrdhrdrcore_logz_scatter_mlmt_error, Time", end - start, "sec, version", __version

    return result
