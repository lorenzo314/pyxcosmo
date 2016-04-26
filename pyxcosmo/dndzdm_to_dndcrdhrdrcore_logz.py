import numpy as np

from scipy.interpolate import interp1d
from scipy import interpolate
import math
from scipy.ndimage.interpolation import map_coordinates

from pyxcosmo.get_cosmo import get_cosmo
from pyxcosmo.icosmo_cosmo import mk_evol
from pyxcosmo.mk_mt import mk_mt
from pyxcosmo.check_param import check_param
from pyxcosmo.mk_lt import mk_lt
from pyxcosmo.icosmo_cosmo import convert_nm
from pyxcosmo.extract_dict import extract_dict
from pyxcosmo.interpol2d import interpol2d

from compute_virtual_index import compute_virtual_index

#import pdb
#import scipy

from pyxcosmo.extract_dict import extract_dict

#   to be checked
#     def compute_probdet(proba, x, y, xnew, ynew):

"""
    Compute count rate and hardness ratio from mass and redshift
	INPUT :
    		nmz: dn / (dmdz)
    			nmz.z, nmz.m
		param

		cr: must be sorted in ascending order and the ratio of consecutive elements must be constant
		hr: must be sorted in ascending order and the ratio of consecutive elements must be constant
		rc: must be sorted in ascending order and the ratio of consecutive elements must be constant

	crtab and hrtab have the same Temperature and Redshift vectors
"""

def dndzdm_to_dndcrdhrdrcore_logz(nmz, param, evol, cosmo, crtab, hrtab, sfunctab, \
    cr, hr, rc, area, threshold=0.005, verbose=False, dndcrdhrdrcore=None):

    if verbose and area == None: print "Area is not defined, I will assume 10000.0 square degrees"
    if verbose and threshold == None: print "Threshold is not defined, I will assume 0.005 c/s"

    area_sqrad = area/math.degrees(1.)/math.degrees(1.)  ### 100*100 square degrees = 3.046 square radians


    n_z = nmz['Z'].size
    dlogz = (np.log(nmz['Z'].max()) - np.log(nmz['Z'].min()))/(n_z-1)
    if verbose: print " n_z", n_z, dlogz

    n_m = nmz['M'].size/n_z
    diff = nmz['M'][0:n_m] - nmz['M'][n_m:2*n_m]
    flag = False
    if (diff.min() < 0)or (diff.max() > 0):
        n_m = nmz['M'].size
        flag = True

    # 03/02/16: RG Do not need to do n_m = nmz['M'].size//n_z
    # if it comes from python
    # if n_m is rebinned (comse from IDL) we divide otherwise no
    #n_m = nmz['M'].size//n_z  #  // for integer division

    dlogm = (np.log(nmz['M'].max()) - np.log(nmz['M'].min())) / (n_m - 1)
    if verbose: print 'n_m', n_m, 'dlogm=', dlogm, 'Python flag=', flag

    nhr = hr.size
    ncr = cr.size
    nrc = rc.size

    dlogcr = (math.log(cr[-1])-math.log(cr[0])) / (ncr-1.0)
    dloghr = (math.log(hr[-1])-math.log(hr[0])) / (nhr-1.0)
    dlogrc = (math.log(rc[-1])-math.log(rc[0])) / (nrc-1.0)
    
    if (verbose):
        print 'ncr', ncr, 'nhr', nhr, 'nrc', nrc, 'dlogcr=', dlogcr,'dloghr=', dloghr, 'dlogrc=', dlogrc

    # From IDL code: mstar:1d14*10.^param.norm_mt/0.7*cosmo.h
    m_star = (1e14*10 ** param['NORM_MT']) / 0.7 * cosmo['H']

    mt_struct1 = {'POW': param['POW_MT'], 'MSTAR': m_star, 'HEVOL': param['HEVOL_MT'],'ZEVOL': param['ZEVOL_MT'],'SCATTER': param['SCATT_MT'], 'STRUCT_NAME':'mt'}

    # From IDL code: lstar:1d44*10.^param.norm_lt
    l_star = 1e44*10 ** param['NORM_LT']
    lt_struct1 = {'POW': param['POW_LT'], 'LSTAR': l_star, 'HEVOL': param['HEVOL_LT'],'ZEVOL': param['ZEVOL_LT'],     'SCATTER': param['SCATT_LT'], 'SCATTER_EVOL': param['SCATTER_LT_EVOL'], 'STRUCT_NAME':'lt'}

    mt_struct2 = check_param(mt_struct1, cosmo, 'mt', verbose=True)
    mt_struct2 ['MSTAR'] = m_star  ###  bug in check_param ??
    lt_struct2 = check_param(lt_struct1, cosmo, 'lt', verbose=True)

    if (verbose):
        print 'm_star', mt_struct2['MSTAR']-m_star, 'l_star', lt_struct2['LSTAR']-l_star

    tab_transposed_cr=crtab['CR_PERL_TIMESD2'].T
    interpolator_cr = interpolate.interp2d(hrtab['Z'], hrtab['T'],  tab_transposed_cr)
    
    tab_transposed_hr=hrtab['CR_PERL_TIMESD2'].T
    interpolator_hr = interpolate.interp2d(hrtab['Z'], hrtab['T'],  tab_transposed_hr)
    
    func_count_rate = sfunctab['COUNTR']
    func_rcore = sfunctab['RCORE']
    func_proba = sfunctab['PROBA']

    da = evol['DA']  ## size n_m e.g. 200
    dl = evol['DL']  ## size n_m e.g. 200
    # Ez = evol.hc / cosmo.h0  code IDL
    Ez = evol['HC'] / cosmo['H0'] ## size n_m e.g. 200
    #  In Nicolas code Ez was computed inside the loop on z

    # initialisation of the result
    if dndcrdhrdrcore is None:
        dndcrdhrdrcore = np.zeros([ncr, nhr, nrc])
    else:
        for i in range(ncr):
            for j in range(nhr):
                for k in range(nrc):
                    dndcrdhrdrcore[i, j, k] = 0.0

    ####
    # z dependency
    ####

    # See below
    evol2 = evol

    rhoc_0_e4 =  2.7889537E7

    for  i in np.arange(0, n_z):
        nm = extract_dict(nmz,i)
        z0 = nm['Z']
 
        # It was in Nicolas' code but it is not clear; Andrea says it is not necessary
        # evol2 = mk_evol(cosmo,ran_z=np.array([0,0.01])+z0,n_z=2)

        #  In Nicolas code Ez was computed here
        rhoc_z = rhoc_0_e4 * (Ez[i]*cosmo['H0']) ** 2.0

        nmconv200c = convert_nm(cosmo, evol2, nm, overdensity=200, rhoc=True)

        for  j in np.arange(0, n_m):
            # compute current point density ??
            nn = nm['DNDMDOMDZ'][j]*nm['M'][j]*dlogm*dlogz*z0*area_sqrad
            if (nn <= 0):
                if verbose: print 'nn zero', i, j, nn
                continue

            # compute current point temperature
            temperature = mk_mt(cosmo, evol2, mt_struct2, mass=nmconv200c['M'][j],redshift=z0)

            # compute current point mass
            # Conversion from m200c to m500c maybe from Arnaud's 2010
            m500c = nmconv200c['M'][j] * 0.714

            # compute current point count rate from luminosity
            luminosity = mk_lt(cosmo, evol2, lt_struct2, temp=temperature, redshift=z0)

            #count_rate = interpolator_cr(z0, temperature)
            count_rate = interpol2d(tab_transposed_cr, crtab['T'], crtab['Z'], temperature, z0 )
            count_rate = count_rate.flatten()
            dimming = luminosity / 1.E44 / (np.square(dl[i]))
            this_cr = count_rate*dimming

            #AV 14/01/2016: apply threshold   nselect=0 or continue?
            #if (this_cr < threshold):  continue
            
            # compute current point hardness ratio hr
            this_hr=interpol2d(tab_transposed_hr, hrtab['T'], hrtab['Z'], temperature, z0 )

            # compute current point radius core rc
            r500c = ( m500c*(.75/math.pi/ 500./rhoc_z / cosmo['H'] ))**(1/3.)

            #  conversion from radian to arcseconds
            this_rc = r500c / da[i] * math.degrees(1) * 3600.0 * param['XC_0']

            # probdet
            flag = True
            virtual_iy = compute_virtual_index(np.log(func_rcore), np.log(this_rc), flag=flag)
            virtual_ix = compute_virtual_index(np.log(func_count_rate), np.log(this_cr), flag=flag)
            coords = np.array([virtual_iy, virtual_ix])
            coords2 = coords.reshape(2,1)
            probdet = map_coordinates(func_proba, coords2, prefilter=True, order=1)

            nselect = nn * probdet
            if(nselect == 0): continue

            #AV 14/01/2016: apply threshold   nselect=0 or continue?
            if (this_cr < threshold):  continue
            
            ii = round((math.log(this_cr)-math.log(cr[0])) / dlogcr)
            if (ii < 0): continue
            if (ii >= ncr): continue

            jj = round((math.log(this_hr)-math.log(hr[0])) / dloghr)
            if (jj < 0): continue
            if (jj >= nhr): continue

            kk = round((math.log(this_rc)-math.log(rc[0])) / dlogrc)
            if (kk< 0): continue
            if (kk >= nrc): continue

            dndcrdhrdrcore[ii, jj, kk] = dndcrdhrdrcore[ii, jj, kk] + nselect

    result = {'DNDCRDHRDRCORE':dndcrdhrdrcore, 'CR':cr, 'HR':hr, 'RCORE':rc, \
                'OMEGA_M':param['OMEGA_M'], 'SIGMA8':param['SIGMA8'],\
                'Z':nmz['Z'], 'M':nm['M']}

    print "dndzdm_to_dndcrdhrdrcore version 111"
    return result
