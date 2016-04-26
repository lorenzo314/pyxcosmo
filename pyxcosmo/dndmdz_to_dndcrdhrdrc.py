from fitsutils import mrdfits
from fitsutils import mwrfits
import numpy as np

from scipy.interpolate import interp1d
from scipy import interpolate
import math

from pyxcosmo.get_cosmo import get_cosmo

from pyxcosmo.icosmo_cosmo import mk_evol
from pyxcosmo.mk_mt import mk_mt
from pyxcosmo.check_param import check_param
from pyxcosmo.mk_lt import mk_lt
from pyxcosmo.icosmo_cosmo import convert_nm

from compute_virtual_index import compute_virtual_index

import pdb
import scipy

from extract_dict import extract_dict

"""
    Compute count rate, hardness ratio, core radius from mass and redshift

    INPUT:
      nmz: dictionary containing the number of clusters for a given mass and redshift
      It comes from make_nm(mass, z)
        z: e.g. one value from 0.05 to 1.8 log scale, 100 steps
        mass: e.g. 1e14 to 1e16 by log scale, 200 steps

      fiducfile: file containing the fiducial cosmological model

      crfile: file containing crtab
        crtab: Table of precomputed values of CR as a function of mass and redshift 

      hrfile: file containing hrtab
        hrtab: Table of precomputed values of HR as a function of mass and redshift
 
      sfuncfile: file containing the selection function

      see dndzdldt_to_dndcrdhr

    For now nmz is computed by a program from Nicolas
"""

def dndmdz_to_dndcrdhrdrc(nmz, fiducfile, crfile, hrfile, sfuncfile, \
  crmin=0.005,crmax=3.5,ncr=64,rcmin=3.5,rcmax=150.0,nrc=64,hrmin=0.1,hrmax=2.0,nhr=64):

    n_z = nmz['Z'].size
    dlogz = (np.log(nmz['Z'].max()) - np.log(nmz['Z'].min()))/(n_z-1)

    param = mrdfits(fiducfile,1)

    # get the cosmo
    ol = 1-param['OMEGA_M']
    oc = param['OMEGA_M'] - param['OMEGA_B']
    (cosmo, evol1) = get_cosmo(wmap5=True, omega_b=param['OMEGA_B'], omega_c=oc, omega_l=ol,
                               sigma8=param['SIGMA8'], w_l=param['W0'], w1_l=param['WA'], h=param['HUBBLE'],
                               n_pk=param['N_PK'], tau=param['TAU'])

    m_star = (1e14*10**param['NORM_MT'])/0.7*0.72

    mt_struct1 = {'POW': param['POW_MT'], 'MSTAR': m_star, 'HEVOL': param['HEVOL_MT'],'ZEVOL': param['ZEVOL_MT'],'SCATTER': param['SCATT_MT'], 'STRUCT_NAME':'mt'}

    lt_struct1 = {'POW': param['POW_LT'], 'LSTAR': 1e44*10**param['NORM_LT'], 'HEVOL': param['HEVOL_LT'],'ZEVOL': param['ZEVOL_LT'], 'SCATTER': param['SCATT_LT'], 'SCATTER_EVOL': param['SCATTER_LT_EVOL'], 'STRUCT_NAME':'lt'}
    mt_struct2 = check_param(mt_struct1, cosmo, 'mt', verbose=True)
    lt_struct2 = check_param(lt_struct1, cosmo, 'lt', verbose=True)

    mt_struct2['MSTAR'] = m_star

    crtab = mrdfits(crfile)
    hrtab = mrdfits(hrfile)

    # LF 07/12/15
    # crtab and hrtab have the same Temperature and Redshift vectors
    # so one may use either crtab or hrtab in both instructions below, but for
    # correctness it's better not to
    tab_transposed_cr=crtab['CR_PERL_TIMESD2'].T
    interpolator_cr = interpolate.interp2d(crtab['Z'], crtab['T'],  tab_transposed_cr)
    
    tab_transposed_hr=hrtab['CR_PERL_TIMESD2'].T
    interpolator_hr = interpolate.interp2d(hrtab['Z'], hrtab['T'],  tab_transposed_hr)
    
    # This can be computed with m_ran and n_m

    n_z = nmz['Z'].size

    # LF 07/12/15: // needed because it must be an integer division
    n_m = nmz['M'].size//n_z
    dlogm = (np.log(nmz['M'].max()) - np.log(nmz['M'].min()))/(n_m-1)
    
    sfunctab = mrdfits(sfuncfile,1)
    func_count_rate = sfunctab['COUNTR']
    func_rcore = sfunctab['RCORE']
    func_proba = sfunctab['PROBA']

    cr_out = np.zeros([n_m, n_z])
    hr_out = np.zeros([n_m, n_z])
    rc_out = np.zeros([n_m, n_z])
    nselect = np.zeros([n_m, n_z])

    nn = np.zeros([n_m, n_z])

    for  i in np.arange(0, n_z):
        print i
        nm = extract_dict(nmz,i)
        z0 = nm['Z']
        evol2 = mk_evol(cosmo,ran_z=np.array([0,0.01])+z0,n_z=2)
        interpolator = interp1d(evol2['Z'], evol2['HC'])
        Ez = interpolator(z0)/cosmo['H0']

        nmconv200c=convert_nm(cosmo, evol2, nm, overdensity=200, rhoc=True)

        # Following Arya
        # m500c is a vector of n_m masses 
        m500c = nmconv200c['M']*0.714

        rhoc_0_e4 =  2.7889537E7
        rhoc_z = rhoc_0_e4*(Ez*cosmo['H0'])**2.
        temperature = mk_mt(cosmo, evol2, mt_struct2, mass=nmconv200c['M'],redshift=z0)
        lum = mk_lt(cosmo, evol2, lt_struct2, temp=temperature, redshift=z0)
        r500c = ( m500c*(.75/math.pi/ 500./rhoc_z / cosmo['H'] ))**(1/3.)
        hr_out[:,i] = (interpolator_hr(z0, temperature)).flatten()

        dim = lum/1.E44/(np.square(evol2['DL'][0]))
        # 0 because i=0 first z?

        count_rate = interpolator_cr(z0, temperature)
        count_rate=count_rate.flatten()
        cr_out[:,i]  = count_rate*dim

        da = evol2['DA'][0]
        #  conversion from radian to arcseconds
        rc_out[:,i] = r500c /da*math.degrees(1)*3600.*param['XC_0']

        nn[:,i] = nm['DNDMDOMDZ']*nm['M']*dlogm*dlogz*z0/math.degrees(1.)/math.degrees(1.)*1e4

    # probdet
    flag = True
    virtual_iy = compute_virtual_index(np.log(func_rcore), np.log(rc_out), flag=flag)
    virtual_ix = compute_virtual_index(np.log(func_count_rate), np.log(cr_out), flag=flag)
    coords = np.array([virtual_iy, virtual_ix])
    probdet = scipy.ndimage.map_coordinates(func_proba, coords, prefilter=True, order=1)
    probdet = probdet.reshape([n_m, n_z])
    nselect = nn*probdet

    # Temporary, for debugging puprposes only
    # result = {'CR':cr_out, 'HR':hr_out, 'RCORE':rc_out, 'Z':nmz['Z'], 'M':nm['M'], 'NSELECT':nselect,'PROBDET':probdet, 'NN':nn }

    # Compute the cube
    dlogcr = (math.log(crmax) - math.log(crmin)) / (ncr - 1.0)
    dloghr = (math.log(hrmax) - math.log(hrmin)) / (nhr - 1.0)
    dlogrc = (math.log(rcmax) - math.log(rcmin)) / (nrc - 1.0)

    dcrs =  math.exp(dlogcr / 2.0)
    dhrs =  math.exp(dloghr / 2.0)
    drcs =  math.exp(dlogrc / 2.0)

    i1 = (cr_out > crmin/dcrs) & (cr_out < crmax*dcrs)
    i2 = (hr_out > hrmin/dhrs) & (hr_out < hrmax*dhrs)
    i3 = (rc_out > rcmin/drcs) & (rc_out < rcmax*drcs)
    iii = i1&i2&i3

    v_ii = (np.log(cr_out[iii]) - math.log(crmin)) / dlogcr
    v_jj = (np.log(hr_out[iii]) - math.log(hrmin)) / dloghr
    v_kk = (np.log(rc_out[iii]) - math.log(rcmin)) / dlogrc

    v_ii = np.round(v_ii)
    v_jj = np.round(v_jj)
    v_kk = np.round(v_kk)

    v_ii.min(), v_ii.max() ,  v_jj.min(), v_jj.max() ,  v_kk.min(), v_kk.max()

    dndcrdhrdrcore = np.zeros([ncr, nhr, nrc])
    weights = nselect[iii]
    n_w = weights.size

    for k in np.arange(0, n_w):
        dndcrdhrdrcore[v_ii[k], v_jj[k], v_kk[k]] += weights[k]

    cr_grid=np.linspace(math.log(crmin), math.log(crmax), num=ncr, endpoint=True)
    cr_grid=np.exp(cr_grid)

    hr_grid=np.linspace(math.log(hrmin), math.log(hrmax), num=nhr, endpoint=True)
    hr_grid=np.exp(hr_grid)

    rc_grid=np.linspace(math.log(rcmin), math.log(rcmax), num=nrc, endpoint=True)
    rc_grid=np.exp(rc_grid)

    result = {'CR':cr_grid, 'HR':hr_grid, 'RCORE':rc_grid, 'DNDCRDHRDRCORE':dndcrdhrdrcore, 'Z':nmz['Z'], 'M':nm['M']}

    return result
