# LF, RG: 02/02/16
import numpy as np

from pyxcosmo.icosmo_cosmo import xgen
from pyxcosmo.icosmo_cosmo import mk_evol
from pyxcosmo.icosmo_cosmo import mk_nm
from pyxcosmo.get_cosmo import get_cosmo

def make_dndmdz(param, z_ran=[0.05,1.8], n_z=100, m_ran=[1e12,1e16], n_m=200):
    ol = 1 - param['OMEGA_M']
    oc = param['OMEGA_M'] - param['OMEGA_B']

    cosmo, evol1 = get_cosmo(wmap5=1, omega_b=param['OMEGA_B'], omega_c=oc, omega_l=ol, \
        sigma8=param['SIGMA8'], w_l=param['W0'], w1_l=param['WA'], h=param['HUBBLE'], \
        n_pk=param['N_PK'], tau=param['TAU'])

    z = xgen(z_ran[0], z_ran[1], n_z, logplot = True)

    # RG, LF 02/02/16: we want to have logarithmic z bins so we recall mk_evol again
    # and store the result in another variable called evol2
    evol2 = mk_evol(cosmo, z=z)

    zbins=evol2['Z']

    nmz = mk_nm(cosmo, evol2, z=np.array(zbins), m_ran=m_ran, n_m=n_m, fit_nm=4, \
        tinkerdelta=200, fit_tk=1, ctype=0, profile=True)

    return nmz, cosmo, evol2
