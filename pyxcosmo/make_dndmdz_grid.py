# from fitsutils import mrdfits
from fitsutils import mwrfits
from make_dndmdz import make_dndmdz
from pyxcosmo.launch_fiduc import launch_fiduc
from pyxcosmo.compare_idlpy import compare_idlpy

# fiducfile = '../../pyxcosmo/test_data/right_Aardvark_cosmology_xc0.24_scatt_fid.fits'

def make_dndmdz_grid(fiducfile, outdir, sigma8, omega_m, verbose=False):
    z_ran=[0.05,1.8]
    n_z=100
    m_ran = [1e12,1e16]
    n_m = 200

    param=launch_fiduc(fiducfile = fiducfile)

    n_omega = omega_m.size
    n_sigma = sigma8.size
    print n_omega, n_sigma

    if outdir[-1] != '/': outdir = outdir + '/'

    for i in range(n_sigma):
        for j in range(n_omega):
            suffix = 'sigma_%2.2d_omega_%2.2d' % (i, j)
            param['SIGMA8']  = sigma8[i]
            param['OMEGA_M'] = omega_m[j]

            dndmdz = make_dndmdz(param, m_ran=m_ran, n_m=n_m, z_ran=z_ran)
            mwrfits(dndmdz, outdir + 'nmz_py_' + suffix+'.fits', clobber = True)
            #mwrfits(param, outdir + 'nmz_py_' + suffix+'.fits')

            if verbose:
                print i, j, 'param', 'omega_m=', omega_m, 'sigma8=', sigma8, ' ', suffix

    print "make_dndmdz_grid finished"
