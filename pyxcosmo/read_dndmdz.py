# Call:
# 	map=read_dndmdz('cosmo_grid_idl_andrea')

import numpy as np
import pyfits

from fitsutils import mrdfits
from calc_likelihood_dndmdz import calc_likelihood_dndmdz

def get_nm(filename):
    f = pyfits.open(filename)
    h = f[1].header
    f.close()

    k = h.values().index('M')
    last = (h.keys()[k])[-1]

    nm = int((h['TFORM'+last])[:-1])

    return nm

def read_dndmdz(outdir, background=1e-12, prefix="nmz_py", postfix='.fits', nn=15, verbose=False):
  if outdir[-1] != '/': outdir = outdir + '/'

  # Read the model; it is at the center
  i = nn / 2
  j = nn / 2
  suffix = 'sigma_%2.2d_omega_%2.2d' % (i, j)
  model_file = outdir + prefix + suffix + postfix
  model = mrdfits(model_file)
  n_m = get_nm(model_file)

  if verbose:
      print 'model', i, j, suffix
      print ' '

  likeit = np.zeros([nn, nn])
  #print nn, type(nn)
  for i in range(nn):
    for j in range(nn):
      suffix = 'sigma_%2.2d_omega_%2.2d' % (i, j)
      model_file = outdir + prefix + suffix + postfix
      dndmdz = mrdfits(model_file, 1)

      likelihood = calc_likelihood_dndmdz(model, dndmdz, background=background, n_m = n_m)
 
      likeit[i,j] = likelihood
      if verbose: print i, j, likelihood

  return likeit
