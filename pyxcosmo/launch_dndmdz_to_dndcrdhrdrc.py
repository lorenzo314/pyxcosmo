from fitsutils import mrdfits
from fitsutils import mwrfits
import numpy as np
import matplotlib.pyplot as plt

import sys
import os

from plotimagediff import plotImageDiff

from pyxcosmo.dndmdz_to_dndcrdhrdrc import dndmdz_to_dndcrdhrdrc

cmap = 'gist_stern' # loadct, 15
cmap = 'rainbow'    # loadct, 34

vmin = 1
vmax = 400

icosmo = os.environ.get('ICOSMO')
if icosmo == None:
    print "Environment variable ICOSMO not defined: I am leaving"
    sys.exit(-1)

indir = icosmo + '/pyxcosmo/test_data/'
outdir = icosmo + '/pyxcosmo/tmp_dir/'
figdir = outdir + '/figures/'

fiducfile = indir + 'right_Aardvark_cosmology_xc0.24_scatt_fid.fits'
crfile    = indir + 'grid_w0_xc0/CRb05-2.fits'
hrfile    = indir + 'grid_w0_xc0/12R051.fits'
sfuncfile = indir + 'grid_w0_xc0/selfunc_official.fits'

crmin=0.005
crmax=3.5
rcmin=3.5
rcmax=150.0
hrmin=0.1
hrmax=2.0
ncr=64
nhr=64
nrc=64

scale = 1.0e5
save = True

nmz = mrdfits(indir + 'nmz_idl.fits')
dndcrdhrdrcore = np.ones([ncr, nhr, nrc])

results_py = dndmdz_to_dndcrdhrdrc(nmz,fiducfile,crfile,hrfile,sfuncfile, \
    crmin=crmin,crmax=crmax,ncr=ncr,rcmin=rcmin,rcmax=rcmax,nrc=nrc,hrmin=hrmin,hrmax=hrmax,nhr=nhr, \
    dndcrdhrdrcore=dndcrdhrdrcore)

for key in results_py.keys(): print key

outfile = outdir + 'dndcrdhrdrcore_py.fits'
if os.path.isfile(outfile): os.remove(outfile)
mwrfits(results_py, outfile)

# Read Andrea's IDL results
results_idl=mrdfits(indir + 'dndcrdhrdrcore_idl.fits')
for key in results_idl.keys(): print key, results_idl[key].shape

cube_py = results_py['DNDCRDHRDRCORE']

# T because IDL and python have different conventions
# and mrdfits does corrects images but not cubes
cube_idl = results_idl['DNDCRDHRDRCORE'].T

diff = cube_py - cube_idl
diff.min(), diff.max()

####### CR HR ##############################
plotImageDiff(cube_py.sum(2), cube_idl.sum(2), cmap, 'CR-HR', scale=scale, save=save, outdir=figdir)

####### CR RC ##############################
plotImageDiff(cube_py.sum(1), cube_idl.sum(1), cmap, 'CR-RC', scale=scale, save=save, outdir=figdir)

####### HR RC ##############################
plotImageDiff(cube_py.sum(0), cube_idl.sum(0), cmap, 'HR-RC', scale=scale, save=save, outdir=figdir)

# Show
plt.show()
