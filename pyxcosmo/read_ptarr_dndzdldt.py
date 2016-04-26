# IPython log file

#get_ipython().magic(u'logstart')
import os.path
from fitsutils import mrdfits

def read_ptarr_dndzdldt(in_dir, pattern='dndzdldt_', verbose=False):
    """
        Compute the virtual indexes.

        INPUTS
        directory

        OUTPUTS
        data_list 

        EXAMPLE
        data_list = read_ptarr_dndzdldt('init_cmp_dndcr_out')

	HISTORY:
        R Gastaud fecit,  8 May 2015

    """
    if not(in_dir.endswith(os.path.sep)): in_dir=in_dir+os.path.sep
    result = []
    all_red_shifts = mrdfits(in_dir+'all_red_shifts.fits.gz',0)
    comb = mrdfits(in_dir+'comb.fits.gz',0)
    nn = comb.size
    result.append(all_red_shifts)
    result.append(comb)

    for i in range(2 ,nn+2):
        suffix = "{0:02d}".format(i)
        if verbose:print 'read', i, in_dir+pattern +suffix+'.fits.gz'
        ref_data =  mrdfits(in_dir+ pattern+suffix+'.fits.gz')
        result.append(ref_data)

    return result
