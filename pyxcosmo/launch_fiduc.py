#!/usr/bin/env python
# 25 March 2014  rene gastaud

import numpy as np
from fitsutils import mrdfits

def launch_fiduc(fiducfile=None,option=[],value=[], verbose=False):
    """ Get the structure fiduc.
    if fiducfile is defined, read the fits files fiducfile.
    else use the default dictionary plus the corrected values given by option and value
    No new tags can be creted, and this option is not valid when fiducfile is given
    """
    if (fiducfile is None):
        if(verbose): print 'fiduc by default'
        fiduc = dict()
        fiduc['OMEGA_B']=     0.043306000
        fiduc['OMEGA_M']=      0.24934400
        fiduc['SIGMA8']=     0.787000
        fiduc['POW_MT']=      1.49000
        fiduc['NORM_MT']=     0.459614
        fiduc['SCATT_MT']=     -1.00000
        fiduc['HEVOL_MT']=     -1.00000
        fiduc['POW_LT']=      2.89000
        fiduc['NORM_LT']=     0.399358
        fiduc['SCATT_LT']=     0.663000
        fiduc['HEVOL_LT']=      1.00000
        fiduc['ZEVOL_LT']=      0.00000
        fiduc['ZEVOL_MT']=      0.00000
        fiduc['XC_0']=     0.100000
        fiduc['XC_H']=      0.00000
        fiduc['XC_Z']=      0.00000
        fiduc['W0']=     -1.00000
        fiduc['WA']=      0.00000
        fiduc['HUBBLE']=     0.724000
        fiduc['TAU']=    0.0890000
        fiduc['N_PK']=     0.961000
        fiduc['SCATTER_LT_EVOL']=      0.00000
        if (len(option) > 0):
            len_option = len(option)
            if (len_option != len(value)):
                raise Exception('launch_fiduc len(option)=%d != len(value) %d'%(len_option, len(value)))
            else:
                i=0
                for key in option:
                    if (key in fiduc.keys()):
                        fiduc[key]= value[i]
                    else:
                        raise Exception('launch_fiduc key %s is not in fiduc'%key)
                    i = i+1
    else:
        if(verbose): print 'read '+fiducfile
        fiduc = mrdfits(fiducfile)
 
    return fiduc
