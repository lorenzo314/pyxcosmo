#!/usr/bin/env python
# 25 March 2014  rene gastaud

import numpy as np

def calc_prior(prior):
    limits = dict()
    limits[ "OMEGA_M"]  =[    0.0700000,      1.00000,      0.00000,      1.07527]
    limits[ "SIGMA8"]   =[    0.0500000,      2.00000,      0.00000,     0.512820]
    limits[ "POW_MT"]   =[      0.00000,      3.00000,      0.00000,     0.333333]
    limits[ "NORM_MT"]  =[     -2.00000,      2.00000,      0.00000,     0.250000]
    limits[ "HEVOL_MT"] =[     -3.00000,      3.00000,      0.00000,     0.166667]
    limits[ "ZEVOL_MT"] =[     -4.00000,      4.00000,      0.00000,     0.125000]
    limits[ "POW_LT"]   =[      0.00000,      4.00000,      0.00000,     0.250000]
    limits[ "NORM_LT"]  =[     -2.00000,      2.00000,      0.00000,     0.250000]
    limits[ "SCATT_LT"] =[      0.00000,      1.00000,      0.00000,      1.00000]
    limits[ "HEVOL_LT"] =[     -3.00000,      3.00000,      0.00000,     0.166667]
    limits[ "ZEVOL_LT"] =[     -5.00000,      3.00000,      0.00000,     0.125000]
    limits[ "XC_0"]     =[      0.00000,     0.900000,      0.00000,      1.11111]
    limits[ "XC_H"]     =[     -3.00000,      3.00000,      0.00000,     0.166667]
    limits[ "XC_Z"]     =[     -3.00000,      3.00000,      0.00000,     0.166667]
    limits[ "SCATT_MT"] =[      0.00000,      1.00000,      0.00000,      1.00000]
    limits[ "N_PK"]     =[     0.800000,      1.10000,      0.00000,      3.33333]
    limits[ "TAU"]      =[    0.0700000,      1.01000,      0.00000,      1.06383]

    nn = len(prior)
    probas = np.zeros(nn)
    i=0
    for key in limits.keys():
        probas[i]=limits[key][2] if ((prior[key] <= limits[key][0])or(prior[key] >= limits[key][1])) else limits[key][3]
        i = i+1

    # special case sccatt_mt
    j = limits.keys().index('SCATT_MT')    
    if(prior[ "SCATT_MT"]==-1.): probas[j]=1.

    proba = probas.prod()
    return np.log(proba)
