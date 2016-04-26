#!/usr/bin/env python
# 25 March 2014  rene gastaud

import numpy as np

def calc_prior(prior, limits):
   
    nn = len(prior)
    probas = np.zeros(nn)
    i=0
    for key in prior.keys():
        if ((prior[key] > limits[key][0])and(prior[key] < limits[key][1])):
            probas[i] = 1./(limits[key][1]-limits[key][0])
        i = i+1

    # special case sccatt_mt
    j = prior.keys().index('SCATT_MT')    
    if(prior[ "SCATT_MT"]==-1.): probas[j]=1.

    proba = probas.prod()
    return np.log(proba)
