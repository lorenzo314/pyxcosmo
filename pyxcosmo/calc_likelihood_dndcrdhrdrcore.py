import numpy as np

import pdb

def calc_likelihood_dndcrdhrdrcore(observable, model, bkg=1.0e-12):
    """
        Compute likelihood

        INPUTS
        observable, dictionary with tag 'DNDCRDHRDRCORE'
        model, dictionary with tag 'DNDCRDHRDRCORE'

        OUTPUTS
        Cash statistic (Cash 1979, ApJ 228, 939) 
          2 * Sum(Mi - DiLog(Mi))

        EXAMPLE

        R Gastaud fecit, 21 january 2016
        L Faccioli 30/02/16
        A Valotti 06/04/16 put bkg
    """

    model_dndcrdhrdrcore = model['DNDCRDHRDRCORE']
    observable_dndcrdhrdrcore = observable['DNDCRDHRDRCORE']

    # RG 31/03/16: CAREFUL ABOUT np.where(model_dndcrdhrdrcore > 0.0)
    # MAYBE IT'S BETTER TO USE A THRESHOLD SPECIFIED HERE
    #q = np.where(model_dndcrdhrdrcore > 0.0)
    #model_dndcrdhrdrcore = model_dndcrdhrdrcore[q]
    #observable_dndcrdhrdrcore = observable_dndcrdhrdrcore[q]

    model = model_dndcrdhrdrcore + bkg
    observable = observable_dndcrdhrdrcore + bkg
    res = 2 * (model - observable * np.log(model)).sum()

    return  res
