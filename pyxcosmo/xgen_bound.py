###
import numpy as np
import __builtin__
import math

def xgen_bound(x1, x2, npoints, log_flag=False):
    """
        Compute the grid

        INPUTS
        x1 = minimum
        x2 = maximum
       
        OPTIONNAL INPUTS
        log_flag


        OUTPUTS
        grid

        EXAMPLE
        vec_cr_bound = xgen_bound(cr_min, cr_max, n_cr, True)

        SEE ALSO
        compute_step, xgen

        R Gastaud fecit,  4 December 2014

    """

    if (log_flag):
        x1 = math.log10(x1)
        x2 = math.log10(x2)
        step = (x2 - x1)/(npoints-1)
        x1 = x1 - step/2.
        x2 = x2 + step/2.
        xx = np.logspace(x1, x2, num=(npoints+1), endpoint=True)
    else:
       step = (x2 - x1)/(npoints-1)
       x1 = x1 - step/2.
       x2 = x2 + step/2.
       xx = np.linspace(x1,x2, num=(npoints+1), endpoint=True)

    return xx
