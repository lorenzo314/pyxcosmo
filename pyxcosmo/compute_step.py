###
import numpy as np
import __builtin__
#import math

def compute_step(x, nn=None, log_flag=False):
    """
        Compute the increasing step of an array x.
        It is the average of the consecutive difference.
        If log_flat set, it is the log of the ratio.

        INPUTS
        the array of values
       
        OPTIONNAL INPUTS
        nn size of the array 
        log_flag


        OUTPUTS
        step

        EXAMPLE
        index = compute_virtual_index(grid, values)

        SEE ALSO
        xgen, xgen_bound

        R Gastaud fecit,  4 December 2014

    """

    nx = x.size
    if (nn != None):
        if (nn != nx):
            print 'error nx != nn', nx, nn

    if (log_flag):
        x = np.log(x)  
        # this is a new vector, so we loose the pointer to the old, 
        #  no return outside the function

    steps = x[1:]-x[0:-1]
    diff = steps.max()-steps.min()
    if (diff > 8e-7):
        print 'error in compute_step', diff
        return -1

    step=steps.mean()  # or median ???

    return  step
