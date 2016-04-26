# get scalar argument  get_scalar_argument

import numpy as np

def get_int_argument(argument):
    if not np.isscalar(argument): 
        # this for the case where shape does not exist
        # don't ask me why we can have no shape with a table of 1 element
        argument = argument.reshape(argument.size)
        argument = argument[0]
    argument = int (argument)
    return  argument


def get_scalar_argument(argument):
    if not np.isscalar(argument): 
        # this for the case where shape does not exist
        # don't ask me why we can have no shape with a table of 1 element
        argument = argument.reshape(argument.size)
        argument = argument[0]
    return  argument
