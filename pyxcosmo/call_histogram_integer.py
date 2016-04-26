###
###  name must be different from histogram_integer.so
###  9 oct 2015 add check of the existance of the file histogram_integer.so

import numpy as np
import os.path

import numpy.ctypeslib as npct
from ctypes import c_int



def call_histogram_integer(arr1, arr2, size1, size2, size3, location):

    ##################
    # input type for the cos_doubles function
    # must be an integer array, with single dimension that is contiguous
    array_1d_int = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')
    
    # check the existence of the library
    filename =  location+"/histogram_integer.so"
    status = os.path.isfile(filename)
    if not(status):
        print filename+' not found by call_histogram_integer'
        return -1
    
    # load the library, using numpy mechanisms
    libcd = npct.load_library("histogram_integer", location)

    # setup the return typs and argument types
    libcd.histogram_integer.restype = None
    libcd.histogram_integer.argtypes = [array_1d_int, array_1d_int, c_int, c_int, c_int,array_1d_int ]

    n_ri = (size1*size2+size3+1)
    ri = np.zeros(n_ri, dtype=np.int64)
    libcd.histogram_integer(arr1, arr2, size1, size2, size3, ri )
    return ri
