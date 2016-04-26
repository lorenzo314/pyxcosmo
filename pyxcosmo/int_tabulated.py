#Date: 19 March 2015
#Author : R Gastaud
# http://stackoverflow.com/questions/14345001/idls-int-tabulate-scipy-equivalent

import numpy as np
from scipy.integrate import newton_cotes 
def int_tabulated(x, f, p=5) :
    def mynewton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        #weights = scipy.integrate.newton_cotes(rn)[0]
        weights = newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1) :
        ret += mynewton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret
