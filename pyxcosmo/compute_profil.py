import numpy as np
import math

def compute_profil(mu, sigma2, xgrid, dlogx):
    """
        Compute the log normal law.

        INPUTS
        mu : the mu parameter of log normal,
        sigma2  : relative variance, variance/mu**2
        xgrid : vector of x for wich the probability density function is computed
        dlogx : step in the log space of xgrid

        OUTPUTS
        log normal probability density function

        EXAMPLE
        np.log(this_rc), scatter_rc, dlogrc, rc.shape
        (array([[ 3.50541506]]), 0.22314355131420976, 0.12122168464114909, (32,))
        
        profil = compute_profil(np.log(this_rc), scatter_rc, rc, dlogrc)

        R Gastaud fecit,  24 oct 2014

    """
    if ( sigma2 > 0):
        profil = np.exp(-0.5*(mu - np.log(xgrid))**2/sigma2)
        profil = profil/math.sqrt(2*math.pi*sigma2)*dlogx
    else:
        index = (np.abs(np.exp(mu)-xgrid)).argmin()
        profil = np.zeros( xgrid.size)
        profil[index] = 1
    

    return  profil

