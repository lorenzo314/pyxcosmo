import numpy as np
from interpol2d import interpol2d
import math
import pdb


def convol_2d_variable(original,errmod, out_cr=None,out_hr=None, accuracy=1e-5,  verbose=False):

    """
    PURPOSE:
    -------------
    convolve by ta model of gassiane error

    INPUTS:
    -------------
    (original) : a dictionnary containing:
    (CR, HR) grid where the function has been evaluated
    must be in increasing order, 2 vectors 1D
    DNDCRDHR= the values of the function on the regular grid, matrix 2D
    ---
    (errmod)


    OUTPUT:
    -------------
    Z: the convolved 

    HISTORY:
    -------------
    Rene Gastaud, 7 february 2015

    """

    in_cr = original['CR']
    in_hr = original['HR']
    n_incr = in_cr.size
    n_inhr = in_hr.size

    assert (n_incr, n_inhr) == original['DNDCRDHR'].shape


    ###
    log_cr = np.log(in_cr)
    steps_logcr =  log_cr[1:] - log_cr[:-1]
    delta_step_logcr = steps_logcr.max() - steps_logcr.min()
    if ( delta_step_logcr > accuracy):
            print "warning non regular grid on log cr ", delta_step_logcr 
    step_logcr =  (log_cr[-1] -log_cr[0])/(n_incr-1.)
    if (verbose): 
        print 'steps_logcr', steps_logcr.min(), step_logcr, steps_logcr.max(), delta_step_logcr
    binsize_cr = in_cr * step_logcr

    ###
    log_hr = np.log(in_hr)
    steps_loghr =  log_hr[1:] - log_hr[:-1]
    delta_step_loghr = steps_loghr.max() - steps_loghr.min()
    if ( delta_step_loghr > accuracy):
            print "warning non regular grid on log hr ", delta_step_loghr 
    step_loghr =  (log_hr[-1] -log_hr[0])/(n_incr-1.)
    if (verbose): 
        print 'steps_loghr', steps_loghr.min(), step_loghr, steps_loghr.max(), delta_step_loghr
    binsize_hr = in_hr * step_loghr

    ###
    #bins = binsize_cr.reshape(n_incr,1)*binsize_hr.reshape(1, n_inhr)
    bins_T = binsize_hr.reshape(n_inhr,1)*binsize_cr.reshape(1, n_incr)

    ## zoom of (1,1)
    if (out_cr == None): out_cr = in_cr
    if (out_hr == None): out_hr = in_hr
    n_outcr = out_cr.size
    n_outhr = out_hr.size

    mid_cr = np.exp(errmod['LOG_CR_TRUE'])
    mid_hr = np.exp(errmod['LOG_HR_TRUE'])
    ### python / IDL matrices are transposed
    mean_log_cr_tab_T = errmod['MEAN_LOG_CR'].T
    mean_log_hr_tab_T = errmod['MEAN_LOG_HR'].T
    var_log_crcr_tab_T = errmod['VAR_LOG_CRCR_TAB'].T
    var_log_hrhr_tab_T = errmod['VAR_LOG_HRHR_TAB'].T
    var_log_crhr_tab_T = errmod['VAR_LOG_CRHR_TAB'].T
    dndcrdhr = original['DNDCRDHR'].T

    totconv = np.zeros([n_outcr,n_outhr])

    (xx,yy) = np.meshgrid(in_cr, in_hr)
    log_xx = np.log(xx)
    log_yy = np.log(yy)

    ######################
    for i in range(0, n_incr):
        for j in range(0, n_inhr):

            cr = in_cr[i]
            hr = in_hr[j]

            isolate_point =dndcrdhr[i,j]*binsize_cr[i]*binsize_hr[j]

            if(isolate_point > 1e-10): 
                #print 'sdfsdf'
                expected_cr = interpol2d(mean_log_cr_tab_T, mid_cr, mid_hr, cr, hr)
                expected_hr = interpol2d(mean_log_hr_tab_T, mid_cr, mid_hr, cr, hr)

                covmatrix = np.zeros([2,2])
                covmatrix[0,0] = interpol2d(var_log_crcr_tab_T, mid_cr, mid_hr, cr, hr)
                covmatrix[1,1] = interpol2d(var_log_hrhr_tab_T, mid_cr, mid_hr, cr, hr)
                covmatrix[0,1] = interpol2d(var_log_crhr_tab_T, mid_cr, mid_hr, cr, hr)
                covmatrix[1,0] = covmatrix[0,1]

                detcov = np.linalg.det(covmatrix)
                invcov = np.linalg.inv(covmatrix)
                vecx = log_xx - expected_cr
                vecy = log_yy - expected_hr
                ugauss = invcov[0,0]*vecx**2 + 2*invcov[0,1]*vecx*vecy + invcov[1,1]*vecy**2
                gauss = 1./(2*math.pi*math.sqrt(detcov)) * np.exp(-.5*ugauss)
                thisconv = isolate_point * (gauss*step_logcr*step_loghr)
                totconv += thisconv
                #if ((i==61) and (j==45)): 
                    #print'debug', detcov,i,j
                    #pdb.set_trace()

    totconv /= bins_T
    result = {'DNDCRDHR':totconv, 'CR':out_cr, 'HR':out_hr}
    return result
