#!/usr/bin/env python
#  15 March 2015  use is None
#  14 March 2015  tau
#  13 March 2015  n_pk
# 2 february 2015 R Gastaud
#
#import numpy as np

from rd_cosmo import rd_cosmo
from icosmo_cosmo import mk_evol

def get_cosmo( h=None,omega_r=0.,omega_c=None,omegac_h2=None,omega_l=None, omega_b=None,omegab_h2=None,w_l=None,w1_l=None,age=None,
               n_pk=None,Pk_feature=None,sigma8=None,gamma=None,tau=None,
               wmap1=None,wmap3=None,wmap5=None,wmap7=None,deltar2=None,kpiv=None,verbose=False,
               ran_z=[0.,5.], n_z=200, z=None):

    status = 1
    myconfig = {'wmap7':0, 'wmap5':1,  'wmap3':2, 'wmap1':3}
    h_tab = [.71, .724, .734, .71]
    omega_b_tab = [0.0449, 0.043306, 0.0432, 0.044]
    omega_c_tab = [0.222,  0.206038, 0.1948, 0.226]
    omega_l_tab = [0.734,  0.750656, 0.762,  0.73]
    n_pk_tab = [.963, .961, .951, .93]
    sigma8_tab = [.801, .787, .744, .84]
    tau_tab = [0.088, 0.089, 0., 0.]

    mode = 'unknown'
    if not(wmap7 is None): mode='wmap7' 
    if not(wmap5 is None): mode='wmap5' 
    if not(wmap3 is None): mode='wmap3' 
    if not(wmap1 is None): mode='wmap1'
    index = myconfig.get(mode)
    if (verbose): print 'mode=', mode, index

    if (index != None): 
        if (h == None): h = h_tab[index] 
        if ((omega_b == None) and (omegab_h2 == None)): omega_b = omega_b_tab[index]
        if ((omega_c == None) and (omegac_h2 == None)): omega_c = omega_c_tab[index]
        if (omega_l == None) : omega_l = omega_l_tab[index]
        if (n_pk == None): n_pk = n_pk_tab[index] 
        if ((sigma8 == None) and (deltar2 == None)): sigma8 = sigma8_tab[index]
        if (tau == None): tau = tau_tab[index] 

    #if(verbose):
    #    print ' yy h, omega_c, omega_b, omega_l, omega_r, n_pk, sigma8, tau'
    #    print  h, omega_c, omega_b, omega_l, omega_r, n_pk, sigma8, tau

    cosmo  = rd_cosmo(h=h,omega_c=omega_c,omegac_h2=omegac_h2,omega_b=omega_b,omegab_h2=omegab_h2,
                      omega_l=omega_l,omega_r=omega_r,w_l=w_l,w1_l=w1_l,age=age,
                      n_pk= n_pk, pk_feature=Pk_feature,sigma8=sigma8,deltar2=deltar2,kpiv=kpiv,gamma=gamma,tau=tau)
    if (verbose): 
        print "cosmo=", cosmo
        print 'ran_z', ran_z, 'n_z=', n_z, 'z=',z

    evol = mk_evol(cosmo, ran_z=ran_z, n_z=n_z, z=z)

    return cosmo, evol

