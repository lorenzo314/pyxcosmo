########
# 25 June 2015 add some comments
#
# 28 April 2015 mk_evol, default value of z is none 
#   chi_a transform a in a numpy array version 121
# 16 April 2015 growth_a verbose version 120
# 23 March 2015 R Gastaud version 119 
#   flags fit_tk and fit_nl put to zero by default and not False
#   mk_del2 if fit_tk != 0 raise an exception
#   mk_nm  use cumtrapz  instead of trapz
# 18 March 2015 R Gastaud version 118 growth_a_ode precision and verbose
# 13 March 2015 R Gastaud version 117 mk_evol np.sin() 
# 10 March 2015 R Gastaud version 116 mk_nm really corrected  (transpose)
#   and mk_del2 (transpose)
#   copied from the version of icosmo_cosmo in brouillon5, i 18 janvier 2012 
#
#  8 March 2015 R Gastaud version 115 with bug in mk_nm corrected  nm
#
# 27 February 2015  R Gastaud version 114 with bug in mk_del2 corrected
#
# 16 March 2016 L Faccioli, R Gastaud: convert_nm now uses extrap1d
# 11/04/2011 N.CLERC
#
########

import numpy as np
from math import *
import math as math
import string
import scipy.integrate
import scipy.interpolate

from pylab import size
import pdb 
import inspect

__version__=122

############
# NAME: hubble_a
#
# HISTORY:
# 	May 05 - output turned into a structure by AR
# 	Jan 05 - Written by A. Refregier
# 	April 2011 - Python : N CLERC
# 	12 May 2011 : Rene Gastaud iron out the little bugs.
#
# PURPOSE: 
# compute the hubble constant H(a) as a function of the
# expansion parameter a and related basic evolution parameters. This routine
# is meant to centralise all the basic dark energy evolution dependence.
#
# INPUT: 
#	a: expansion parameter, here a nympy 1d array  
#		BEWARE a must be doubles, see a.dtype
#       cosmo: cosmological parameter structure
#
# OUTPUT: 
#	hubble_a: structure containing:
# 		H(a): Hubble constant [km/s/Mpc]
#		w_a: dark energy equation of state parameter w(a) [1]  ==> -1 ??
#		qevol: evolution of dark energy density    ==> 1 ???
#                                  rho_l(a)/rho_l(0)   [1]
#             	h2evol: (H(a)/H_0)^2  [1]
#            	dlnh_dlna:  d(ln H)/d(ln a)  [1]
#
############

def hubble_a(cosmo,a):

    # w(a) in Linder-Polarsky parametrisation
    w_a = cosmo["W_L"]+(1.-a)*cosmo["W1_L"]

    # recast a scalar in numpy array
    a = np.asarray(a)

    # rho_l(a)/rho_l(0)
    # first exponent is zero a**0=1 and second exponent is also zero exp(0)=1
    # result is [1,1,]
    if (cosmo["W_L"] == -1) and (cosmo["W1_L"] == 0):
         qevol = np.ones(a.size)
    else:
        qevol=a**(-3.*(1. +cosmo["W_L"]+cosmo["W1_L"]))*np.exp(-3. *cosmo["W1_L"]*(1. -a))  

    # (H(a)/H_0)^2  [1]
    # I introduce ia to speed up the code of a factor 3.6, so now as quick as bad idl code 
    ia  = 1./a
    ia2 = ia*ia
    ia3 = ia2*ia
    ia4 = ia2*ia2
    h2evol  = cosmo["OMEGA_R"]*ia4   + cosmo["OMEGA_M"]*ia3 
    h2evol += cosmo["OMEGA_L"]*qevol + cosmo["OMEGA_K"]*ia2
    
    # H(a) [km/s/Mpc] 
    hc = cosmo["H0"]*np.sqrt(h2evol)
    
    # ; d(ln H)/d(ln a)  [1]

    #dlnh_dlna=-.5/h2evol*(4.*cosmo["OMEGA_R"]*a**(-4.)+3.*cosmo["OMEGA_M"]*a**(-3.)+3.*cosmo["OMEGA_L"]*qevol*(1.+w_a)+2.*cosmo["OMEGA_K"]*a**(-2.))
    dlnh_dlna  = 4.*cosmo["OMEGA_R"]*ia4            + 3.*cosmo["OMEGA_M"]*ia3
    dlnh_dlna += 3.*cosmo["OMEGA_L"]*qevol*(1.+w_a) + 2.*cosmo["OMEGA_K"]*ia2
    dlnh_dlna = -.5/h2evol*dlnh_dlna

    hubble = {"W_A":w_a,"HC":hc,"QEVOL":qevol,"H2EVOL":h2evol,"DLNH_DLNA":dlnh_dlna}
    return hubble

def growth_a(cosmo, a, unnorm=False, type='ODE'):
    """ 
    NAME: growth_a
    HISTORY:
            Sep. 05 - Written by F. Pacaud using A. Refregier sub-procedures
            April 2011 Python : N CLERC
            12 May 2011 : Rene Gastaud iron out the little bugs.

    PURPOSE: 
            compute the linear growth factor D as a function of the
            expansion parameter a, for a given cosmology. This is done using fitting
            functions published in the litterature. The default normalisation is such 
            that D(z=0)=1, but it can be set such that D(a)=a in a matter dominated 
            universe (omega_m=omega=1)

    INPUT: 
            a: expansion parameter
            cosmo: cosmological parameter structure

    KEYWORD: 
            unnorm: normalise such that D(a)=a for omega_m=omega=1
                      universe instead default of D(z=0)=1

    gtype: growth type, used by growth_a, legal values:
            'ODE' for growth_a_ode   resolving differential equation
            'INT' for d=growth_a_int numerical integration of the Heath equation     
            'ANA' for growth_a_ana   analytical formulae

    OUTPUT: 
            D(a): growth factor as a function of a [1]
    """

    if (abs(cosmo["OMEGA_M"]-1.) < 1.e-3) and (abs(cosmo["OMEGA_K"]) < 1.e-3):
        type='INT'
    mtype= type.upper()
    #print 'growth_a '+mtype
    done = False
    if (mtype == 'ODE'):
        result=growth_a_ode(cosmo, a, unnorm=unnorm)
        done = True
    elif (mtype == 'INT'):
        result=growth_a_int(cosmo, a, unnorm=unnorm)
        done = True
    elif (mtype == 'ANA'):
        result=growth_a_ana(cosmo, a, unnorm=unnorm)
        done = True
    
    if (not(done)):
        print "method %s not implemented" % mtype
        raise Exception('growth_a type=%s not yet implemented'%mtype)
        
    return result

################
# NAME: growth_a_ana
#
# HISTORY:
#   Jan 05 - Written by A. Refregier
#   April 2011 Python : N CLERC
#   12 May 2011 : Rene Gastaud iron out the little bugs.
#
# PURPOSE: 
#   compute the linear growth factor D as a function of the
#   expansion parameter a, for a given cosmology. This is done using fitting
#   functions published in the litterature. The default normalisation is such 
#   that D(z=0)=1, but it can be set such that D(a)=a in a matter dominated 
#   universe (omega_m=omega=1)
#
# INPUT: 
#   a: expansion parameter
#   cosmo: cosmological parameter structure
#
# KEYWORD:
#   unnorm: normalise such that D(a)=a for omega_m=omega=1 universe instead default of D(z=0)=1
#
# OUTPUT: 
#   D(a): growth factor as a function of a [1]
################
def growth_a_ana(cosmo, a, unnorm=False):
    if(cosmo["W1_L"] != 0.):
        print "model with varying w not supported"
        sys.exit(1)

    hubble = hubble_a(cosmo,a)
    h_h0=hubble["HC"]/cosmo['H0']    # H/H_0
    omega_m_a=cosmo["OMEGA_M"]*a**(-3)*h_h0**(-2)
    omega_l_a=cosmo["OMEGA_L"]*hubble["QEVOL"]*h_h0**(-2)
    omega_a=omega_m_a+omega_l_a
    omega_k_a=1.-omega_a

    # Linear growth factor from analytical formulae (for w=constant only)
    # (from Ma 1998, APJL, 508, L5, who cites lahav et al. 1991 and 
    # Caroll et al. 1992)
    if (cosmo["OMEGA_L"] == 0. or (cosmo['W_L'] == -1)):  # not a Q model
        # at z =0:
        omm=cosmo["OMEGA_M"] 
        oml=cosmo["OMEGA_L"]
        g0=2.5*omm/( omm**(4./7.) -oml + (1.+.5*omm)*(1.+oml/70.) )

        # at z<>0
        omma=omega_m_a
        omla=omega_l_a
        g=2.5*omma/( omma**(4./7.) -omla + (1.+.5*omma)*(1.+omla/70.) )
        d_an=a*g/g0          # normalised, i.e. D(z=0)=1
        d_an_u=a*g           # unnormalised

    else:
    # w=cste model, use results from Ma et al. 1999 ApJL 521,1
        # check if the model is supported
        if (cosmo["W1_L"] != 0. or cosmo["OMEGA"] != 1.):
            print 'mk_evol: error: varying w or non flat quintessence mode not supported'
            sys.exit(1)

        # first, calculate growth factor for associated Lambda model
        omm=cosmo["OMEGA_M"]
        oml=cosmo["OMEGA_L"]                        # z=0
        g0_l=2.5*omm/( omm**(4./7.) -oml + (1.+.5*omm)*(1.+oml/70.) )    
        hc_l=cosmo["H0"]*sqrt(cosmo["OMEGA_M"]*a**(-3)+cosmo["OMEGA_L"]+cosmo["OMEGA_K"]*a**(-2))
        h_h0_l=hc_l/cosmo["H0"]    # H/H_0
        omma_l=cosmo["OMEGA_M"]*a**(-3)*h_h0_l^(-2)
        omla_l=cosmo["OMEGA_L"]*h_h0_l**(-2)
        g_l=2.5*omma_l/( omma_l**(4./7.) -omla_l + (1.+.5*omma_l)*(1.+omla_l/70.) )

        # apply quintessence corrections from Eq. (5) in Ma et al. 1999
        w=cosmo['W_L']
        t0=-(.255+.305*w+.0027/w)*(1.-omm)-(.366+.266*w-.07/w)*np.log(omm)
        g0=g0_l*(-w)^t0
        omma=omega_m_a
        t=-(.255+.305*w+.0027/w)*(1.-omma)-(.366+.266*w-.07/w)*np.log(omma)
        g=g_l*(-w)^t
        d_an=a*g/g0    # normalised (D(z=0)=1)
        d_an_u=a*g     # unnormalised

    # return growth factor with desired normalisation
    if(unnorm):
        d=d_an_u
    else:
        d=d_an

    return d

###########################
# NAME: t0_tilde
#
# HISTORY:
#	Jan 05 - Written by A. Refregier
# 	April 2011 Python : N CLERC
# 	12 May 2011 : Rene Gastaud iron out the little bugs.
#
# PURPOSE:
#
# INPUT:
# 	Q an image
#	alpha scalar
#	beta scalar
# OUPUT
#	an image
#
def t0_tilde(Q, ALPHA, BETA):
	ee=exp(1.)
	#print ee, log(ee)
	c=14.2/ALPHA+386./(1.+69.9*Q**1.08)
	mylog=np.log(ee + 1.8*BETA*Q)
	sortie = mylog/(mylog+c*Q**2)
	return sortie

#########################
# NAME: tk.py 
#
# HISTORY:
# 	12 May 2011 : Rene Gastaud iron out the little bugs.
# 	Apr 2011 : Python NCLERC
#
# 	Nov 2007 - Updated by AR to compute the Eisenstein & Hu transfer function
#          for a treatement of the baryons
# 	March 1999 - Written by A. Refregier
#
# PURPOSE: 
# compute the (unnormalised) linear transfer function T(k) for a
# given cosmological model and including baryonic corections. This is computed 
# using the analytical fits of Eisenstein and Hu (1998, ApJ, 496, 605)
#
# INPUT: 
#	k: comoving wave number [h Mpc^-1] a 1d or 2d vector
#	cosmo: cosmological parameter (a dictionary)
#
# OPTIONAL INPUT: 
#	fit_tk: fitting function to use: 
#		0: E&H without wiggles (default), 
#		1: E&H with wiggles, 
#		2: BBKS as summarised by Peacock & Dodds (1997) 
#
# OUTPUT:
# 	T(k) from Eisenstein & Hu (1998)
#
#########################
def tk(k,cosmo,fit_tk=0):
    if (fit_tk == 2):
        # use P&D
        # transfer function from Peacock and Dodds (1996, mnras, 280, L19)
        # as summarized by Peacock (1997, mnras, 284, 885)
        q_pd=k/cosmo['GAMMA']
        tk_out=np.log(1.+2.34*q_pd)/(2.34*q_pd)*( 1.+3.89*q_pd+(16.1*q_pd)**2+(5.46*q_pd)**3+(6.71*q_pd)**4 )**(-.25)

    else:
        # use E&H
        # declarations
        t_cmb=2.728   # [K] from Fixsen et al. (1996)
        theta_27=t_cmb/2.7
        omm=cosmo['OMEGA_M']
        omb=cosmo['OMEGA_B']
        h=cosmo['H']

        # compute various scales
        kk=k*h                 # k in [Mpc^-1]
        z_eq=2.50e4*omm*h**2.*theta_27**(-4.)  # matter/radiation equ. [1]
        k_eq=7.46e-2*omm*h**2/theta_27**2.    # [Mpc^-1]
        b1=0.313*(omm*h**2)**(-0.419)*(1.+0.607*(omm*h**2)**0.674)
        b2=0.238*(omm*h**2)**0.223
        z_d=1291.*(omm*h**2)**0.251 / (1.+0.659*(omm*h**2)**0.828)*(1.+b1*(omb*h**2)**b2)
        r_eq=31.5*omb*h**2./theta_27**4.*1e3/z_eq
        r_d =31.5*omb*h**2./theta_27**4.*1e3/z_d
        s=2./(3.*k_eq)*np.sqrt(6./r_eq)* np.log((np.sqrt(1.+r_d)+np.sqrt(r_d+r_eq))/(1.0+np.sqrt(r_eq)))
        k_silk=1.6*(omb*h**2)**0.52*(omm*h**2)**0.73 * (1.+(10.4*omm*h**2)**(-0.95)) # [Mpc^-1]

        if (fit_tk == 1):# E&H with wiggles
            # transfer functions parameters
            a1=(46.9*omm *h**2)**0.670*(1.0+(32.1*omm*h**2)**(-0.532))
            a2=(12.0*omm*h**2)**0.424*(1.0+(45.0*omm*h**2)**(-0.582))
            alpha_c=a1**(np.double(-omb/omm))*a2**(-(omb/omm)**3.)
            bb1=0.944/(1.0+(458.*omm*h**2)**(-0.708))
            bb2=(0.395*omm*h**2)**(-0.0266)
            beta_c=1./(1.0+bb1*((cosmo['OMEGA_DM']/omm)**bb2-1.0))
            
            # compute transfer function
            q=kk/(13.41*k_eq)   # dimensionless wave number

            # CDM part of T(k)
            f=1./(1.0+(kk*s/5.4)**4.)  # interpolation factor
            tk_c=f*t0_tilde(q,1.,beta_c)+(1. -f)*t0_tilde(q,alpha_c,beta_c)
      
            # Baryonic part of T(k)
            beta_node=8.41*(omm*h**2)**0.435
            s_tilde=s/(1.0+(beta_node/(kk*s))**3.)**(1./3.)
            y=(1.0+z_eq)/(1.0+z_d)
            gy=y*(-6.*sqrt(1.0+y)+(2.0+3.*y)* np.log((np.sqrt(1.0+y)+1.0)/(np.sqrt(1.0+y)-1.)))
            alpha_b=2.07*k_eq*s*(1.0+r_d)**(-3./4.)*gy
            beta_b=0.5+omb/omm+(3.0-2.*omb/omm)*np.sqrt((17.2*omm*h**2.)**2. +1.) 
            # rg decomposed  q is an image, t0_tilde(q,1.,1.) is an image
            #print "kk.shape", kk.shape, " s", s.shape, " alpha_b", alpha_b.shape, " beta_b", beta_b.shape
            #print "k_silk", k_silk.shape, " s_tilde", s_tilde.shape
            #  kk (100, 200)  s (1,)  alpha_b (1,)  beta_b (1,)  k_silk (1,)  s_tilde (100, 200)
            tk_b=(t0_tilde(q,1.,1.) / (1.0+(kk*s/5.2)**2.) + alpha_b/(1.0+(beta_b/kk/s)**3.)*np.exp(-(kk/k_silk)**1.4))
            tk_b *= np.sin(kk*s_tilde)/(kk*s_tilde)
            ##bid=check_math(mask=32)
   
            # total transfer function
            tk_out=omb/omm*tk_b+cosmo['OMEGA_DM']/omm*tk_c
        else:  # E&H without wiggles
            ee=exp(1.)
            alpha_gamma=1.0-0.328*log(431.*omm*h**2)*omb/omm+0.38*log(22.3*omm*h**2)*(omb/omm)**2
            gamma_eff=omm*h*(alpha_gamma+(1.0-alpha_gamma)/(1.0+(0.43*kk*s)**4))
            q=k*theta_27**2/gamma_eff    # rescaled wave number
            c0=14.20+731./(1.0+62.5*q)
            l0=np.log(2.*ee+1.8*q)
            tk_out=l0/(l0+c0*q**2)
	
    return tk_out

#################################
# chi_a and chi_a_int
#
# Apr 11 : N.CLERC Python
# Jan 05   - modified by AR to call hubble_a.pro
# Dec 04   - modified by AR to include varying dark energy eq. of state
# March 99 - Written by A. Refregier
#
# PURPOSE: compute the comoving distance Chi(a) as a function of the
# expansion parameter a, for a given cosmology. This is done by numerical
# integration
# INPUT: a: expansion parameter
#        cosmo: cosmological parameter structure
# OUTPUT: chi(a): comoving distance in unit of the hubble radius
#                 i.e. chi/R_0 is given rather than chi
##################################
def chi_a_int(a, cosmo):
    #  cosmo["W_L"]
    #  cosmo["W1_L"]
    #  cosmo['H0']
    #  cosmo["OMEGA_M"]
    #  cosmo["OMEGA_L"]
    #  cosmo["OMEGA_K"]
    a = np.asarray(a)
   
    #  rho_l(a)/rho_l(0)
    # first exponent is zero a**0=1 and second exponent is also zero exp(0)=1
    #  result is [1,1,]
    if (cosmo["W_L"] == -1) and (cosmo["W1_L"] == 0):
         #print "zero"
         qevol = np.ones(a.size)
         result = np.sqrt(cosmo["OMEGA_M"] + cosmo["OMEGA_L"]*a*a*a + cosmo["OMEGA_K"]*a)
    else:
        #print "1"
        qevol=a**(-3.*(1. +cosmo["W_L"]+cosmo["W1_L"]))*np.exp(-3. *cosmo["W1_L"]*(1. -a))  
        result = np.sqrt(cosmo["OMEGA_M"] + cosmo["OMEGA_L"]*qevol*a*a*a + cosmo["OMEGA_K"]*a)

    result = np.sqrt(cosmo["OMEGA_M"] + cosmo["OMEGA_L"]*qevol*a*a*a + cosmo["OMEGA_K"]*a)

    result = -cosmo["SQRTK"]/np.sqrt(a)/result

    return result

def chi_a(cosmo,a):
    a = np.asarray(a) # for the loop to compute integrale1 27 April 2015
    cosmo2 = copy_cosmo(cosmo)
    myargs = (cosmo2,) ## now myargs is a tuple containing as argument only one argument, the dictionary
    integrale1 = [scipy.integrate.romberg(chi_a_int, 1., x, args=myargs, vec_func=True) for x in a]
    # now transform the list in a numpy array
    integrale1 = np.asarray(integrale1)
    integrale1 = integrale1.ravel()

    return integrale1

#################################
# t_a
#
# Apr 11 : N.CLERC Python
# Jan 05 - modified by AR to call hubble_a.pro
# Dec 04   - modified by AR to include varying dark energy eq. of state
# March 99 - Written by A. Refregier
#
# PURPOSE: compute the look back time as a function of the
# expansion parameter a, for a given cosmology. This is done by numerical
# integration of t=int_1^a da'/(a'H(a')) 
# INPUT: a: expansion parameter
#        cosmo: cosmological parameter structure
# OUTPUT: t(a): look back time as a function of a in units of [Ho^-1]
#################################
def t_a_int(a, cosmo):
    #  cosmo["W_L"]
    #  cosmo["W1_L"]
    #  cosmo['H0']
    #  cosmo["OMEGA_M"]
    #  cosmo["OMEGA_L"]
    #  cosmo["OMEGA_K"]

    a = np.asarray(a)
   
    #  rho_l(a)/rho_l(0)
    # first exponent is zero a**0=1 and second exponent is also zero exp(0)=1
    #  result is [1,1,]
    if (cosmo["W_L"] == -1) and (cosmo["W1_L"] == 0):
         #print "zero"
         qevol = np.ones(a.size)
         result = np.sqrt(cosmo["OMEGA_M"] + cosmo["OMEGA_L"]*a*a*a + cosmo["OMEGA_K"]*a)
    else:
        #print "1"
        qevol=a**(-3.*(1. +cosmo["W_L"]+cosmo["W1_L"]))*np.exp(-3. *cosmo["W1_L"]*(1. -a))  
        result = np.sqrt(cosmo["OMEGA_M"] + cosmo["OMEGA_L"]*qevol*a*a*a + cosmo["OMEGA_K"]*a)
    #
    result = np.sqrt(a)/result
    return result

######################################
# t_a
######################################
def t_a(cosmo,a):
    cosmo2 = copy_cosmo(cosmo)
    myargs = (cosmo2,) ## now myargs is a tuple containing as argument only one argument, the dictionary
    # we want a tabulated primitive, not a definite integral
    #t = -scipy.integrate.romberg(t_a_int,1.,a,tol=1.e-6,vec_func=True)

    if np.isscalar(a):
         myintegrales =  -scipy.integrate.romberg(t_a_int, 1., a, args=myargs, vec_func=True)
    else:
        myintegrales = np.zeros(a.size)
        for i, x in enumerate(a):
            myintegrales[i] = -scipy.integrate.romberg(t_a_int, 1., x, args=myargs, vec_func=True)
    
    return myintegrales

######################################
# growth_a_intgd
# Apr 11: Python N CLERC
# Dec 04 - Written by A. Refregier
#
# PURPOSE: compute the linear growth factor D as a function of the
#  expansion parameter a, for a given cosmology. This is done by numerical
# integration of the Heath equation D= H*int_1^a da'/(a'H(a'))^3. 
#  http://ned.ipac.caltech.edu/level5/Carroll/Carroll3_5.html
#  Heath 1977 
# The default
# normalisation is such that D(z=0)=1, but it can be set such that
# D(a)=a in a matter dominated universe (omega_m=omega=1).
# NOTE: this equation is only valid for w=-1 models
# http://ned.ipac.caltech.edu/level5/Carroll/Carroll3_5.html
# INPUT: a: expansion parameter
#        cosmo: cosmological parameter structure
# KEYWORD: unnorm: normalise such that D(a)=a for omega_m=omega=1
#                  universe instead default of D(z=0)=1
# OUTPUT: D(a): growth factor as a function of a [1]
#
#######################################
def growth_a_intgd(a, cosmo):
    #  cosmo["W_L"]
    #  cosmo["W1_L"]
    #  cosmo['H0']
    #  cosmo["OMEGA_M"]
    #  cosmo["OMEGA_L"]
    #  cosmo["OMEGA_K"]

    a = np.asarray(a)
 
    #  rho_l(a)/rho_l(0)
    # first exponent is zero a**0=1 and second exponent is also zero exp(0)=1
    #  result is [1,1,]
    if (cosmo["W_L"] == -1) and (cosmo["W1_L"] == 0):
         #print "zero"
         qevol = np.ones(a.size)
         result = np.sqrt(a)/(cosmo['H0']*np.sqrt(cosmo["OMEGA_M"]+cosmo["OMEGA_L"]*a*a*a+cosmo["OMEGA_K"]*a))
    else:
        #print "1"
        qevol=a**(-3.*(1. +cosmo["W_L"]+cosmo["W1_L"]))*np.exp(-3. *cosmo["W1_L"]*(1. -a))  
        result = np.sqrt(a)/(cosmo['H0']*np.sqrt(cosmo["OMEGA_M"]+cosmo["OMEGA_L"]*qevol*a*a*a+cosmo["OMEGA_K"]*a))
    #
   
    result = result*result*result

    return result

#######################################
###   growth_a_int
#######################################
def growth_a_int(cosmo,a,unnorm=False):
    
    cosmo2 = copy_cosmo(cosmo)
    myargs = (cosmo2,) ## now myargs is a tuple containing as argument only one argument, the dictionary
     
    hubble = hubble_a(cosmo2,a)

    #integrale1 = np.zeros(a.size)
    #for i, x in enumerate(a):
    #    integrale1[i] = scipy.integrate.romberg(growth_a_intgd, 0., x, args=myargs, vec_func=True)

    integrale1 = [scipy.integrate.romberg(growth_a_intgd, 0., x, args=myargs, vec_func=True) for x in a]
    # now transform the list in a numpy array
    integrale1 = np.asarray(integrale1)
    integrale1 = integrale1.ravel()

    if(unnorm):  
        growth=5.*cosmo2["OMEGA_M"]/2.*hubble["HC"]*integrale1*cosmo2["H0"]**2
        # normalised so that D(a)=a for om_m=om_tot=1 model
    else:
        integrale2 = scipy.integrate.romberg(growth_a_intgd,0.,1.,args=myargs, vec_func=True)
        growth=hubble["HC"]*integrale1/(cosmo2["H0"]*integrale2)
        # normalised so that D(z=0)=1

    return growth

#######################################
###  beware dy/dt = f(y,t) and not f(t,y) 
def growth_a_derivs( y, x,  cosmo):
    # x can be a scalar
    # print "x=", x, "y=",y
    x = np.asarray(x)
    #  x=ln(a)
    a = np.exp(x)  
    ## a = exp(x)  
    hubble = hubble_a(cosmo, a)  
    loga = np.log10(a)
    lacube  = -3. * loga
    grandeur = lacube.round()
    lacube   -= grandeur
    temph2evol = np.log10(hubble["H2EVOL"]) - grandeur

    lomega_m_a = lacube - temph2evol
    omega_m_a  = 10.**lomega_m_a
    omega_m_a *= cosmo["OMEGA_M"]

    #Numerically stable when called by odeint
    # Old was :
    #omega_m_alternatif = cosmo2.omega_m*a^(-3)/hubble.h2evol
    #print, abs(omega_m_alternatif-omega_m_a)/omega_m_alternatif

    #  compute derivatives of y1=G and y2=dG/dlna
    #  f1=dy1/dx=y2=dG/dlna
    f1=y[1]  
    # ; f2=dy2/dx=d^2G/d(lna)^2=-(4.+dlnH/dlna)y2-(3+dlnH/dlna-1.5Omega_m(a))y1                 
    f2=-(4. + hubble["DLNH_DLNA"])*y[1] - (3. + hubble["DLNH_DLNA"] - 1.5*omega_m_a)*y[0]
    #print "f1.shape", f1.shape, "f2.shape", f2.shape                            
    return np.array([f1,f2])
    #return f2

#######################################
### growth_a_ode
#######################################
#def growth_a_ode(cosmo, a, unnorm=False, a_min=1.e-6, eps=0.05, hwant=0.02 ):
def growth_a_ode(cosmo, a, unnorm=False, a_min=1.e-6, eps=0.003, hwant=0.001 ):
    # this is for debugging = who is calling me?
    #frm = inspect.stack()
    #print 'frm[2] ', frm[2]
    # print 'frm[3] ', frm[3]
    #pdb.set_trace()
    # http://stackoverflow.com/questions/1095543/get-name-of-calling-functions-module-in-python
    #print ' growth_a_ode', __version__, ' integration accuracy eps=', eps, ' hwant=', hwant   
 
    cosmo2 = copy_cosmo(cosmo)
    myargs = (cosmo2,) ## now myargs is a tuple containing as argument only one argument, the dictionary
    
    # initial condition for integration
    # min takes two arguments, np.min one argument, which is a numpy.array or a scalar.
    a_start = min(np.min(a), a_min)   
    x_start = log(a_start)
    x_end = 0. ####  and not log(np.max(a))   x_end=0 ==> 1/(1+z)=1 ==> z=0 ==> today (now) and not the past
    #  initial conditions: y1=G=1, y2=dG/dlna=0.
    y_start = [1., 1e-10]
    #eps     = .05          # integration accuracy
    #hwant   = .2           # desired (i.e.) maximum step in x
    # RG 8 december 2011 hwant=0.02 instead of 0.2
    h1      = hwant        # guess for first step in x
   
    #; compute G(a)=D(a)/a from a=a_start to a=1 by numerical integration for 
    #; the ODE system:
    #; x=ln(a), y1=G, y2=dG/dlna
    #; dy1/dx=dG/dlna=y2
    #; dy2/dx=d^2G/d(lna)^2=-(4+dlnH/dlna)y2-(3+dlnH/dlna-1.5*Omega_m(a))y1
    #; with initial conditions: y1=G=1 and y2=dG/dlna=0 for a->0

    # odeint,y_start,x_start,0.d,eps,h1,0.d,nok,nbad,'growth_a_derivs',xp,yp,hwant=hwant
    #print 'y_start=', y_start
    #print  'x_start=', x_start, 'x_end=', x_end
    #print "eps ", eps, h1
    t_array = np.linspace(x_start, x_end, ceil((x_end-x_start)/hwant))
    y_trajectory = scipy.integrate.odeint(growth_a_derivs, y_start, t_array, args=myargs)
    
    #print "y_trajectory0 shape",y_trajectory[:,0].shape
    #print "y_trajectory1 shape",y_trajectory[:,1].shape
    #print "t_array", t_array.shape
    #print "a", a.shape
    # interpolate to desired values and set normalisation
    linear_interp = scipy.interpolate.interp1d(t_array, y_trajectory[:,0])
    d = linear_interp(np.log(a)) *a
    # D(a) normalised so that D=a in a matter dominated universe (omega_m=omega=1)
    #print " normalisation", y_start[0], unnorm, y_trajectory[-1,0] 
    if not (unnorm): 
       d=d/y_trajectory[-1,0] 
       #print "norm"
    return d
    # D(a) normalised so that D

#######################################
def copy_cosmo(cosmo):
    
    cosmo2 = cosmo.copy()
    noms = cosmo2.keys()
    
    if not ('W_L' in noms):
        cosmo2['W_L'] = -1.

    # another way to do it is with has_key
    if not (cosmo2.has_key('W1_L')):
        cosmo2['W1_L'] = 0.

    return cosmo2

################
# NAME: cmp_deltac
#
# HISTORY:
#	Sep. 05 - Written by F. Pacaud using A. Refregier sub-procedures
# 	April 2011 Python : N CLERC
# 	12 May 2011 : Rene Gastaud iron out the little bugs.
#
# PURPOSE: 
# compute the density threshold delta_c and the mean
# halo density ratio Delta_c as function of z, omega_m and omega_l.
# This is done using appendix A in Kitayama & Suto, 1996, ApJ, 469, 480.
#
# INPUT: 
#	cosmo: cosmological parameter structure
#
# OPTIONAL INPUT: z: redshift (default z=0)
# KEYWORD: 
#	fof_length:
#
# OUTPUT: 
# deltac structure: 
#	dc: density threshold delta_c(z,Omega_m,omega_l) 
#                   (nu=dc/sigma(M,z) in P&S mass function)
#	ddc: density contrast Delta_c(Omega_m,omega_l)=rho(r<r_vir)/rho_bar
#
################
def cmp_deltac(cosmo, z=0., fof_length=0.2):
    omm = cosmo["OMEGA_M"]
    oml = cosmo["OMEGA_L"]
    z1 =  1.+z
    z2 = z1*z1
    z3 = z2*z1
    #omf = omm*z3/( omm*z3 + (1.-omm-oml)*z2 + oml )
    # omf = omm*(1.+z)^3/( omm*(1.+z)^3+(1.-omm-oml)*(1.+z)^2+oml )
    iomf = ( omm*z3 + (1.-omm-oml)*z2 + oml )/(omm*z3)

    pi2 = math.pi*math.pi
    done = False
    if (abs(omm-1.) < 1.e-3) and (abs(oml) < 1.e-3):
        dc  = 3.*(12.*math.pi)**(2./3.)/20.
        ddc = 18.*pi2
        dc  = repeat(dc, z.shape)
        ddc = repeat(ddc, z.shape)
        done= True
        method='1'
    elif (abs(cosmo["OMEGA"]-1.) < 1.e-3) and (oml > 0.):
        dc  = 3.*(12.*math.pi)**(2./3.)/20.*(1. - .0123*np.log10(iomf))
        wf  = iomf-1.
        ddc = 18.*pi2*(1.+.4093*wf**.9052)
        done= True
        method='2'
    elif (cosmo["OMEGA"] < 1.):
        xx   = 2.*iomf-1.
        etaf = np.log(xx+np.sqrt(xx*xx-1.))
        sinh_etaf = np.sinh(etaf)
        cosh_etaf = np.cosh(etaf)
        dc   = 1.5*( 3.*sinh_etaf*(sinh_etaf-etaf)/(cosh_etaf-1.)**2 -2. )
        dc   *= ( 1.+(2.*math.pi/(sinh_etaf-etaf))**(2./3.) )
        ddc  = 4.*pi2*(cosh_etaf-1.)**3/(sinh_etaf-etaf)**2
        done = True
        method='3'
    if not(done):
        print 'problem'
    #else:
        #print 'method=', method
       
    ddc_fof = 3./2./fof_length**3
    ddc_fof = np.repeat(ddc_fof, z.size)

    deltac = {'Z':z, 'DC':dc, 'DDC':ddc, 'DDC_FOF':ddc_fof}
    return deltac

################
def enclosed_mass(r):
    """
    Provides the enclosed Mass at radius r in the NFW profile

    from  M(r) = 4*!pi*r^3*rho_s * f(rs/r) 
    """

    r = np.asarray(r)
    mass = r*r*r*( np.log(1.+1./r) - 1./(1.+r) )
    return mass

def convert_nm(cosmo,evol,nm,overdensity=200,rhoc=False):
    """
    Convert nm from M200b to some other mass convention assuming
    a NFW profile with using the method from Hu & Kravtsov (2003)
    """

    #; parameters from Hu & Kravtsov fit to x(f)
    a1=.5116
    a2=-.4283
    a3=-3.13e-3
    a4=-3.52e-5
    #  define mass convention
    ratio=overdensity/200.
    if(rhoc):
        linear_interp = scipy.interpolate.interp1d(np.array(evol["Z"]), np.array(evol["OMEGA_M_A"]))
        #linear_extrap = extrap1d(linear_interp)
        ratio  =  ratio/linear_interp(np.array(nm["Z"]))

    newnm = nm.copy()
    # new f value taking into account concentration c
    c200=np.double(nm["C"])
    c200i3=1./c200**3
    fnew = ratio * ( np.log(1+c200) - c200/(1+c200)) *c200i3

    # new x (=> c) value using Hu & Kravtsov's fit to x(f)
    p      = a2 + a3*np.log(fnew) + a4*(np.log(fnew)**2)
    myroot = np.sqrt(a1*fnew**(2*p)+ 9./16.) 
    ix=myroot/(1+2*fnew*myroot)
    newnm["C"] = ix

    # new M value
    newnm["M"] = nm["M"]*ratio*c200i3*(ix*ix*ix)

    # new dn/dm value
    newnm["DNDM"] = derivee(newnm["M"], nm["M"])*nm["DNDM"]
  
    # new dn/dm/domega/dz value
    newnm["DNDMDOMDZ"]= nm["DNDMDOMDZ"] /nm["DNDM"] * newnm["DNDM"] 

    return newnm
    #   {m:newnm.m,dndm:newnm.dndm,nm:newnm.nm,c:newnm.c,c_vir:newnm.c_vir,
    #    M_vir:newnm.M_vir,r_s:newnm.r_s,dndmdomdz:newnm.dndmdomdz}

#---------------------- General Function:  derivee -----------------------
def derivee(x, y):
    """
    Calculate the derivative along a single dimension.

    Calling Sequence:
        Result = derivee(x, y)

    Positional Input Arguments:
    * x  Abscissa values of y to take with respect to. 
      Numeric array of same shape and size as 
      y.  Must be monotonic and with no duplicate values.
      First positional argument out of two.

    * y  Ordinate values, to take the derivative with respect 
      to.  Numeric array vector of rank 1.  Required.  Second posi-
      tional argument.

    Output Result:
    * Derivative of y with respect to x .  
      Numeric array of same shape 
      and size as y. 

    References:
    * Press, W. H., et al. (1992):  Numerical Recipes in Fortran 
      77:  The Art of Scientific Computing.  New York, NY:  Cambridge
      University Press, pp. 180-184.

    * Wang, Y. (1999):  "Numerical Differentiation," Introduction to 
      MHD Numerical Simulation in Space, ESS265: Instrumentation, 
      Data Processing and Data Analysis in Space Physics (UCLA).
      URL:  http://www-ssc.igpp.ucla.edu/personnel/russell/ESS265/
      Ch10/ylwang/node21.html.

    Example:
    >>>x = np.arange(8)/(2.*np.pi)
    >>>y = np.sin(x)
    >>>dydx = deriv(x, y)
    >>>['%.7g' % dydx[i] for i in range(5)]
    ['0.9957836', '0.9831985', '1e+20', '0.8844179', '1e+20']
    >>>true = np.cos(x)       #- Compare with exact solution
    >>>['%.7g' % true[i] for i in range(5)]  
    ['1', '0.9873616', '0.9497657', '0.8881628', '0.8041098']
    """

    #import numpy as np

    # check of the arguments not done
    dx1=x[:-2]
    dx2=x[2:]

    dx=np.ones(x.size)
    dx[1:-1]=dx2-dx1

    dy1=y[:-2]
    dy2=y[2:]
    dy=np.zeros(y.size)
    dy[1:-1]=dy2-dy1
    derivee=dy/dx

    # edges
    derivee[0]=derivee[1]-(derivee[2]-derivee[1])/(x[2]-x[1])*(x[1]-x[0])

    derivee[-1]=derivee[-2]+(derivee[-2]-derivee[-3])/(x[-2]-x[-3])*(x[-1]-x[-2])

    return derivee

#---------------------- General Function: xgen -----------------------
def xgen(x1, x2, npoints=100, logplot=False):
    """
    Calculate an  x vector array for plotting.
    """

    # check of the arguments not done

    if (logplot):
         xx = np.logspace(math.log10(x1), math.log10(x2), num=npoints, endpoint=True)
    else:
         xx = np.linspace(x1,x2, num=npoints, endpoint=True)

    return xx

#---------------------- Dedicated Function:  pk0_raw for sigma8_norm -----------------------
def pk0_raw(cosmo, k, fit_tk=0, cmpnorm=False):
  """
  Computes something ??? 
  
  INPUTS: 
  cosmo: cosmological parameter structure
  k: comoving wave number in h-corrected physical coordinates i.e. k=k_comoving/R_0/h   [h Mpc^-1]

  Keyed Input:
      fit_tk:
       	0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 

       Call the function tk
  """

  if (not (cmpnorm) and (cosmo.has_key('PK_FEATURE'))):
         input = [-2, -1, np.arange(cosmo['PK_FEATURE']+2)]
         kf = 0.657/1.47**(input)
         F_kf = [1., 1., cosmo['PK_FEATURE'], 1., 1.]
	 # F = exp(interpol(alog(F_kf),alog(kf),alog(k))) 
         # IDL print, interpol(x, y, 1.5) 
         # python linear_interp = scipy.interpolate.interp1d(x, y) 
         linear_interp = scipy.interpolate.interp1d( log(kf), log(F_kf) )
         F = exp(linear_interp(log(k)))
  else:
         F = 1.0

  F =  F * k**cosmo['N'] * tk(k, cosmo, fit_tk)**2

  return F

#---------------------- Dedicated Function:  sigma8_norm for del2_lin -----------------------
def sigma8_norm(cosmo, fit_tk=0):
    """
    Computes sigma ??? 

    INPUTS: 
    cosmo: cosmological parameter structure

    Keyed Input:
    fit_tk:
       	0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 
  
       Call the functions xgen and  pk0_raw
    """

    #wavenumber range
    k_ran = [1.e-4, 100.]
    n_k = 1000 
    k = xgen(k_ran[0], k_ran[1], npoints=n_k, logplot=True)
   
    # compute power spectrum
    del2 = k**3 / 2. / math.pi**2 * pk0_raw(cosmo, k, fit_tk=fit_tk, cmpnorm=True)

    #  compute sigma8  do not use math.sin or math.cos
    w8 = 3./(k*8.)**3*( np.sin(k*8.) - np.cos(k*8.)*(k*8.) ) # for 8 h^-1 Mpc
    #sig8=sqrt( int_tabulated(alog(k), del2 *w8**2) )
    # np.trapz(y,x) and int_tabulated(x,y)
    sig8 = sqrt( np.trapz(del2 *w8**2, np.log(k)) )

    return sig8

#---------------------- Dedicated Function: del2_lin for ??  -----------------------
def del2_lin(cosmo, k, z=0., gtype='ODE', fit_tk=0):
    """
    compute the linear normalised density variance spectrum
    Delta^2(k)=k^3*P(k)/(2pi^2) for a given cosmological model.
    This quantity gives the contribution to the variance per log k.
    This is computed using the analytical fits quoted in Peacock and Dodds  (1996, mnras, 280, L19),
    as summarized by Peacock (1996, mnras, 284, 885). 
    The corrections for the quintessence model are derived from MA et al. (1999, APJL, 521, 1)

    INPUTS: 
    cosmo: cosmological parameter structure
    k: comoving wave number in h-corrected physical coordinates, i.e. k=k_comoving/R_0/h   [h Mpc^-1]

    OPTIONAL INPUT: 
    z: redshift (z=0 by default)
    gtype: growth type, used by growth_a, legal values:
        'ODE' for growth_a_ode   resolving differential equation
        'INT' for d=growth_a_int numerical integration of the Heath equation     
        'ANA' for growth_a_ana   analytical formulae
    fit_tk:
       	0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 

    Call the functions xgen and  pk0_raw
    """
 
    # get rid of dimension 2
    if (z.ndim == 2):
        z1d=np.double(z[:,0])
    else:
        z1d=np.double(z)

    if (k.ndim == 2):
        k1d=np.double(k[0,:])
    else:
        k1d=np.double(k)

    if (cosmo['NORM'] == 1):
          norma = (cosmo['SIGMA8']/sigma8_norm(cosmo, fit_tk=fit_tk))**2.   # sigma8 normalisation
    else:
          kk = kk=cosmo['KPIV']/cosmo['H']
          
          dnorma = 7.7374384e-14 * kk**3 / 2. / math.pi**2 * pk0_raw(cosmo, kk, fit_tk=fit_tk, cmpnorm=True)
          dnorma = dnorma*(cosmo['OMEGA_M'] / tk(kk, cosmo, fit_tk=fit_tk)/growth_a_ode(cosmo, 1., unnorm=True)/(kk*kk))
          norma = cosmo['DELTAR2']/ (dnorma*dnorma)

    # normalised power spectrum P(k) at z=0
    delta2_z0 = norma * k1d**3 / 2. / math.pi**2 * pk0_raw(cosmo, k1d, fit_tk=fit_tk)
    #print "norma", norma

    # compute the growth factor 
    # if (n_elements(z) eq 0) then d=1. else d=growth_a(cosmo,1./(1.+z),type=type)
    # the pylab function size() works for both built-in type float and numpy.array

    #if not (iterable(z)):   or np.isscalar(z)
    #    z=numpy.array(z)    

    #print 'taille ', size(z)
    if ((size(z)==1) and (z==0.)):
        d1d = 1.
    else:
        d1d = growth_a(cosmo, 1./(1.+z1d), type=gtype)

    # compute P(k,z)
    #delta2 = delta2_z0*(d*d)
    delta2 = np.outer((d1d*d1d),delta2_z0)

    return delta2

#---------------------- Dedicated Function: mk_evol  -----------------------
def mk_evol(cosmo, ran_z=[0.0,5.0], n_z=200, z=None,  gtype='ODE'):
    """
    Make an evolution table consisting of the value of various
    parameters (a,chi,g) as a function of redshift z
    Using interpol.pro this table can be used to compute a parameter as a function of another
    For eg., interpol(evol.chi, evol.z, 0.1) returns the comoving distance chi(z = 0.1)

    INPUTS: 
    cosmo: cosmological parameter structure
    k: comoving wave number in h-corrected physical coordinates, i.e. k=k_comoving/R_0/h [h Mpc^-1]

    OPTIONAL INPUT: 
    z: redshift (default: None)
    ran_z: range of z values (default: [0,5])
    n_z: number of z values (default: 200)

    gtype: growth type, used by growth_a, legal values:
        'ODE' for growth_a_ode   resolving differential equation
        'INT' for d=growth_a_int numerical integration of the Heath equation     
        'ANA' for growth_a_ana   analytical formulae
  
    OUTPUT: evol: evolution structure containing
           z: redshift
           a: expansion parameter (a=1/(1+z))
           hc: hubble 'constant' H as a function of a  [km/s/Mpc]
           omega_a,omega_m_a,_l_a,k_a: Omega_*(a)   [1]
           chi: comoving radial distance in units of R_0       [1] 
           sk: sk(chi) comoving angular radius in units of R_0 [1]
           da,dl: angular-diameter and luminosity distances    [Mpc]
           d: linear growth factor D normalised so that D(z=0)=1    [1]
           g,g0: unnormalised growth factor as defined in Ma 1998,
                    defined so that D=a*g/g_0)       [1]
           w: dark energy equation of state parameter w(a) [1]
           t: look back time [Gyr]
           th0: t*Ho look back time in unit of Ho [1]

    Call the functions xgen, etc...
    """

    if z is None: z = xgen(ran_z[0], ran_z[1], npoints=n_z)

    # expansion parameter a=R/Ro=1/(1+z)   [1]
    a = 1.0 / (1.0 + z)

    # hubble constant and dark energy quantities
    hubble = hubble_a(cosmo, a)  
    # structure contains:
    #      hc: hubble constant H [km/s/Mpc]
    #      qevol: dark energy density evolution rho_l(a)/rho_l(0) [1]
    #      w_a: dark energy equation of state parameter w(a)   [1]
    w_a = hubble['W_A']

    # evolution of the density parameters
    h_h0 = hubble['HC']/cosmo['H0']    # H/H_0
    omega_r_a = cosmo['OMEGA_R']*a**(-4)*h_h0**(-2)
    omega_m_a = cosmo['OMEGA_M']*a**(-3)*h_h0**(-2)
    omega_l_a = cosmo['OMEGA_L']*hubble['QEVOL']*h_h0**(-2)

    omega_a = omega_m_a + omega_l_a + omega_r_a
    omega_k_a = 1.0 - omega_a

    # comoving distance in units of R_0  [1]
    chi = chi_a(cosmo, a)

    # look back time
    
    th0 = t_a(cosmo, a)
    # t*Ho [1]
    
    age = t_a(cosmo,1.e-10)/cosmo['H0']*3.085678e19/3.15581498e7/1e9
    # age of the universe [Gyr]
    
    t = th0/cosmo['H0']*3.085678e19/3.15581498e7/1e9
    # t [Gyr]

    # sk(chi): radial comoving distance in units of R_0  [1]
    # np.sin() 13 March 2015
    if cosmo['OMEGA'] < 1. :
        sk = np.sinh(chi)
    elif cosmo['OMEGA'] == 1. :
        sk = chi
    else :
        sk = np.sin(chi)

    # angular-diameter and luninosity distance  [Mpc]
    da = cosmo['R0']*sk/(1.+z)
    dl = da*(1.+z)**2

    # compute the linear growth factor by numerical ODE solving

    d = growth_a(cosmo, a, type=gtype)
    # normalised to D(z=0)=1 

    # New (07/2011) compute volume

    c = 299792458.0
    # speed of light 

    d2Vdomdz = c/1e3/ hubble['HC'] * da**2. *(1.+z)**2.

    # store in structure evol
    evol = {'Z':z, 'A':a, 'HC': hubble['HC'], 'CHI':chi, 'SK':sk, 'DA':da, 'DL':dl, 'D':d,
      'W_A':w_a, 'OMEGA_A':omega_a, 'OMEGA_M_A':omega_m_a, 'OMEGA_R_A':omega_r_a, 'T':t, 'TH0':th0,
      'OMEGA_L_A':omega_l_a, 'OMEGA_K_A':omega_k_a,'N_Z':n_z, 'AGE':age,'Z_MAX':ran_z[1],
      'D2VDOMDZ': d2Vdomdz  }

    return evol

#---------------------- Dedicated Function: mk_sigmam  -----------------------
def mk_sigmam(cosmo, m_ran=[1e13,1e16] , n_m=50, k_ran=[.001,10.] , n_k=200, z=0.,
              plotit=False, gtype='ODE', fit_tk=0, fit_nl=0, verbose=False ):
    """  
    Compute the linear rms density contrast in a sphere of radius R 
    corresponding to a mass M at redshift z. 
    This is to be input in Press-Schechter calculation of the halo mass function.

    INPUT: 
      cosmo: structure produced by rd_cosmo.pro

    OPTIONAL INPUT:
      m_ran: mass M range [h^-1 M_sun], default=[1e13,1e16]
      n_m: number of M values (default=50)
      k_ran: range of k values [h Mpc^-1], default=[.001,10.]
      n_k: number of k values  (default=100)
      z: redshift (default=0)
      plotit: plot results if set
      fit_tk:
       	0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 

      fit_nl: fitting function for non-linear corrections
        0: Peacock and Dodds default
        1: Ma et al., 2: Smith et al. not yet implemented

     gtype: growth type, used by growth_a, legal values:
        'ODE' for growth_a_ode   resolving differential equation
        'INT' for d=growth_a_int numerical integration of the Heath equation     
        'ANA' for growth_a_ana   analytical formulae


    OUTPUT: sigmam structure:
      m: M [h^-1 M_sun]
      r: R(M) [h^-1 Mpc]
      sigma: sigma(M)    [1]
      dlnsdlnm: |d(ln sigma)/d(ln M)|  [1]
    """

    if (verbose): print 'mk_sigmam '+gtype

    rhoc=2.77536627e11
    # critical density [h^2 M_sun Mpc^-3]

    n_z = np.size(z)

    # compute linear power spectrum at z (default z=0.)
    # mk_del2,cosmo,z,del2,kl_ran=k_ran,n_k=n_k,gtype=gtype,fit_tk=fit_tk
    # fit_nl missing
    del2 = mk_del2(cosmo, z, kl_ran=k_ran, n_k=n_k, 
         gtype=gtype, fit_tk=fit_tk, fit_nl=fit_nl)
    k = del2['K_L']

    # compute M and R(M)
    m = xgen(m_ran[0],m_ran[1], npoints=n_m, logplot=True)
    r = (3.0 / (4.0 * math.pi) / rhoc / cosmo['OMEGA_M'] * m)**(1.0 / 3.0)
    # [h^-1 Mpc]
    
    # compute sigma(M) for each M
    mysigma = np.zeros([n_z, n_m])
    x=np.log(k)
    k3=k*k*k

    #print r[0].shape, k.shape, k3.shape, del2['PK_L'].shape
    #r[0] = scalar, k and k3 = array 200 = n_m-1,    del2['PK_L'].shape=  (64, 200) = (n_z, n_m-1)
    #pdb.set_trace()

    for i in range(n_m):
        k_t_r = k*r[i]
        # k_t_r shape 200 w same shape

        w = 3.0 / (k_t_r)**2 * (np.sin(k_t_r) / (k_t_r) - np.cos(k_t_r))

        for j in range(n_z):
            y=k3*del2['PK_L'][j,:]*w**2
            mysigma[j,i] = np.trapz(y,x)

    mysigma = np.sqrt(1./(2*math.pi**2) * mysigma)

    # compute d(ln sigma)/d(ln M) numerically
    dlnsdlnm =  np.zeros([n_z, n_m])
    for j in range(n_z):
        dlnsdlnm[j,:] = derivee(np.log(m),np.log(mysigma[j,:]))

    sigmam = {'SIGMA':mysigma, 'R':r, 'Z':z, 'M':m, 'DLNSDLNM':np.abs( dlnsdlnm ) }

    return sigmam

#---------------------- Dedicated Function: mk_del2  -----------------------
def mk_del2(cosmo, z, kl_ran=[.01,10.], n_k=100, fit_nl=0, gtype='ODE', fit_tk=0, verbose=False):
    """  
    Compute the nonlinear (and linear) normalised density variance
    spectrum Delta^2(k)=k^3*P(k)/(2pi^2) for a given cosmological model.
    This quantity gives the contribution to the variance per log k.
    This is computed from the non-linear variance spectrum according to the
    prescrition of Peacock and Dodds (1996, mnras, 280, L19),
    as summarized by Peacock (1997, mnras, 284, 885). 
    Optionally, the Ma et al. (1999, ApJ, 521, L1) or Smith et al. 
    (2002, astro-ph/0207664) formula is used.

    INPUT: 
      cosmo: cosmological parameter structure
      z redshift vector

    OPTIONAL INPUT:
      kl_ran: range of linear k values to consider [h Mpc^-1]
        (h-corrected physical units; default=[.01,100.] h Mpc^-1). 
        Note that the non-linear k range will be shifted from that.
      n_k: number of k values to consider (default=100)
      fit_nl: fitting function for non-linear corrections
        0: Peacock and Dodds default
        1: Ma et al., 2: Smith et al. not yet implemented
      fit_tk:
        0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 

      gtype: growth type, used by growth_a, legal values:
        'ODE' for growth_a_ode   resolving differential equation
        'INT' for d=growth_a_int numerical integration of the Heath equation     
        'ANA' for growth_a_ana   analytical formulae

    OUTPUT: 
      del2 structure containing:
      k: (nonlinear) comoving wave number in h-corrected physical units
      i.e. k=k_comoving/R_0/h   [h Mpc^-1]
      del2: Delta^2(k) (nonlinear) [1]
      pk: Power spectrum P(k) (nonlinear) [h^-3 Mpc^3]
      k_l: comoving wave number for the linear spectrum
      in the same units [h Mpc^-1]
      del2_l: linear Delta^2(k_nl)  [1]
      pk_l: linear power spectrum [h^-3 Mpc^3]

    CALLED FUNCTIONS: xgen, del2_lin, extrap1d, growth_a

    HISTORY: R. Gastaud, 27 February 2015, correct the resize bug 
            beware resize changes the number of elemements:
            if the new array is larger than the original array, then the new array is filled 
            with repeated copies of a
            and we have to transpose
    """

    kl = xgen(kl_ran[0], kl_ran[1], npoints=n_k, logplot=True)
    # RG read IDL fits file
    #print 'kl.shape', kl.shape
    #ref_kl = mrdfits( 'kl.fits',0)
    #print 'ref_kl.shape', ref_kl.shape
    #kl = ref_kl[0,:]   
    #print 'replaced kl by idl reference \n', kl.shape

    aa = 1./(1.+z)

    #  compute the linear variance and power spectra
    if (verbose): print ' mk_del2 call del2_lin'
    del2_l = del2_lin(cosmo, kl, z, gtype=gtype, fit_tk=fit_tk)

    # RG read IDL fits file
    #print 'del2_l.shape', del2_l.shape
    #ref_del2_l = mrdfits( 'del2_l.fits',0)
    #print 'ref_del2_l.shape', ref_del2_l.shape
    #del2_l = ref_del2_l
    #print 'replaced del2_l by idl reference \n',  del2_l.shape

    pk_l = 2.*math.pi**2/(kl**3)*del2_l
    # compute non-linear power spectrum using P&D formula unless otherwise

    if (fit_nl == 0):
        # compute the effective spectral index by computing the logarihtmic derivative of the power spectrum
        # (see eq. 25 in Peacock 97)
        neff = -3 + np.log(del2_l[0,2:]/del2_l[0,0:-2])/np.log(kl[2:]/kl[0:-2])
        kl_neff = 2.*kl[1:-1] # the factor 2 is in equation 22

        #  to check
        #nneff = np.zeros(n_k-2)
        #print "neff shape", neff.shape, "n_k", n_k, " kl", kl.shape, " del2_l", del2_l.shape
        #for i in range(1, n_k-1):
        #    nneff[i-1]=-3.+np.log(del2_l[0,i+1]/del2_l[0,i-1]) / np.log(kl[i+1]/kl[i-1])
        #diff = neff-nneff
        #print " *** diff neff", diff.min(), diff.max()
    
        linear_interp = scipy.interpolate.interp1d(kl_neff, neff)
        linear_extrap = extrap1d(linear_interp)
        n = linear_extrap(kl)
        #print "n.shape=",n.shape
        #print "n=", n[0:9], n[100]
        #from pylab import *
        #plot(n)
        n3 = (1.+n/3.)
        a  = 0.482*(n3**-.947)
        b  = 0.226*(n3**-1.778)
        alpha = 3.310*(n3**-.244)
        beta  = 0.862*(n3**-0.287)
        v     = 11.55*(n3**-0.423)
 
        ## print "*** v, alpha, beta, b, a", v, alpha, beta, b, a
  
        # compute the nonlinear variance spectrum
        # formula 21 of Peacock 1996
        d = growth_a(cosmo, aa, unnorm=True, type=gtype)

        # RG read IDL fits file
        #print "d.shape", d.shape
        #ref_d = mrdfits('d.fits', 0)
        #print "ref_d.shape", ref_d.shape
        #d = ref_d[:,0]
        #print "d changed to idl, d.shape \n", d.shape
        
        g = d/aa
        x = del2_l
        numerator = (1.+b*beta*x+(a*x)**(alpha*beta))
        #pdb.set_trace()
        ## horribile visu ! 
        #gg = np.resize(g,[100,200])
        #gg = np.reshape(g,[g.shape[0],1])
        gg = np.resize(g,[kl.size, z.size]).transpose()
        deno = 1.+((a*x)**alpha*gg**3/(v*np.sqrt(x)))**beta
        del2_nl = x*(numerator/deno) **(1./beta)

        #  for debugging write a fits file with intermediate variables
        # tel = {'v': v, 'alpha': alpha, 'beta':beta, 'a':a, 'b':b, 'n':n, 'gg':gg,'x':x}
        # from fitsutils import mwrfits
        # mwrfits(tel, 'tel_py.fits')
        # print "write tel_py.fits"
        # end of debugging

        # compute the non-linear wave number
        k_nl = (1.+del2_nl)**(1./3.)*kl

    else:
        print 'mk_del2 fit_nl='+string(fit_nl)+' not yet implemented'
        raise Exception('mk_del2 fit_nl='+string(fit_nm)+' not yet implemented')
        # 
    pass

    pk_nl = 2.*math.pi**2/k_nl**3*del2_nl
    del2 = {'PK_L':pk_l, 'DEL2_L':del2_l, 'K':k_nl, 'DEL2': del2_nl, 'K_L':kl, 'PK':pk_nl, 'Z':z }
    
    return del2

#---------------------- Dedicated Function: mk_nm -----------------------
def mk_nm(cosmo, evol, z=0., m_ran=[1e8,1e16] , n_m=60, fit_nm=2, 
          biasmodel=-1, profile=True, ctype=0, fit_tk=0, fit_nl=0, gtype='ODE',
          verbose=False, tinkerdelta=None, tinkerrhoc=False, n_k = 200):
    """  compute the  halo mass function and related quantities
    INPUT: 
    cosmo: cosmological parameter structure
    evol: evolution structure produced by mk_evol.pro

    OPTIONAL INPUT: 
    z: redshift (default:z=0)
    m_ran: mass M range [h^-1 M_sun]
    n_m: number of M values (default=50)
    fit_nm: fitting function for the mass function 
        0: Press & Schechter
        1: Shev & Tormen
        2: Jenkins et al (default)
        3: Warren
        4: Tinker
    biasmodel: fitting function for the bias
        0: Mo & White
        1: Sheth, Mo & Tormen
        -1: default: not computed
    profile: default True, flag to compute the profile, use the parameter cytpe
    ctype: type of the profile, default 0
        0:Bullock et al 2001
        else : ???

    n_k: number of k values to consider (default=200)
    fit_nl: fitting function for non-linear corrections
        0: Peacock and Dodds default
        1: Ma et al., 2: Smith et al. not yet implemented
    fit_tk:
       	0: E&H without wiggles (default), 
        1:E&H with wiggles, 
        2: BBKS as summarised by Peacock & Dodds (1997) 

    gtype: growth type, used by growth_a, legal values:
        'ODE' for growth_a_ode   resolving differential equation
        'INT' for d=growth_a_int numerical integration of the Heath equation     
        'ANA' for growth_a_ana   analytical formulae

     
    verbose: flag, default False
    tinkerdelta: default interpolation on omega_m, etc...
    tinkerrhoc: flag, default False, to be set only if tinkerdelta set by the user
    n_k = 200

    OUTPUT:
    nm structure:
           rho_bar: comoving matter density [h^2 M_sun Mpc^-3] 
           z: redshift
           m: M [h^-1 M_sun]
           sigma: sigma(M) at z   [1]
           deltac: delta_c density threshold (z-independent) [1]
           ddeltac: Delta_c density contrast (still experimental) [1]
           nu: normalised density conrtrast, nu=delta_c/sigma(M)  [1]
           dndm: differential counts dn/dM  [1/(h^-1 Mpc)^3/(h^-1 M_sun)]
           nm: integrated counts n(>M) [1/(h^-1 Mpc)^3] 
           m_star: M*, the nonlinear mass scale, for which nu(M*,z)=1 [h^-1 M_sun]
           c: NFW compactness parameter for our masses (c=R_200/R_s) [1]
   	   cvir: usual NFW compactness parameter (cvir=R_v/R_s form  [1] 
           Bullock et al 2001, onlyvalid for LCDM; see below)
           Mvir: virial mass of the halo (derived from c and cvir) [h^-1 M_sun]
           r_s:  scaling radius of the NFW profile  [h^-1 Mpc]
           r_200: radius for 200 overdensity   [h^-1 Mpc]
           r_s: NFW scale radius Rs=R_200/c   [h^-1 Mpc]
           delta_bar: NFW density amplitude [1]
           b: bias parameter b(nu) (from Mo & White or Sheth & Tormen) [1]
           f: 
           dndmdomdz: sky-projected differential mass function 
                      dN/d(omega)/dz/dM  [1/(h^-1 M_sun)/srad/deltaz]

    called functions: extrap1d, mk_sigmam, cmp_deltac
    """

    n_z = np.size(z)
    #n_k = 200
    k_ran = np.array([5.e-4, 200.])
    rhoc = 2.77536627e11      #  critical density at z=0 [h^2 M_sun Mpc^-3]

    if (tinkerdelta == None) and tinkerrhoc :
        print 'tinkerdelta == None) and tinkerrhoc set to true, forbidden'
        return None

    if not np.array_equal(evol['Z'], z):
        if verbose: print ' z not equals evol.z ==> interpolation'
        linear_interp = scipy.interpolate.interp1d(evol["Z"], evol["OMEGA_M_A"])
        linear_extrap = extrap1d(linear_interp)
        rho_bar  = rhoc*linear_extrap(z)
        linear_interp = scipy.interpolate.interp1d(evol["Z"], evol["HC"])
        linear_extrap = extrap1d(linear_interp)
        rho_crit  = linear_extrap(z)
    else:
        if verbose: print ' z equals evol.z no interpolation'
        rho_bar  = rhoc* evol["OMEGA_M_A"]
        rho_crit  = evol["HC"]

    rho_crit = rhoc*(rho_crit/cosmo["H0"])**2

    # creates a modified mass range to obtain required upper mass bound
    #   when last element will be removed 
    fake_ran = [ m_ran[0] , 10**(n_m/float(n_m-1) *
                 (np.log10(m_ran[1]) - np.log10(m_ran[0])) + np.log10(m_ran[0]))]

    fake_n_m = n_m + 1

    #  compute sigma(M) 
    #sigmam = mk_sigmam(cosmo, m_ran=m_ran, n_m=n_m, k_ran=k_ran, n_k=n_k, z=z, fit_tk=fit_tk)
    sigmam = mk_sigmam(cosmo, m_ran=fake_ran,n_m=fake_n_m, k_ran=k_ran, n_k=n_k, z=z, fit_tk=fit_tk, fit_nl=fit_nl, gtype=gtype, verbose=verbose)
    m=sigmam['M']
    if verbose: print 'sigmam.sigma=', sigmam['SIGMA'].shape

    # compute delta_c without growth factor
    deltac = cmp_deltac(cosmo, z=z)
    
    # compute nu = delta_c / sigma(M)
    mydc = np.reshape(np.tile(deltac['DC'], n_m+1), [n_m+1,n_z])
    mydc = mydc.transpose()
    if verbose: print 'mydc=', mydc.shape
    #nu = mydc[:, :-1]/sigmam['SIGMA']
    nu = mydc/sigmam['SIGMA']

    # compute M*, the nonlinear mass scale, for which nu(M*,z)=1
    m_star = np.zeros(n_z)
    log10_m = np.log10(m)
    one = np.array([1.])  ### extrap1d accept only np.array as input
    if verbose: print 'nu=', nu.shape, ' m=',m.shape, ' m_start=', m_star.shape
    for i in range(n_z):
        linear_interp = scipy.interpolate.interp1d(nu[i,:],  log10_m)
        linear_extrap = extrap1d(linear_interp)
        mypower = linear_extrap(one)
        m_star[i] = 10.**mypower

    #  compute f(nu) with different simulation fits, default is 2 Jenkins 
    if fit_nm == 0:  #   Press-Schechter
        f = np.sqrt(2./math.pi)*np.exp(-nu**2/2.)
        mf_delta = np.tile( 200 , n_z) ##  F.o.F 200 instead of tinkerdelta ?
    elif fit_nm == 1: # Sheth & Tormen 1999
        # parameters from Robinson 2000 (astro-ph/0004023?) or Seljak 2000 (astro-ph/0001493)
        a_st  = .707
        q_st  = .3
        aa_st = .3222*math.sqrt(2.*a_st/math.pi)
        f = aa_st * (1. + (a_st**(-q_st))*(nu**(-2.*q_st)))  * np.exp(-a_st*nu**2/2.)
        mf_delta = np.tile( 200 , n_z) ##  F.o.F 200 instead of tinkerdelta ?
    elif fit_nm == 2: # Jenkins et al. (2001)
        norm_J  = 0.315
        shift_J = 0.61
        pow_J   = 3.8
        Ln_InvertSig = np.log(1./sigmam['SIGMA'])
        f = norm_J*np.exp(-np.abs(Ln_InvertSig+shift_J)**pow_J)/nu
        if  (Ln_InvertSig.min() <  -1.2) or ( Ln_InvertSig.max() > 1.05):
            print ' Ln_InvertSig out of bonds'
        if  Ln_InvertSig.min() <  -1.2:
            print ' Ln_InvertSig lower than -1.2 ', Ln_InvertSig.min()
        if  Ln_InvertSig.max() > 1.05:
            print ' Ln_InvertSig gerater than 1.5', Ln_InvertSig.max()
        
        mf_delta = np.tile( 200 , n_z) ##  F.o.F 200 instead of tinkerdelta ?
    elif fit_nm == 4: # Tinker
        Tink_del_arr = np.array([ 200.,  300.,  400.,  600.,  800., 1200., 1600., 2400., 3200.])
        Tink_aa_arr  = np.array([0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260])
        Tink_a_arr   = np.array([1.47,  1.52,  1.56,  1.61,  1.87,  2.13,  2.30,  2.53,  2.66])
        Tink_b_arr   = np.array([ 2.57,  2.25,  2.05,  1.87,  1.59,  1.51,  1.46,  1.44,  1.41])
        Tink_c_arr   = np.array([ 1.19,  1.27,  1.34,  1.45,  1.58,  1.80,  1.97,  2.24,  2.44])
        # ; default spherical overdensity for Tinker M.F. = M200c
        # remark no need for interpolation if evol['Z'] equals z
        if tinkerdelta == None :
             if  not np.array_equal(evol['Z'], z):
                 linear_interp = scipy.interpolate.interp1d( evol['Z'], evol['OMEGA_M_A'])
                 linear_extrap = extrap1d(linear_interp)
                 tinkerdelta = 200./linear_extrap(z)
             else:
                 tinkerdelta = 200./ evol['OMEGA_M_A']
        if  np.size(tinkerdelta) == 1: 
             mf_delta = np.tile(tinkerdelta, n_z) 
        else:
             mf_delta = tinkerdelta
        if  np.size(mf_delta) != n_z: 
             print 'Invalid mf_delta size',np.size(mf_delta), 'n_z=', n_z 
             return None
        if tinkerrhoc:
            if  not np.array_equal(evol['Z'], z):
                 linear_interp = scipy.interpolate.interp1d( evol['Z'], evol['OMEGA_M_A'])
                 linear_extrap = extrap1d(linear_interp)
                 mf_delta = mf_delta/linear_extrap(z)
            else:
                 mf_delta = mf_delta/evol['OMEGA_M_A']
 
        log10_mf_delta = np.log10(mf_delta)
        log10_Tink_del_arr = np.log10(Tink_del_arr)
        # size of Tink_del_arr is 9, size of mf_delta is n_z, here 100, need to interpolate
        # Tink_aa0 will be size n_z
        #if  not np.array_equal(Tink_del_arr, mf_delta):
            
        Tink_aa0 = rg_interpolator(log10_Tink_del_arr, Tink_aa_arr, log10_mf_delta)
        Tink_a0  = rg_interpolator(log10_Tink_del_arr, Tink_a_arr,  log10_mf_delta)
        Tink_b0  = rg_interpolator(log10_Tink_del_arr, Tink_b_arr,  log10_mf_delta)                           
        Tink_c   = rg_interpolator(log10_Tink_del_arr, Tink_c_arr,  log10_mf_delta)

        Tink_alpha = 10**( - ( 0.75/np.log10(mf_delta/75.) )**1.2)

        zplus1 = z + 1. 
        Tink_aa  = Tink_aa0 *zplus1**(-0.14)
        Tink_a   = Tink_a0  *zplus1**(-0.06)
        Tink_b   = Tink_b0  *zplus1**(-Tink_alpha)

        my_ones = np.tile(1.,n_m+1)
        faa = np.outer(Tink_aa, my_ones)
        fa  = np.outer(Tink_a,  my_ones)
        fb  = np.outer(Tink_b,  my_ones)
        fc  = np.outer(Tink_c,  my_ones)

        # rg print faa.shape, fa.shape, fb.shape, fc.shape, fa.__class__
        # rg print sigmam['SIGMA'].shape, nu.shape

        f   = faa * ( (sigmam['SIGMA']/fb)**(-fa) + 1. ) * np.exp(-fc/sigmam['SIGMA']**2.)/nu 
    else:
         print 'not yet implemented fit_nm', fit_nm

    # Mass function [dn/dm](z)
    dndm = cosmo['OMEGA_M']*rhoc/m**2*nu*sigmam['DLNSDLNM']*f

    #compute cummulative counts n(>M)
    nm = np.zeros([n_z, n_m])
    # rg print 'n_m=', n_m, '  n_z=', n_z, ' m.shape', m.shape, ' dndm.shape', dndm.shape
    for j in range(n_z):
        tampon = scipy.integrate.cumtrapz(m[::-1]*dndm[j, ::-1], np.log(m[::-1]))
        #print 'tampon.shape', tampon.shape,nm[j, 0:n_m-1].shape 
        nm[j,:] = -tampon[::-1]
   
    # compute sky projected differential number counts
    myones = np.tile(1.,n_m+1)
    #a = 1.(1. + outer(z, myones))

    #  Correction 07/2011 'D2VDOMDZ'
    dndmdomdz =  rg_interpolator(evol['Z'], evol['D2VDOMDZ'], z)
    dndmdomdz = np.outer(dndmdomdz, myones)
    if (verbose): print ' dndm shape', dndm.shape, dndmdomdz.shape
    dndmdomdz = dndmdomdz*dndm*cosmo['H']**3

    nm_str =  {'Z':z, 'RHO_BAR':rho_bar, 'RHO_CRIT':rho_crit, 'DELTAC':deltac['DC'], 'DDELTAC':deltac['DDC'] ,
               'MF_DELTA':mf_delta, 'M_STAR':m_star, 'M':m[:-1], 'SIGMA':sigmam['SIGMA'][:, :-1], 'NU':nu[:, :-1],  
               'F':f[:, :-1], 'DNDM':dndm[:, :-1], 'NM':nm, 'DNDMDOMDZ': dndmdomdz[:, :-1] }

    # concentration model
    if profile:
        if verbose: print ' concentration model ', ctype
        if ctype == 0: #  Bullock et al 2001
            # compute M*(z=0), to be used in the Bullock et al. formula for c 
            if ((size(z)==1) and (z==0.)):  # z not set
                Mstar0 = m_star[0]
            else: # z set
                # compute sigma(M) at z=0
                #print 'gtype=', gtype, 'fit_tk=', fit_tk, cosmo.keys()
                zero = np.array([0.])  # del2_lin does not work with scalar
     	        sigma0 = mk_sigmam(cosmo, z=zero, gtype=gtype, n_m=500, m_ran=[1.e10,1.e14], fit_tk=fit_tk, fit_nl=fit_nl, verbose=verbose)
    	        # compute delta_c at z=0
                delta0 = cmp_deltac(cosmo, z=zero) # does not work with scalar
    	        # compute nu at z=0
     	        nu0 = delta0['DC']/sigma0['SIGMA']
                # Finally M*(z=0)
                nu0 = nu0.flatten()
     	        Mstar0 = rg_interpolator( nu0, np.log10(sigma0['M']),  np.array([1.]))
     	        Mstar0 = 10.** Mstar0
                if (verbose): print 'nu0', nu0.shape, '  Mstar0 =',  Mstar0
 
            # compute the concentration parameter for the NFW profile. 
            # This is done by adopting the expression from Bullock 2001 (MNRAS, 321, 359)
            # Translation from Mv to M200 is obtained  using the appendix of Hu & Kravtsov 2003 (ApJ 584, 702)
            # WARNING: this is only valid for a LCDM model !!!
            # parameters of the bullock et al. relation: 
            # c(Mv)=Bnorm/(1+z)*(Mv/Mstar)^(Bexp)
            Bnorm = 9.
            Bexp = -0.13
            cv_test = (np.arange(10000.)+1.)/100.
            cv = np.zeros([n_z, n_m+1])
            c200 = np.zeros([n_z, n_m+1])
            for i in range(n_z):
                z_test = z[i]
                ddc_test = deltac['DDC'] [i]
                # invert_f_nfw replaced by inv_f_nfw
                m200_test = (Bnorm/(1+z_test))**(-1/Bexp) * Mstar0 * 200./ddc_test * inv_f_nfw( 200./ddc_test * f_nfw(1/cv_test))**(-3.) * cv_test**(-(3*Bexp-1)/Bexp)

                m200_test = m200_test.flatten()
                #print ' debug ', cv[i,:].shape, m200_test.shape, cv_test.shape,' m.shape',  m.shape
                #   m and not m[i,:] because m is one dimensionnal
                cv[i,:]   = rg_interpolator(m200_test, cv_test, m)
                c200[i,:] = inv_f_nfw(200./ddc_test*f_nfw(1/cv[i,:]))

            c200 = 1./c200
            Mv = m*(np.outer(deltac['DDC'], myones ))/200.*(cv/c200)**3.
        else: ###  ctype > 0
            print ' not yet implemented ctype=', ctype
            c200 = np.zeros([n_z, n_m])
	    cv   = np.zeros([n_z, n_m])
	    Mv   = np.zeros([n_z, n_m])
            #for i in range(n_z):
                #cnfw = cmp_cnfw(cosmo,nm[i].m,type=ctype,/m200,z=nm[i].z
	    	#c200[*,i]=cnfw.c200
	    	#cv[*,i]=cnfw.c_vir
	    	#Mv[*,i]=cnfw.m_vir

    #  add profile parameters
    #   idl 0:n_m-1  or in python [:, :-1]
    nm_str['C']     =  c200[:, 0:n_m]
    nm_str['C_VIR'] =  cv[:, 0:n_m]
    nm_str['M_VIR'] =  Mv[:, 0:n_m]              
    
    return nm_str
#---------------------- END of Dedicated Function: mk_nm -----------------------

##########
# extrap1d
# http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-a-an-extrapolated-result-beyond-the-input-rang
###########
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike

def rg_interpolator(x0, y0, x1):
    ### extrap1d accept only np.array of []  as input
    if np.size(x1) == 1 :
        x1 = np.array([x1])

    if  not np.array_equal(x0, x1):
        diff = x0[1::] - x0[0:-1]
        dmin = diff.min()
        dmax = diff.max()
        if (dmin*dmax <= 0):
            print "error in rginterpolator input x non monotonic"
            raise Exception("rg_interpolator input x non monotonic")
            return -1

        if (dmin < 0):
           #print "warning in rginterpolator input x decreasing"
           x0=x0[::-1]
           y0=y0[::-1]
 
        linear_interp = scipy.interpolate.interp1d(x0, y0)
        linear_extrap = extrap1d(linear_interp)
        y1 = linear_extrap(x1)
    else :
        y1 = y0

    return y1

#---------------------- Dedicated Function: sigma8  -----------------------
def sigma8(cosmo, kl_ran=[.0001,2000], n_k=5000, gtype='ODE', fit_tk=0):
    """  
    Compute  sigma8
    INPUT: 
          cosmo: cosmological parameter structure
       
    OPTIONAL INPUT: 
    kl_ran: range of linear k values to consider [h Mpc^-1]
    (h-corrected physical units; default=  [.01,100.] h Mpc^-1). Note that the non-linear
    k range will be shifted from that.           
    n_k: number of k values to consider (default=100)
    fit_tk: 

    OUTPUT: 
    sigma8 a scalar

    called functions: xgen, del2_lin
    """

    # wavenumber range
    k = xgen(kl_ran[0], kl_ran[1], npoints=n_k, logplot=True)
  
    #  compute power spectrum  for z=0 (now)
    z0 = z=np.array(0.)  # del2_lin does not work with scalar
    del2 = del2_lin(cosmo, k, z0, gtype=gtype, fit_tk=fit_tk)
    #compute sigma8
    k8=k*8
    w8 = 3./(k8**3)*( np.sin(k8) - np.cos(k8)*k8 ) # for 8 h^-1 Mpc
    x = np.log(k)
    y = del2 *w8**2
    sig8 = math.sqrt(np.trapz(y, x))
    #sig8=sqrt( int_tabulated(math.log(k), del2*w8**2) )

    return sig8
pass

#---------------------- Dedicated Function: f_nfw  -----------------------
def f_nfw(x):
    """ 
    f_nfw compute  
    INPUT: x      
    OPTIONAL INPUT: x.
    OUTPUT: y

    called functions: none
    """

    y = x**3 *(np.log(1.+1./x) -1./(1.+x))
    return y
#---------------------- Dedicated Function: inv_f_nfw  -----------------------
def inv_f_nfw(y):
    """
    f_nfw compute  
    INPUT: y
    OPTIONAL INPUT: x.
    OUTPUT: x

    called functions: none
    """

    # a = np.array([ 0.53842618, -0.41760202, -0.0018323937, 1.1329062e-05,
    #-1.3502223e-05, 0.0094997755, -0.0014941493,6.4596344e-05])

    a = np.array([0.5116, -0.4283, -3.13e-3, -3.52e-5])

    pp = (a[1] + a[2]*np.log(y) + a[3]*(np.log(y))**2.)*2
    x  = ( a[0]*y**pp  + (3./4.)**2. )**(-1./2.) + 2*y

    #x  = x + a[4] + a[5]*y + a[6]*y**2 + a[7]*y**3
    
    return x
