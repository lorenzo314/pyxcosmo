##########
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

  CALLED FUNCTIONS: xgen, del2_lin, growth_a

  HISTORY:  R. Gastaud, 27 February 2015, correct the resize bug 

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

