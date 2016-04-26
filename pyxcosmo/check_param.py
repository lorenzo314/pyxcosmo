###########################################################################
###       check_param
###
###       version 5 mai 2015 lt struct broken
###       version 4 RG 28 February 2015
###       add check_param_tables to version 3
###
###########################################################################
from fitsutils import mrdfits
import os.path
import numpy as np

__version__ = 5

#---------------------- Dedicated Function: check param
#---------------------- -----------------------

def check_param_tables(struct, tables, verbose=False):
    status=1
    for name in tables:
        if not(struct.has_key(name)):
            if(verbose):print ' missing ', name
            status = -1
        else:
            fstatus = os.path.isfile(struct[name])
            print 'check if file exists', fstatus, name, struct[name]
            if(name[-4:] == 'FITS') and fstatus:
                table = mrdfits(struct[name])
                struct.update({name[:-4:]:table})
                if(verbose):print ' add table ', name[:-4:]
            if(name == 'CRTABFITS_2') and fstatus:
                table = mrdfits(struct[name])
                struct.update({'CRTAB_2':table})
                if(verbose):print ' add table ', name

    return status

#---------------------- -----------------------
#---------------------- -----------------------

def check_param_default(struct, default, verbose=False):
    for key in default.keys():
        #print key
        if not(struct.has_key(key)):
            struct[key] = default[key]
            if (verbose):print 'add key', key, ' value', default[key]
        pass
    return

#---------------------- -----------------------


def check_param_xamin(fstruct, cosmo, type, verbose):
    if (verbose): print ' function check_param_xamin'


    model_default =   { 'SELTYPE': 'selfunc', 'SFUNCTAB': 0, 'CRTAB': 0}
    model_tables  =    ['CRTABFITS', 'SFUNCTABFITS'] 

    mydict = fstruct.copy() 

    check_param_default(mydict, model_default, verbose=verbose)

    status = check_param_tables(mydict, model_tables,  verbose=verbose) 

    return mydict

#---------------------- -----------------------
#---------------------- -----------------------


def check_param_model(fstruct, cosmo, type, verbose):
    if (verbose): print ' function check_param_model'
    #  not defined 'NPROC',
    #  defined after :
    #     'n_m',ceil(total(alog10(stru.m_ran)*[-1,1]) * 50.)
    #   cr_min_2, cr_max_2, n_cr_2

    model_default =  {'N_CR':100, 'N_Z': 200, 'Z_PEIGNE':1, 'N_T':300, 'TLIM':0.5,
     'N_L':300,'TIME_EXPO':1e4,'XC_Z':0,'LIBDIR':'','Z_RAN': [0.05,1.8],
     'XC_H':0., 'CR_MAX':1e1, 'M_RAN':[1.e12,1.e16],'XC_0':0.303,'CRLIM':1e-6,'PHYS_RC':-1,
     'CR_MIN':1e-3}

    mydict = fstruct.copy() 

    check_param_default(mydict, model_default, verbose=verbose)

    if not(mydict.has_key('N_M')):
        n_m = np.ceil((np.log10(mydict['M_RAN'])*[-1,1]*50.).sum())
        mydict['N_M']= np.int32(n_m)

    if not(mydict.has_key('N_CR_2')):
        mydict['N_CR_2'] = mydict['N_CR'] 

    if not(mydict.has_key('CR_MIN_2')):
        mydict['CR_MIN_2'] = mydict['CR_MIN'] 

    if not(mydict.has_key('CR_MAX_2')):
        mydict['CR_MAX_2'] = mydict['CR_MAX'] 

    if (mydict.has_key('MLIMZ_M')):
        print 'Adding a M(z) function...'
        if not(mydict.has_key('MLIMZ_Z')):
            print ' miss a vector MLIMZ_Z'
            #mydict['CR_MAX_2'] = mydict['MLIMZ_M'] 

    return mydict

#---------------------- -----------------------
#---------------------- -----------------------
#---------------------- -----------------------


def check_param_fluxmes(fstruct, cosmo, type, verbose):
    if (verbose): print ' function check_param_fluxmes'
    #
    compulsory = ['CRTABFITS', 'CRTABFITS_2', 'INDEXCRTABFITS', 'OBSTYPE', 'PROFTABLE']
    for key in compulsory:                                                    
        if not(fstruct.has_key(key)): raise Exception("check_parama_fluxmes"+key)

    mydict = fstruct.copy() 
    filenames = ['CRTABFITS', 'CRTABFITS_2', 'INDEXCRTABFITS']
    check_param_tables(mydict, filenames, verbose=False)

    # default value
    if not(mydict.has_key('MES_ERROR')):
        mydict['MES_ERROR']=0.1

    return mydict

#---------------------- -----------------------
#---------------------- -----------------------


def check_param_lt(lt_struct, cosmo, type, verbose):
    if (verbose): print ' function check_param_lt'
    #  5 mai 2015
    if not(lt_struct.has_key('POW')): 
        lt_struct['POW'] = 2.88
    lpower = lt_struct['POW']
    #
    if not(lt_struct.has_key('LSTAR')): 
        raise Exception("check_parama_lt LSTAR make_evol")

    if lt_struct.has_key('AUTHOR'): 
        if (lt_struct['AUTHOR'] == 'markevitch'): 
            lpower = 2.64
            raise Exception("check_parama_lt LSTAR make_evol")
            #mk_evol,cosmo,ev1,z=0.05
        if (lt_struct['AUTHOR'] == 'pratt2009'): lpower = 2.70
   
    #
    lt_default =  {'POW':lpower, 'BROKEN1':lpower,'BROKEN2':lpower, 'BROKEN3':lpower,\
          'BROKEN4':lpower, 'HEVOL':1., 'ZEVOL':0., 'TEVOL':0., 'SCATTER_EVOL':0., 'SCATTER':0.6}

    mydict = lt_struct.copy() 

    check_param_default(mydict, lt_default, verbose=verbose)

    return mydict

#---------------------- -----------------------

def check_param_mt(mt_struct, cosmo, type, verbose):
    if (verbose): print ' function check_param_mt'
    #if (verbose): print 'struct', struct
    hfac = 0.7
    m200 = 5.74e14*(4./5.)**1.49*hfac
    m500 = 4.14e14*(4./5.)**1.49*hfac
    
    mt_default={'MTYPE':'m200', 'POW':1.49, 'HEVOL':-1.0, 'ZEVOL':0, 'SCATTER':-999.}
    #dict.update(struct)  5 May 2015
    mydict = mt_struct.copy()
    check_param_default(mydict, mt_default, verbose=verbose)

    if (mydict['MTYPE'] == 'm500'):
        Mstar = m500
        mtype = 'm500'
    else:
        Mstar = m200
        mtype = 'm200'
    #dict['MSTAR']   = Mstar
    #dict = dict.fromkeys('MSTAR', Mstar)
    power = mydict['POW']

    if (mydict.has_key('AUTHOR')):
        author = mydict['AUTHOR']
    else:
        author = 'none'

    done = False
    if (author == 'all500'):
        # Uses Arnaud et al. M500-T
        mtype = 'm500'
        power = 1.71
        Mstar =  3.84e14*(4./5.)**(mydict['POW'])*hfac
        done = True
        
    if (author == 'all200'):
        # Uses Arnaud et al. M200-T
        mtype = 'm200'
        power = 1.72
        Mstar =  5.34e14*(4./5.)**(mydict['POW'])*hfac
        done = True


    if (author == 'vikhlinin'):
        # Uses  Vikhlinin et al 2006 M500-T
        mtype = 'm500'
        power = 1.58
        Mstar =  1.26e14*(4./5.)**(mydict['POW'])*hfac
        done = True

    if (author == 'sun'):
        # Uses  Sun et al. 2009 M500-T
        mtype = 'm500'
        power = 1.65
        Mstar =  2.89e14*(4./5.)**(mydict['POW'])*hfac
        done = True
   
    #if (done):
    mydict['MTYPE']   = mtype
    mydict['MSTAR']   = Mstar
    mydict['BROKEN1'] = power
    mydict['BROKEN2'] = power
    mydict['BROKEN3'] = power
    mydict['BROKEN4'] = power
    mydict['POW']     = power

    return mydict

#-------------------------------------------------------------------------

def check_param_ml(ml_struct, cosmo, type, verbose):
    if (verbose): print ' function check_param_ml'

    hfac = 0.7
    m200 = 5.74e14*(4./5.)**1.49*hfac
    m500 = 4.14e14*(4./5.)**1.49*hfac
    
    ml_default={'MTYPE':'m200', 'POW':0.52, 'HEVOL':-1.5, 'ZEVOL':0, 'SCATTER':-999.}
    mydict = ml_struct.copy()
    check_param_default(mydict, ml_default, verbose=verbose)

    #if (mydict['MTYPE'] == 'm500'):
    #    Mstar = m500
    #    mtype = 'm500'
    #else:
    Mstar = m200
    mtype = 'm200'

    power = mydict['POW']

    mydict['MTYPE']   = mtype
    mydict['MSTAR']   = Mstar
    mydict['BROKEN1'] = power
    mydict['BROKEN2'] = power
    mydict['BROKEN3'] = power
    mydict['BROKEN4'] = power
    mydict['POW']     = power

    return mydict

def check_param(struct, cosmo, type, verbose=False):
    mytypes = {'mt': check_param_mt, 'lt': check_param_lt , 'fluxmes': check_param_fluxmes, 'model':check_param_model, 'xamin': check_param_xamin }
    # mytypes[type]()
    if mytypes.has_key(type):
        new_struct = mytypes[type](struct, cosmo, type, verbose)
    else:
        print 'error check_param type:', type, ' unknown'
        new_struct = struct

    return new_struct
