
import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
import os.path
import sys
import pdb
import asciitable  # to be replaced by astropy.io.ascii
import os  # to decompose the absolute file name
import imp  # Access the import internals  useful to get the path of the sources

from fitsutils import mrdfits
from fitsutils import mwrfits

from get_cosmo import get_cosmo

#mk_nm_array, init_cmp

#from make_selfunc import make_selfunc
from make_selection_function import make_selection_function

from make_detcrtab import make_detcrtab  # brouillon 12
from make_mesindexcrtab  import make_mesindexcrtab

# make_meshrtab  ?? brouillon12
from cmp_dndcrdhr import cmp_dndcrdhr  # brouillon17
from combine_dndzdldt import combine_dndzdldt  # /brouillon14
from init_cmp_dndcr import init_cmp_dndcr # brouillon26

from read_ptarr_dndzdldt import read_ptarr_dndzdldt  ## use for debugging

__version__ = 4

def read_table_rg(filename):
    # astropy.io.ascii.convert_numpy
    converters = { 'col1': [asciitable.convert_list(str)],
    'col2': [asciitable.convert_list(float)],
    'col3': [asciitable.convert_list(float)],
    'col4': [asciitable.convert_list(float)],
    'col5': [asciitable.convert_list(float)],
    'col6': [asciitable.convert_list(str)],
    'col7': [asciitable.convert_list(str)],
    'col8': [asciitable.convert_list(str)],
    'col9': [asciitable.convert_list(str)],
    'col10': [asciitable.convert_list(float)]}


    names = ('POINT_NAME', 'RA_PNT',  'DEC_PNT', 'BACK','NH', 'MOS1FILTER', 'MOS2FILTER', 'PNFILTER', 'EXPOTIME')

    data = asciitable.read(filename, delimiter='|', converters=converters)
    npoint = data.size

    name = data['col1'].tolist()

    back = data['col4']
    nh   = data['col5']

    m1filt = data['col6'].tolist()
    m2filt = data['col7'].tolist()
    pnfilt = data['col8'].tolist()

    timestamp = data['col9']
    area = data['col10']

    return (name, back, nh, m1filt, m2filt, pnfilt, timestamp, area, npoint)

#####################################################################
####### beginning of the function calculate_xclass_dndcrdhr  ########

def calculate_xclass_dndcrdhr(rootdir, param, table, dirdir, outfile, 
     mesband='05-2', crspec=[.001,1.,64], meshr='1-2R05-1', hrspec=[.25,4,64],
     aperture=[30.,60.,120.], save_all=False, nproc=1, z_ran=[0.05,1.8],
     n_z=64, phys_rc=125, n_t=150, n_l=180, libdir=None, verbose=False, debug=None,
     proftable='Profils/profile_table_betamodel_10.0000_2000.00_80_0.0100000_0.900000_40.fits',
     in_selfunc = 'selec_func/selfunc_backevol/criteres_50-5_40-5/SaveProb_EXT_C2-0-600_summary_CLUSPNTback_',
        gtype='ODE'):

    """
    Compute an observable XCLASS (dn/dcr/dHR)
    INPUT : 
         rootdir :  A directory where are the input data (count rate tables),
                   the profile table, the selection function table
         param : A structure containing : 
                      * the cosmology
                      * the scaling relations
                      * NO SURVEY SIZE (added for each pointing)
                      * NO SELECTION function (added for each pointing)
         table : A standard [updated] XCLASS pointing table (with NH,
                 background, area...)
         dirdir : a directory where we store the results and
               intermediate files
         outfile : the name of the output file
         save_all : allows to save intermediate dN/dCR/dHR for each
                    pointing from the list
         libdir: does not exist in idl version
         proftable: does not exist in idl version
         in_selfunc: does not exist in idl version
         debug: does not exist in idl version
         gtype : does not exist in idl version, type of growth_a, combine_dndzdldt, init_cmp_dndcr
         aperture: size of the aperture, scalar, list or numpy array
         
    OUTPUT : this function returns a list containing 
            1) the array of all readshifts
            2) an index for select some redshfits
            3) the observables on the selected redshifts



	HISTORY:
        RG 13 October 2015  version 4 : use sys.modules to get the directory libdir, 
            aperture is tranformed in numpy array (so we can use size)
        RG 09 October 2015 version 3 : add gtype, see init_cmp_dndcr, combine_dndzdldt
        RG 12 May 2015 version 2 : create dirdir if it does not exist.
        RG 11 May 2015 version 1
        

    """

    aperture = np.array(aperture)
    if not(rootdir.endswith(os.path.sep)): rootdir=rootdir+os.path.sep
    if not(dirdir.endswith(os.path.sep)): dirdir=dirdir+os.path.sep
    if not os.path.exists(dirdir):
        os.makedirs(dirdir)
        print 'create the directory '+dirdir
    if (libdir == None):
        #toto =imp.find_module('get_cosmo')
        #head, tail = os.path.split(toto[1])
        mypath = (sys.modules['pyxcosmo.get_cosmo']).__file__
        head, tail = os.path.split(mypath)
        libdir = head + os.path.sep
    #
    print 'libdir', libdir
    print 'proftable', proftable
    print 'rootdir=',rootdir
    print 'rootdir+table', rootdir+table
    print 'param:', param, 
    print ''
    print 'dirdir=', dirdir, 'outfile=', outfile
    print nproc, z_ran, n_z, phys_rc, n_t, n_l
    print
    #/Users/gastaud/icosmo/idldirs_26012012/icosmo_v0.1/blah.test

    ##-------------------------------------------
    ## I/ Read the pointing table 
    ##-------------------------------------------
    if (verbose): print '  I/ Reading the XCLASS pointing table'
    (name, back, nh, m1filt, m2filt, pnfilt, timestamp, area, npoint) = read_table_rg(rootdir+table)

    timefloat = float(timestamp[0][:-2]) ### remove the ks at the end

    if (verbose): 
       	print "Simulate XCLASS survey :"
       	print "  Npoint = ", npoint
       	print "  Median Nh = ", np.median(nh)
       	print "  Median back = ", np.median(back)
       	print "  timefloat = ", timefloat


    ##-------------------------------------------
    ## II/ Initialize the pointing-independent quantities (dndzdldt)
    ##-------------------------------------------
    if (verbose): print '  II/ Setting up pointing-independent quantities'

    n_cr = crspec[2]
    n_hr = hrspec[2]

    n_ap = aperture.size
    xclassdndcrdhr = np.zeros([n_ap, n_hr, n_cr])
    #xclassdndcrdhr = np.zeros([n_cr, n_hr, n_ap])
    # xclassdndcrdhr = dblarr(n_cr,n_hr,n_ap)
    #print 'xclassdndcrdhr.shape', xclassdndcrdhr.shape
    if (verbose): print 'n_cr=', n_cr, 'n_hr=', n_hr, 'n_ap=', n_ap
    #get the cosmo
    ol = 1-param['OMEGA_M']
    oc = param['OMEGA_M'] - param['OMEGA_B']

    (cosmo, evol) = get_cosmo(wmap5=True, omega_b=param['OMEGA_B'],  \
      omega_c=oc, omega_l=ol, sigma8=param['SIGMA8'],  \
    w_l=param['W0'], w1_l=param['WA'], h=param['HUBBLE'],  \
    n_pk=param['N_PK'], tau=param['TAU'])
    #################################
    mt_struct = {'POW': param['POW_MT'], 'MSTAR': 1e14*10**param['NORM_MT'],  \
               'HEVOL': param['HEVOL_MT'],'ZEVOL': param['ZEVOL_MT'], \
                     'SCATTER': param['SCATT_MT'], 'STRUCT_NAME':'mt'}

    lt_struct = {'POW': param['POW_LT'], 'LSTAR': 1e44*10**param['NORM_LT'],  \
                  'HEVOL': param['HEVOL_LT'],'ZEVOL': param['ZEVOL_LT'], \
        'SCATTER': param['SCATT_LT'], 'SCATTER_EVOL': param['SCATTER_LT_EVOL'], 'STRUCT_NAME':'lt'}

    full_proftable = rootdir+proftable

    fluxmes_struct = {'CRTABFITS':"", 'CRTABFITS_2':"", 'INDEXCRTABFITS':"", 'OBSTYPE':"crhr",
                      'APERTURE':aperture, 'PROFTABLE': full_proftable, 'STRUCT_NAME':'fluxmes'}

   
    z_peigne = 2
    time_expo = 10000.
    nproc = 1
    #pdb.set_trace()

    model_struct = { 'PHYS_RC':phys_rc, 'N_Z':n_z, 'Z_PEIGNE':z_peigne, 'N_T':n_t, 
                     'N_L':n_l, 'TIME_EXPO':time_expo, 'CR_MIN':crspec[0], 'CR_MAX':crspec[1],
                     'N_CR':n_cr,  'CR_MIN_2':hrspec[0], 'CR_MAX_2':hrspec[1], 'N_CR_2':n_hr,
                     'XC_0':param['XC_0'], 'XC_H':param['XC_H'], 'XC_Z':param['XC_Z'],
                     'NPROC':nproc, 'LIBDIR' : libdir, 'Z_RAN':z_ran, 'STRUCT_NAME':'model'} 


    #model_struct = {'N_CR':n_cr, 'N_Z':n_z, 'Z_PEIGNE':z_peigne, 'N_T':n_t,  'TLIM':tlim,
    #                'N_L':n_l, 'TIME_EXPO':time_expo, 'XC_Z': xc_z,
    #                'LIBDIR' : libdir, 'Z_RAN':z_ran, 'XC_H':xc_h,'NPROC':nproc,
    #                'CR_MAX':crmax, 'CR_MAX_2':crmax_2, 'M_RAN':m_ran, 'XC_0':xc_0,
    #                'CRLIM':crlim, 'PHYS_RC':phys_rc, 'N_CR_2':n_cr_2, 'CR_MIN_2':cr_min_2,
    #                'CR_MIN':cr_min}

    #pdb.set_trace()  in idl data_list is replaced by zvector which is an array of pointers
    #data_list = init_cmp_dndcr(cosmo, MT_struct, LT_struct, FLUXMES_struct, MODEL_struct, verbose=False, directory=None):
    
    
    #pdb.set_trace()
    if(debug):
        data_list = read_ptarr_dndzdldt('init_cmp_dndcr_out')
        print '*********  debug read data_list init_cmp_dndcr_out  ******************'
    else:
        data_list = init_cmp_dndcr(cosmo, mt_struct, lt_struct, fluxmes_struct, model_struct, gtype=gtype)

    # idl  combine_dndzdldt,directory,ptarr,cosmo,MODEL_struct,FLUXMES_struct,verbose=verbose,clean=clean
    #combined = combine_dndzdldt(dir, zvector, cosmo, Model_struct, FLUXMES_struct, verbose=verbose)
    
    if(debug):
        combined = mrdfits('combined_idl.fits')
        print '*********  debug read combined  ******************'
    else:
        combined = combine_dndzdldt (data_list, cosmo, model_struct, fluxmes_struct, gtype=gtype, verbose=verbose)


    ##-------------------------------------------
    ## III/ Process each pointing independently
    ##-------------------------------------------
    if (verbose): print '  III/ Simulating each pointing separately'

    for i in np.arange(npoint):
        if (verbose): print i, '### Pointing : '+name[i]+' ###'

        if(area[i] <  1e-6):
            print 'area too small', area[i]  
            continue # next i

        ###  Construct the good selection function
        outselfunc = dirdir+'/'+table+'_'+name[i]+'_xclass_selfunc.fits'

        if (os.path.isfile(outselfunc)):
            print 'selection table ', outselfunc, ' already exists'
        else:
            # make_selfunc.pro ==> make_selection_function.py
            # make_selection_function(template, background, out=None, verbose=False)
            #full_in_selfunc = rootdir+'selec_func/selfunc_backevol/criteres_50-5_40-5/SaveProb_EXT_C2-0-600_summary_CLUSPNTback_'+timestamp[i]
            full_in_selfunc = rootdir+in_selfunc+timestamp[i]
            if(verbose):print 'read ', full_in_selfunc
            selfunctable =  make_selection_function( full_in_selfunc, back[i], verbose=verbose)
            mwrfits( selfunctable, outselfunc)  # can be done inside make_selection_function with option out
            #raise Exception("calculate_xclass_dndcrdhr make_selfunc not yet implemented")

        ###  Construct the good detection CR table
        middlefix = '05-2' ## FIX  we detect in b2  band!
        filtab = [m1filt[i], m2filt[i], pnfilt[i]]
        outdetcrtab = dirdir+'/'+table+'_'+name[i]+'_xclass_detcrtab.fits'
        # /Users/gastaud/icosmo/idldirs_26012012/tmp_test1/blah.test_0001930101_10ks_xclass_detcrtab.fits
        if (os.path.isfile(outdetcrtab)):
            print 'Detection CR table ', outdetcrtab, ' already exists'
        else:
            #  TTAB2D  (150,64) = (n_t, n_z)   is replaced by TEMP (150) =n_t
            #  ZTAB2D  (150,64) = (n_t, n_z)   is replaced by Z     (64) = n_z
            if(verbose):print '  Construct the good detection CR table', outdetcrtab
            detcrtable = make_detcrtab(rootdir+'CRtabs_NHvariable/', middlefix, filtab=filtab, nh=nh, Ttab2D=combined['TEMP'],
                                                  Ztab2D=combined['Z'], verbose=verbose)
            mwrfits( detcrtable, outdetcrtab)
        if (verbose): print 'Detection CR table '+outdetcrtab

        ### Construct the good observable CR table
        #    Force it to be in the Thin/Thin/Thin filter reference frame
        middlefix = mesband 
        tabdir = rootdir+'CRtabs_NHvariable/'
        filtab = ['Thin1','Thin1','Thin1']
        outmescrtab = dirdir+'/'+table+'_'+name[i]+'_'+middlefix+'_xclass_mescrtab.fits'
        if (os.path.isfile(outmescrtab)):
            print 'Measured CR table ', outmescrtab, ' already exists'
        else:
            #  TTAB2D  (150,64) = (n_t, n_z)   is replaced by TEMP (150) =n_t
            #  ZTAB2D  (150,64) = (n_t, n_z)   is replaced by Z     (64) = n_z
            mescrtable = make_detcrtab(rootdir+'CRtabs_NHvariable/', middlefix, filtab=filtab, nh=nh[i], Ttab2D=combined['TEMP'],
                                                  Ztab2D=combined['Z'], verbose=verbose)
            if (verbose):print 'write Measured CR table fits file ', outmescrtab
            mwrfits( mescrtable, outmescrtab)
        if (verbose): print 'Measured CR table '+outmescrtab
        
        ## Construct the good observable CR table indexed in the
        ## vec_cr_bounds vector
        outmescrindextab = dirdir+'/'+table+'_'+name[i]+'_'+middlefix+'_xclass_mescrtab.fits.index'
        if (os.path.isfile(outmescrindextab)):
                print 'Indexed Measured CR table ', outmescrindextab, ' already exists'
        else:
                #  brouillon19
                #raise Exception("make_mesindexcrtab not yet implmented")
                mesindexcrtab = mescrindextable = make_mesindexcrtab(outmescrtab, model_struct['CR_MIN'],model_struct['CR_MAX'],model_struct['N_CR'],
                                                                     save=False, verbose=verbose)
                mwrfits( mesindexcrtab, outmescrindextab)
        if (verbose): print 'Indexed Measured CR table ', outmescrindextab	


        ### Construct the good observable HR table
        #    Force it to be in the Thin/Thin/Thin filter reference frame
        middlefix =  meshr 
        tabdir = rootdir+'HRtabs_NHvariable/'
        filtab = ['Thin1','Thin1','Thin1']
        outmeshrtab = dirdir+'/'+table+'_'+name[i]+'_'+middlefix+'_xclass_meshrtab.fits'
        #print 'good observable HR table=', outmeshrtab
        if (os.path.isfile(outmeshrtab)):
                print 'Measured HR table ', outmeshrtab, ' already exists'
        else:

            #  TTAB2D  (150,64) = (n_t, n_z)   is replaced by TEMP (150) =n_t
            #  ZTAB2D  (150,64) = (n_t, n_z)   is replaced by Z     (64) = n_z
            print 'debug rene ', rootdir+'CRtabs_NHvariable/', middlefix, 'filtab=',filtab, 'nh=', nh[i]
            # for hrtable do not give filtab
            meshrtable = make_detcrtab(rootdir+'CRtabs_NHvariable/', middlefix, nh=nh[i], Ttab2D=combined['TEMP'],
                                                  Ztab2D=combined['Z'], verbose=verbose)
            mwrfits( meshrtable, outmeshrtab)

        if (verbose): print ' good observable HR table ', outmeshrtab	

        #+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        ## Calculate the expected result
        ## Construct the dictionnary (reset XAMIN and FLUXMES) but not the others

        #  xamin good detection
        xamin_struct = {'CRTABFITS':outdetcrtab, 'SFUNCTABFITS':outselfunc, 'SELTYPE':'selfunc', 'STRUCT_NAME':'xamin'}
        #pdb.set_trace()
        # good observable
        new_fluxmes_struct = fluxmes_struct.copy()
        new_fluxmes_struct['CRTABFITS']   = outmescrtab
        new_fluxmes_struct['CRTABFITS_2']   = outmeshrtab
        new_fluxmes_struct['INDEXCRTABFITS']   = outmescrindextab
        new_fluxmes_struct['MES_ERROR']   = -1

        #+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        #compute the result
        # brouillon17 for cmp_dndcrdhr
        result = cmp_dndcrdhr(combined, cosmo, xamin_struct, new_fluxmes_struct, model_struct, verbose=10)
        #print 'result[DNDCRDHR].shape=', result['DNDCRDHR'].shape,'i=', i, area.shape, area[i]
        result['DNDCRDHR'] = area[i]*result['DNDCRDHR']
        #print 'result[DNDCRDHR].shape=', result['DNDCRDHR'].shape, xclassdndcrdhr.shape
        xclassdndcrdhr +=  result['DNDCRDHR']
        #pdb.set_trace()
        filename = dirdir+'/'+table+'_'+name[i]+'_xclass_dCR-'+mesband+'_dHR-'+meshr+'.fits'
        if(save_all): 
            if (os.path.isfile(filename)):
                print 'result file ', filename, ' already exists, overwrite'
                os.remove(filename)
            mwrfits(result, filename )
            if (verbose): print 'result stored in ', filename

    # endfor

    xclasscr = result['CR'] # should be the same for all the pointings
    xclasshr = result['HR']

    # 'NPOINT10',floor(total(timestamp EQ '10ks'))
    # "NPOINT20",floor(total(timestamp EQ '20ks'))

    output = {'CR':xclasscr, 'HR':xclasshr, 'DNDCRDHR':xclassdndcrdhr, 
              'XCLASTAB':table, 'NPOINT':npoint,  "MED_NH":np.median(nh), "MED_BACK":np.median(back), "APERTURE":aperture }

    return output


