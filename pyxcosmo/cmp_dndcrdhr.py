
from  check_param import check_param
from dndzdldt_to_dndcrdhr import dndzdldt_to_dndcrdhr
import numpy as np
import math

__version__ = 3

# LF 08/12/15: see cmp_dndcrdhr in /dsm/sappccosmos0/work/lfacciol/faites_par_nicolas
def cmp_dndcrdhr(dndzdldt, cosmo, xamin_struct, fluxmes_struct, model_struct, 
		 verbose=False):
        """
        Short Summary
        -------------
        Computes dndcrdhr from dndzdldt, mainly calling
        dndzdldt_to_dndcrdhr

        Inputs
        ----------
        dndzdldt: 
        cosmo:
	xamin_struct:
        fluxmes_struct:
        model_struct:

	Needed Files:
	crtab, crtab_2, sfunctab, indextab, and the histo_integer.so

        LF: crtab_2 MEANS HR !!!!!!!

        Returns
        -------
        dndcrdhr: a dictonnary containing a matrix  dndcrdhr 
        versus a vector CR and a vector HR

	HISTORY:
	7 May 2015 : conversion of seltype in string only if needed, version 3
	4 May 2015 use check_param and not check_paral2 ,
	remove MT_struct and LT_struct version 2
	check_param2 returns a status
	check_param returns an update dictionary
        """

	##  check the items and read the tables crtab, crtab_2, sfunctab, indextab
	# xamin read tables
	xamin_struct = check_param(xamin_struct, cosmo, type='xamin', 
				   verbose=verbose)

	fluxmes_struct = check_param(fluxmes_struct, cosmo, type='fluxmes', 
				     verbose=verbose)

	model_struct = check_param(model_struct, cosmo, type='model', 
				   verbose=verbose)


	## extract the parameters from the dictionnaries

	crtab = xamin_struct['CRTAB']

	crtab_2 =  fluxmes_struct['CRTAB_2']
	mes_hr = crtab_2['CR_PERL_TIMESD2']

	mesindexcrtab = fluxmes_struct['INDEXCRTAB']
	index_mes_cr = mesindexcrtab['INDEX_CR_PERL_TIMESD2']

	# 7 May 2015 : conversion of seltype in string only if needed, version 3
	seltype = xamin_struct['SELTYPE']
	if (isinstance(seltype, np.ndarray)):
	    seltype = np.ndarray.tostring(seltype)
	
	aperture =  fluxmes_struct['APERTURE']

	sfunctab = xamin_struct['SFUNCTAB']


	result = dndzdldt_to_dndcrdhr(dndzdldt, cosmo['H'], seltype, crtab, 
		      sfunctab, model_struct, aperture, index_mes_cr, mes_hr)

	result = result*math.radians(1)*math.radians(1)

	new = {'HR':dndzdldt['VEC_HR'], 'CR':dndzdldt['VEC_CR'], 
	       'DNDCRDHR':result }

	return new
