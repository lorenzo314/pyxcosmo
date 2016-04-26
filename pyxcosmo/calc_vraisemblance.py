import numpy as np
import math
import pdb

from scipy.interpolate import RectBivariateSpline

from fitsutils import mrdfits
from fitsutils import mwrfits
from convol_2d_variable import convol_2d_variable
from interpol2d import interpol2d
from int_tabulated_2d import int_tabulated_2d

def calc_vraisemblance(observable,model,errormodel=None, verbose=False):

   """
   PURPOSE:
   -------------
   compute the likelihood

   INPUTS:
   -------------
    (observable) : a dictionnary containing:
      * the catalogue of CR{i}
      * the catalogue of HR{i}
      * a grid (CRgrid,HRgrid) on which is defined the mean bias
          int(  1-f+f*4pi/Surfin*alpha*mu_in * P(f) df)
          for each point of the grid

    (model)


    (errormodel): a dictionnary or a string, name of the fits file contaiting the dictionnary


   OUTPUT:
   -------------
    Z: likelihood 

   HISTORY:
   -------------
    Rene Gastaud, 7 february 2015 
   """
   # Convolve the model with the measurement error
   model_cr = model['CR']
   model_hr = model['HR']

   log_obs_cr = np.log(observable['CR'])
   log_obs_hr = np.log(observable['HR'])

   if(errormodel != None):
      print '*** calc_vraisemblance : convolving error model... : '
      if (errormodel.__class__ == str):
         errmod = mrdfits(errormodel,1)
      else:
         errmod = errormodel
      proba_Xmes = convol_2d_variable(model,errmod)
   else:
      proba_Xmes = model

   #----------------------------------------
   #
   # First term of the likelihood : the sum of the probabilities
   n_amas = observable['L4SDBID'].size

   #for i in np.range(0, n_amas):
   #  proba_amas, expected_bias, proba_biased, like_amas are 1D vector, size n_amas  
   proba_amas = interpol2d(proba_Xmes['DNDCRDHR'].T, np.log(proba_Xmes['CR']), np.log(proba_Xmes['HR']), log_obs_cr , log_obs_hr )
   if (verbose): print 'proba amas', proba_amas.min(), proba_amas.max()
   proba_amas = proba_amas.clip(min=0)

   expected_bias = interpol2d(observable['EXPECTED_BIAS'].T, np.log(observable['CRGRID']), np.log(observable['HRGRID']), log_obs_cr , log_obs_hr )
   # Multiply by the bias
   proba_biased = expected_bias * proba_amas
   like_amas    = np.log(proba_biased)
   # now we sum
   sum = like_amas.sum()

   #
   # Second term of the likelihood : the constant term
   #
   # Attention, on n'integre que la partie qui a ete observee
   # (i.e. si on enleve plein de points a gauche et a droite en CR, on
   # ne les compte pas, bien sur) en plus biais mal connu...

   max_obs_cr = observable['CRMAX']
   min_obs_cr = observable['CRMIN']

   max_obs_hr = observable['HRMAX']
   min_obs_hr = observable['HRMIN']

   if (verbose):print min_obs_cr, max_obs_cr, max_obs_hr, min_obs_cr

   (imincr, imaxcr) = np.digitize(np.array([min_obs_cr[0], max_obs_cr[0]]),model_cr) -1
   (iminhr, imaxhr) = np.digitize(np.array([min_obs_hr[0], max_obs_hr[0]]),model_hr) -1

   interieur_cr = model_cr[imincr+1:imaxcr+1] # the last index is excluded in python
   restricted_cr_mes = np.concatenate((min_obs_cr, interieur_cr), axis=0)
   if(model_cr[imaxcr] != max_obs_cr): restricted_cr_mes =  np.append(restricted_cr_mes, max_obs_cr)

   interieur_hr = model_hr[iminhr+1:imaxhr+1] # the last index is excluded in python
   restricted_hr_mes = np.concatenate((min_obs_hr, interieur_hr), axis=0)
   if(model_hr[imaxcr] != max_obs_hr): restricted_hr_mes =  np.append(restricted_hr_mes, max_obs_hr)

   log_restricted_cr_mes = np.log(restricted_cr_mes)
   log_restricted_hr_mes = np.log(restricted_hr_mes)

   # here RectBivariateSpline better and quicker than interpol2d
   #   to be checked (cloked)
   ## 
   interp_all = RectBivariateSpline(np.log(proba_Xmes['HR']), np.log(proba_Xmes['CR']),proba_Xmes['DNDCRDHR'], kx=1, ky=1)

   proba_all = interp_all(log_restricted_hr_mes, log_restricted_cr_mes)

   #
   # can not use RectBivariateSpline because we have extrapolation
   (xx,yy) = np.meshgrid(restricted_cr_mes, restricted_hr_mes)
   log_xx = np.log(xx)
   log_yy = np.log(yy)
   all_expected_bias = interpol2d(observable['EXPECTED_BIAS'].T, np.log(observable['CRGRID']), np.log(observable['HRGRID']), log_xx, log_yy )
   
   all_proba_biased = all_expected_bias * proba_all
   ### integration
   n_cr = restricted_cr_mes.size
   n_hr = restricted_hr_mes.size
   data = all_proba_biased*restricted_cr_mes.reshape(1, n_cr)*restricted_hr_mes.reshape(n_hr,1)
   constant_term = int_tabulated_2d(log_restricted_hr_mes, log_restricted_cr_mes, data)
   if (verbose):print 'constant_term', constant_term
   #The likelihood
   res = sum - constant_term

   return -2*res
