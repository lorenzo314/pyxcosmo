
import numpy as np
import pdb   ### debugger

#from scipy.ndimage import map_coordinates
from scipy.interpolate import RectBivariateSpline


def int_tabulated_2d(x, y, data,  verbose=False):
	
	"""
	PURPOSE:
	-------------
	double integral with sampled data python
	
	ALGORITHM:
	-------------        
	use np.trapz
	see http://stackoverflow.com/questions/24875722/double-integral-with-function-and-sampled-data-python

	INPUTS:
	-------------
	(x, y) : the regular grid where the function has been evaluated
	must be in increasing order, 2 vectors 1D
	data = the values of the function on the regular grid, matrix 2D	
	
	OUTPUT:
	-------------
	Z: the interpolated values for  (xq, yq), table 2 D
	
	HISTORY:
	-------------
	Rene Gastaud
	10 february 2015 
	
	"""
	nx, ny = x.size, y.size
	assert (nx, ny) == data.shape
	if (verbose): print nx, ny

	x_steps = x[1:]-x[:-1]
	delta_x_step = x_steps.max()-x_steps.min()
	if (verbose): print 'x_steps', x_steps.min(), x_steps.max(), delta_x_step
	## ascending order
	assert ( x_steps.min() > 0)
	
	y_steps = y[1:]-y[:-1]
	delta_y_step = y_steps.max()-y_steps.min()
	if (verbose): print 'y_steps', y_steps.min(), y_steps.max(), delta_y_step 
	## ascending order
	assert ( y_steps.min() > 0)

	# now integrate in x ==i
	integral_x = np.zeros( ny )

	for j in range(ny):
		integral_x[j] = np.trapz( data[:,j], x)
	# then an integral over the result   
	integral = np.trapz( integral_x, y )

	#
	return 	integral
