
import numpy as np
import pdb   ### debugger

#from scipy.ndimage import map_coordinates
from scipy.interpolate import RectBivariateSpline


def interpol2d_2d(data, xd, yd, xq, yq, accuracy=1e-5, kx=1, ky=1, verbose=False):
	
	"""
	PURPOSE:
	-------------
	2D interpolation on a grid of a function defined on a 
	regular grid
	
	
	ALGORITHM:
	-------------        
	use SI.RectBivariateSpline.
	see https://github.com/scipy/scipy/issues/3164

	INPUTS:
	-------------
	(xd, yd) : the regular grid where the function has been evaluated
	must be in increasing order, 2 vectors 1D
	data = the values of the function on the regular grid, matrix 2D
	(xq, yq) = coordinates of the random points where to interpolate, 2 vectors 1D
	
	
	OUTPUT:
	-------------
	Z: the interpolated values for  (xq, yq), table 2 D
	
	HISTORY:
	-------------
	Rene Gastaud
	10 february 2015 
	
	"""
	nx, ny = xd.size, yd.size
	assert (nx, ny) == data.shape
	
	# regular grid
	regular_flag = True
	x_steps = xd[1:]-xd[:-1]
	delta_x_step = x_steps.max()-x_steps.min()
	if (verbose): print 'x_steps', x_steps.min(), x_steps.max(), delta_x_step
	if (delta_x_step > accuracy):
		if (verbose):print "warning non regular grid on x", delta_x_step
	## ascending order
	assert ( x_steps.min() > 0)
	
	y_steps = yd[1:]-yd[:-1]
	delta_y_step = y_steps.max()-y_steps.min()
	if (verbose): print 'y_steps', y_steps.min(), y_steps.max(), delta_y_step 
	if (delta_y_step > accuracy):
		if (verbose):print "warning non regular grid on y", delta_y_step 
	## ascending order
	assert ( y_steps.min() > 0)

	# to to: test xq, yq

	print ' kkk', kx, ky
	print 'shapes', data.shape, xd.shape, yd.shape, xq.shape, yq.shape
	interp_fun = RectBivariateSpline(yd, xd, data, kx=kx, ky=ky)
	z = interp_fun.ev(yq, xq)
	#
	return z
