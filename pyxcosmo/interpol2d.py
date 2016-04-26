
import numpy as np
import pdb   ### debugger

from scipy.ndimage import map_coordinates



def interpol2d(data, xd, yd, xq, yq, accuracy=1e-5, extrapol=True, verbose=False, **kwargs):
	
	"""
	PURPOSE:
	-------------
	2D interpolation on random points of a function defined on a 
	regular grid
	
	
	ALGORITHM:
	-------------        
	from IDL interpol2d.pro or interpolate: 
	bilinear interpolation
	We have the same result with ndscipy.interpolate.RectBivariateSplin,
	and scipy.ndimage.map_coordinates
	Both are QUICK because they take into account the fact that 
	the input grid is regular.
	
	INPUTS:
	-------------
	(xd, yd) : the regular grid where the function has been evaluated
	must be in increasing order, 2 vectors 1D
	data = the values of the function on the regular grid, matrix 2D
	(xq, yq) = coordinates of the random points where to interpolate, 2 vectors 1D
	
	
	OUTPUT:
	-------------
	Z: the interpolated values for  (xq, yq), vector 1D
	
	HISTORY:
	-------------
	Rene Gastaud, 6 February 2015 works on a non regulard grid using np.digitize
	7 february 2015 better check the ascending order
	9 february 2015 numpy.digitize(x, bins)  with x and bins array like np.atleast_1d
	
	"""
	nx, ny = xd.size, yd.size
	assert (nx, ny) == data.shape
	
	# regular grid
	regular_flag = True
	x_steps = xd[1:]-xd[:-1]
	delta_x_step = x_steps.max()-x_steps.min()
	if (verbose): print 'x_steps', x_steps.min(), x_steps.max(), delta_x_step
	if (delta_x_step > accuracy):
		#if (verbose):print "warning non regular grid on x", delta_x_step
		regular_flag = False
	## ascending order
	assert ( x_steps.min() > 0)
	
	y_steps = yd[1:]-yd[:-1]
	delta_y_step = y_steps.max()-y_steps.min()
	if (verbose): print 'y_steps', y_steps.min(), y_steps.max(), delta_y_step 
	if (delta_y_step > accuracy):
		#if (verbose):print "warning non regular grid on y", delta_y_step 
		regular_flag = False
	## ascending order
	assert ( y_steps.min() > 0)

	## ascending order
	#assert (xd[-1] > xd[0]) and (yd[-1] > yd[0])

	
	# Translate to a unit-cell coordinate system (so that indexes are
	# coordinates)
	# This mapping requires an ascending grid for (xd,yd)

	if (regular_flag):
		#  Compute ip, jp assuming uniform sampling step, and increasing.
		#  x_step = (xd[-1]-xd[0])/(nx-1), 
		xp = (xq-xd[0])*(nx-1)/(xd[-1]-xd[0])
		yp = (yq-yd[0])*(ny-1)/(yd[-1]-yd[0])
		ip = (np.floor(xp)).astype('int')
		jp = (np.floor(yp)).astype('int')
		# we do not use xp, yp
	else:
		# inds = np.digitize(x, bins)
		xq = np.atleast_1d(xq)
		yq = np.atleast_1d(yq)
		ip = np.digitize(xq, xd) -1
		jp = np.digitize(yq, yd) -1
	#
	#pdb.set_trace()
	if (extrapol):
		#if(ip.min()<0): print 'warning extrapolation because xp minimum < 0'
		#if(jp.min()<0): print 'warning extrapolation because yp minimum < 0'
		#if(ip.max()>(nx-2)): print 'warning extrapolation because xp maximum nx-2'
		#if(jp.max()>(ny-2)): print 'warning extrapolation because yp maximum ny-2'
		ip=ip.clip(0, (nx-2))
		jp=jp.clip(0, (ny-2))
	else:
		# check we are inside, no extrapolation
		assert(ip.min()>0)
		assert(jp.min()>0)
		assert(ip.max()<(nx-1))
		assert(jp.max()<(ny-1))
	#

	if (regular_flag):
		t = xp -ip
		u = yp -jp
	else:
		# unequal steps	
		t= (xq-xd[ip])/(xd[ip+1]-xd[ip])
		u= (yq-yd[jp])/(yd[jp+1]-yd[jp])

	res=(1-t)*(1-u)*data[ip,jp] + t*(1-u)*data[ip+1,jp] + (1-t)*u*data[ip,jp+1] + t*u*data[ip+1,jp+1]
	#
	#coord = np.vstack([yp,xp])
	#zq = map_coordinates(data, coord, **kwargs)
	# interp_fun = SI.RectBivariateSpline(xp, yp, proba_Xmes['DNDCRDHR'].T)
	# z = interp_fun.ev(xq, yq)
	#
	return res
