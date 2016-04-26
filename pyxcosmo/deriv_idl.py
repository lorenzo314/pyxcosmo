
import numpy as np

def deriv_idl( x, y):

    """
    Perform numerical differentiation using 3-point, Lagrangian 
    interpolation.

    df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
    Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
  
    Copied from IDL deriv.pro

    """
    n = x.size
    if (x.size != y.size):
        print 'x and y not same size', x.size,y.size
        return -1
    if (n < 3):
        print 'arrays too small', x.size,y.size
        return -1

    x12 = x - np.roll(x,-1)            
    x01 = np.roll(x,1) - x          
    x02 = np.roll(x,1) - np.roll(x,-1) 

    # Middle points
    d = np.roll(y,1) * (x12 / (x01*x02)) + \
          y * (1./x12 - 1./x01) - \
          np.roll(y,-1) * (x01 / (x02 * x12))

    # Formulae for the first and last points:

    # First point
    d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - \
             y[1] * x02[1]/(x01[1]*x12[1]) + \
             y[2] * x01[1]/(x02[1]*x12[1])

    n2 = n-2

    # Last point
    d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) + \
               y[n-2] * x02[n2]/(x01[n2]*x12[n2]) - \
               y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2])

    return d
