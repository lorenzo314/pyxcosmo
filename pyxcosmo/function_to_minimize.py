from make_dndmdz import make_dndmdz
from calc_likelihood_dndmdz import calc_likelihood_dndmdz

# Function to maximize the likelihood
def function_to_minimize(x, param, model, verbose=False):
    # Compute dndmdz
    param['OMEGA_M'] = x[0]
    param['SIGMA8']  = x[1]
    dndmdz = make_dndmdz(param)

    if verbose:
        print "omega_m", x[0]
        print "sigma8", x[1]

    # Compute likelihood
    likelihood = calc_likelihood_dndmdz(model, dndmdz)

    return -likelihood
