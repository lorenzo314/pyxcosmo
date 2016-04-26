# IPython log file

def extract_dict(dict, i, verbose=False):
    """
        extract a sub dictionnary
        
        Arguments
        ---------
        dict: dictonnary
        n_z: int
        n_m: int
        verbose: bool
        
        Returns
        -------
        A sub dictionary of arrays
        
        Update:
        ------
        Rene Gastaud
        Samedi 30 janvier 2016
        
        """
    if verbose:
        print 'my extension = ', i
    n_z = dict['Z'].size
    n_m = dict['M'].size/n_z
    diff = dict['M'][0:n_m] - dict['M'][n_m:2*n_m]
    flag = False
    if (diff.min() < 0)or (diff.max() > 0):
        n_m = dict['M'].size
        flag = True
    #print "extract_dict6", n_z, n_m, flag
    if verbose: print 'nz',n_z, 'n_m', n_m, flag

    # load data
    subdict = {}
    for key in dict.keys():
        if (verbose):print key, dict[key].shape
        if (dict[key].size == n_m*n_z):
            if (verbose):print ' n_m n_z'
            if(dict[key].ndim == 2):
                subdict[key]=dict[key][i, :]
            else:
                subdict[key]=dict[key][i*n_m:(i+1)*n_m]
        if (dict[key].size == n_z):
            if (verbose):print ' n_z', n_z
            subdict[key]=dict[key][i]
        if (dict[key].size == n_m):
            if (verbose):print ' n_m', n_m
            subdict[key]=dict[key]

    return subdict
