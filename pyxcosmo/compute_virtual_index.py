def compute_virtual_index(grid, values, flag=True, threshold=1e-6):
    """
        Compute the virtual indexes.

        INPUTS
        grid position of each sample

        OUTPUTS
        indexes 

        EXAMPLE
        index = compute_virtual_index(grid, values)

        R Gastaud fecit,  24 oct 2014

    """
    x = grid
    nx = x.size
    steps = x[1:]-x[0:-1]
    diff = steps.max()-steps.min()
    if (diff > threshold):
        print 'error in compute_virtual_index', diff, theshold
        return -1

    step=steps.mean()
    origin = x[0]
      
    virtual_i = (values-origin)/step
    virtual_i =  virtual_i.clip(min=0, max=nx-1)

    return  virtual_i

