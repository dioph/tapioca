import numpy as np
from astroquery.mast import Observations

def get_tpf_kepler(target, quarter=None, campaign=None):
    """
    returns table of TPFs from Kepler or K2 for a given target
    """
    
    if 0 < target < 2e8:
        name = 'kplr{:09d}'.format(target)
    elif 2e8 < target < 3e8:
        name = 'ktwo{:09d}'.format(target)
    else:
        # implement error handling function
        pass
        
    obs = Observations.query_criteria(target_name=name, project=['Kepler', 'K2'])
    products = Observations.get_product_list(obs)
    suffix = 'lpd-targ'
    mask = np.array([suffix in fn for fn in products['productFilename']])
    when = campaign if campaign is not None else quarter
    if when is not None:
        mask &= np.array([desc.lower().endswith('q{}'.format(when)) or
                            desc.lower().endswith('c{:02}'.format(when)) or
                            desc.replace('-','').lower().endswith('c{:03d}'.format(when))
                            for desc in products['description']])
    return products[mask]
    
def get_lc_kepler(target, quarter=None, campaign=None):
    """
    returns table of TPFs from Kepler or K2 for a given target
    """
    
    if 0 < target < 2e8:
        name = 'kplr{:09d}'.format(target)
    elif 2e8 < target < 3e8:
        name = 'ktwo{:09d}'.format(target)
    else:
        # implement error handling function
        pass
        
    obs = Observations.query_criteria(target_name=name, project=['Kepler', 'K2'])
    products = Observations.get_product_list(obs)
    suffix = 'Lightcurve Long'
    mask = np.array([suffix in fn for fn in products['productFilename']])
    when = campaign if campaign is not None else quarter
    if when is not None:
        mask &= np.array([desc.lower().endswith('q{}'.format(when)) or
                            desc.lower().endswith('c{:02}'.format(when)) or
                            desc.replace('-','').lower().endswith('c{:03d}'.format(when))
                            for desc in products['description']])
    return products[mask]
    
def acf(y):
    N = len(y)
    m = np.mean(y)
    s = np.square(y-m).sum()
    R = np.zeros(N)
    for h in range(N):
        a = y[h:]
        b = y[:N-h]
        R[h] = ((a-m)*(b-m)).sum() / s
    return R
    
    
