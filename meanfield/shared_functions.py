import numpy as np

def w_function(x,w_0,w_1,w_sigma):
    """Gaussian"""
    return w_0 + w_1 * np.exp(- x**2 /( 2. * w_sigma**2 ))

def g_function(x,g_0,g_1,g_sigma,g_r):
    """Generalized gaussian with shape factor g_r"""
    return g_0 + g_1 * np.exp(-(abs(x)/g_sigma)**g_r)