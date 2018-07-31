import numpy as np
from scipy.special import binom

'''
Implement formulas for the moments of the number of leaves subtended by branches
at different levels from Fu (1995) as reproduced in Rafajlovic et al. (2014).
'''

def marginal_leaf_prob(n):
    '''
    Calculate the probability that a branch at level k has i leaves.
    Return as an (n-1)x(n-1) array.
    '''
    k = np.arange(2, n+1)
    i = np.arange(1, n)
    num = (k[None,:]-1) * binom(n-k[None,:], i[:,None]-1)
    den = i * binom(n-1, i)
    return num / den[:,None]

def joint_leaf_prob(n):
    '''
    Calculate the joint probability that a branches at level k and m
    have i and j leaves respectively.
    Return as an (n-1)^4 array.
    '''
    p_kimj = np.zeros((n-1,n-1,n-1,n-1))
    return p_kimj
