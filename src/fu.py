import numpy as np

'''
Implement formulas for the moments of the number of leaves subtended by branches
at different levels from Fu (1995) as reproduced in Rafajlovic et al. (2014).
'''

def marginal_leaf_prob(n):
    '''
    Calculate the probability that a branch at level k has i leaves.
    Return as an (n-1)x(n-1) array.
    '''
    p_ki = np.zeros((n-1,n-1))
    return p_ki

def joint_leaf_prob(n):
    '''
    Calculate the joint probability that a branches at level k and m
    have i and j leaves respectively.
    Return as an (n-1)^4 array.
    '''
    p_kimj = np.zeros((n-1,n-1,n-1,n-1))
    return p_kimj
