import numpy as np
import fu
from scipy.special import binom

class TwoSFS():
    '''
    Calculate the two-site frequency spectrum for perfectly-linked sites.
    Uses formulas from Fu (1995), Eriksson et al. (2010), and Rafajlovic et al. (2014).
    Currently only implemented for exponential growth.
    '''

    # Class variable that stores parameters p_ki and p_kimj
    # These are independent of the demographic model.
    fu_tensors = {}

    def __init__(self, sample_size, growth_rate=None, alpha=None):
        self.n = sample_size

        if growth_rate and alpha:
            raise ValueError('May not specify both growth rate and alpha.')

        # For exponential growth, use Rafajlovic eqns B.1-B.3
        elif growth_rate is not None:
            # If already calculated, retreive fu tensors
            # If not, calculate them and store them.
            try:
                p_ki, p_kimj = TwoSFS.fu_tensors[self.n]
            except KeyError:
                print('Calculating...')
                p_ki = fu.marginal_leaf_prob(self.n)
                p_kimj = fu.joint_leaf_prob(self.n)
                TwoSFS.fu_tensors[self.n] = (p_ki, p_kimj)
            # Calculate first and second moments of coalescence times.
            t_k = time_first_moments(self.n, growth_rate)
            t_km = time_second_moments(self.n, growth_rate)

            # Rafajovic B.1
            k = np.arange(2, self.n+1)
            self.sfs = np.dot(p_ki, k * t_k) / 2
            # Rafajovic B.2
            count_factor = k[:,None]*k[None,:] - np.diag(k)
            # TODO: Check order of dot product
            self.two_sfs = np.dot(p_kimj, count_factor * t_km) / 4

        # For beta coalescent, use Birkner eqns ???
        elif alpha is not None:
            # TODO:
            raise ValueError('Beta-coalescent not implemented yet.')
        else:
            raise ValueError('Must specify either growth_rate or alpha.')

    def get_sfs(self, folded=False):
        if not folded:
            return self.sfs
        else:
            return fold_sfs(self.sfs)

    def get_2sfs(self, folded=False):
        if not folded:
            return self.two_sfs
        else:
            return fold_two_sfs(self.two_sfs)


def fold_sfs(sfs):
    '''Fold SFS for unidentified ancestral state.'''
    # TODO:
    return

def fold_two_sfs(two_sfs):
    '''Fold 2SFS for unidentified ancestral state.'''
    # TODO:
    return

def time_first_moments(n, r):
    '''
    Calculate expected coalescence times using Eriksson eq 9.
    applied to an exponentially growing population with rate r
    '''
    k = np.arange(2,n+1)
    if r == 0:
        return 1.0 / binom(k,2)
    else:
        # TODO:
        return np.zeros(n-1)

def time_second_moments(n, r):
    '''
    Calculate second moments of coalescence times using Eriksson eq 14.
    applied to an exponentially growing population with rate r
    '''
    # TODO:
    return np.zeros((n-1,n-1))

if __name__ == '__main__':
    n = 10
    k = np.arange(2,n+1)
    i = np.arange(1,n)
    tsfs = TwoSFS(n, growth_rate=0)
    print(tsfs.get_sfs() * i)
