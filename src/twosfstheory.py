import numpy as np
import fu

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
                p_ki, p_kimj = fu_tensors[self.n]
            except KeyError:
                p_ki = fu.marginal_leaf_prob(self.n)
                p_kimj = fu.joint_leaf_prob(self.n)
                fu_tensors[self.n] = (p_ki, p_kimj)
            # Calculate first and second moments of coalescence times.
            t_k = time_first_moments(self.n, growth_rate)
            t_km = time_second_moments(self.n, growth_rate)
            # TODO: Check order of dot product
            self.sfs = np.dot(p_ki, t_k)
            self.two_sfs = np.dot(p_kimj, t_km)

        # For beta coalescent, use Birkner eqns ???
        elif alpha is not None:
            # FIXME:
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
        return np.binom(n,k)
    else:
        return

def time_second_moments(n, r):
    '''
    Calculate second moments of coalescence times using Eriksson eq 14.
    applied to an exponentially growing population with rate r
    '''
    return
