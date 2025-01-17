import numpy as np
from scipy.special import binom, loggamma
from numpy.polynomial.laguerre import laggauss
from math import factorial
from functools import partial

# Degree of laguerre polynomials for numerical integration
LAGDEGREE = 40

'''
Impliments the formulas in Zivkovic et al. for calculating second moments of coalescent times
'''

def sigma_i(n, mode, **params):
    '''
    Calculate the expected site frequency spectrum.

    Parameters
    ----------
    n : int
        Sample size
    mode : string
        Type of model. For now may be "constant", "exponential", or "two-epoch".
    params : {str : float}
        Dictionary of parameters, which depend on mode.
        "constant" takes no parameters.
        "exponential" takes a growth rate parameter 'g'.
        "two-epoch" takes two parameters: growth time "tau" and inverse growth factor "f".
    '''
    alpha = zivkovic_alpha(n)
    p_nki = marginal_leaf_prob(n)
    K = np.arange(n+1)
    B = binom(K,2)
    if mode == 'constant':
        I1_j = 1 / B
    else:
        I1_j = laguerre_integral(lambda_inv(mode, **params), B)

    I1_j[np.isinf(I1_j)] = 0
    I1_j[np.isnan(I1_j)] = 0
    sign_jk = diagonal_signs(n)
    alpha_jk = alpha[n]
    T_k = I1_j @ (sign_jk*alpha_jk)
    T_k[:2] = 0
    p_ki = p_nki[n]
    return (K*T_k) @ p_ki

def sigma_ij(n, mode, **params):
    '''
    Calculate the covariance matrix site frequency spectrum.

    Parameters
    ----------
    n : int
        Sample size
    mode : string
        Type of model. For now may be "constant", "exponential", or "two-epoch".
    params : {str : float}
        Dictionary of parameters, which depend on mode.
        "constant" takes no parameters.
        "exponential" takes a growth rate parameter 'g'.
        "two-epoch" takes two parameters: growth time "tau" and inverse growth factor "f".
    '''

    if n > 4 and mode=="two-epoch":
        raise ValueError('Two-epoch growth has not been tested for n>4!')
    Sigma_i = sigma_i(n, mode, **params)
    Ett_kpk = time_second_moments(n, mode, **params)
    Sigma_ij1 = sigma_ij1(n, Ett_kpk)
    Sigma_ij2 = sigma_ij2(n, Ett_kpk)
    Sigma_ij3 = sigma_ij3(n, Ett_kpk)

    Sigma_ij = Sigma_ij1 + Sigma_ij2 + Sigma_ij3
    Sigma_ij += Sigma_ij.T
    Sigma_ij -= np.outer(Sigma_i, Sigma_i)
    Sigma_ij[np.diag_indices(n+1)] = 0
    return Sigma_ij


def zivkovic_alpha(n):
    r = np.arange(n+1)
    vec_n = r[:,None,None]
    vec_j = r[None,:,None]
    vec_k = r[None,None,:]
    ln_alpha = np.log(2*vec_j - 1) + np.log(vec_n) + 2*loggamma(vec_n) \
                + loggamma(vec_k + vec_j - 1) - loggamma(vec_j - vec_k + 1) \
                - np.log(vec_k) - 2*loggamma(vec_k) - loggamma(vec_n - vec_j + 1) \
                - loggamma(vec_n + vec_j)
    alpha = np.exp(np.real(ln_alpha))
    alpha[np.isnan(alpha)] = 0
    return alpha

def lambda_inv(mode, **params):
    funcs = {"constant": lambda x: x,
            "exponential": lambda_inv_exponential,
            "two-epoch": lambda_inv_2epoch}
    return partial(funcs[mode], **params)

def lambda_inv_eq2(mode, **params):
    func = lambda_inv(mode, **params)
    return lambda x, y: (func(y+x) - func(x))**2

def lambda_inv_eq3(mode, **params):
    func = lambda_inv(mode, **params)
    return lambda x, y: func(x) * (func(y+x) - func(x))

def lambda_inv_sq(mode, **params):
    func = lambda_inv(mode, **params)
    return lambda x:  func(x)**2

def lambda_inv_exponential(x, g=0):
    if g == 0:
        return x
    else:
        return np.log(1+x*g)/g

def lambda_inv_2epoch(x, tau=0, f=1):
    return np.piecewise(x, [x>tau], [lambda y: tau + f*(y-tau), lambda y: y])

def laguerre_integral(f, x_scale, degree=LAGDEGREE):
    lag_x, lag_w = laggauss(degree)
    if np.isscalar(x_scale):
        scaled_x = lag_x / x_scale
    elif len(x_scale.shape) == 1:
        scaled_x = lag_x[:,None] / x_scale[None,:]
    return lag_w @ f(scaled_x)

def laguerre_double_integral(f, x_scale, degree=LAGDEGREE):
    lag_x, lag_w = laggauss(degree)
    if np.isscalar(x_scale):
        scaled_x = lag_x / x_scale
        F = f(scaled_x[:,None], scaled_x[None,:])
        print(F)
    else:
        scaled_x = lag_x[None,:] / x_scale[:,None]
        F = f(scaled_x[:,:,None,None], scaled_x[None,None,:,:])
    return lag_w @ (F @ lag_w)

def marginal_leaf_prob(n):
    r = np.arange(n+1)
    vec_n = r[:,None,None]
    vec_k = r[None,:,None]
    vec_i = r[None,None,:]
    p = binom(vec_n-vec_i-1, vec_k-2) / binom(vec_n-1, vec_k-1)
    p[np.isnan(p)] = 0
    p[np.isinf(p)] = 0
    p[:,:,0] = 0
    return p

def pstar(n):
    p_nki = marginal_leaf_prob(n)

    r = np.arange(n+1)
    bin_kkp = binom(r[None,:]-1, r[:,None]-1)
    bin_kkp[np.isnan(bin_kkp)] = 0

    p1_ki = p_nki[n]

    p2_kkpj = np.zeros_like(p_nki)
    for k in range(1,n+1):
        p2_kkpj[k,(k-1):,:] = p_nki[n-k+1,:(n+1)-(k-1),:]

    pstar_kkpij = bin_kkp[:,:,None,None] \
                * p1_ki[:,None,:,None] \
                * p2_kkpj[:,:,None,:]
    return pstar_kkpij

def diagonal_signs(n):
    r = np.arange(n+1)
    return (-1)**(r[:,None]+r[None,:])

### <T_k' T_k> ###
def time_second_moments(n, mode, **params):
    alpha = zivkovic_alpha(n)
    moments = np.zeros((n+1, n+1))
    moments[np.tril_indices(n+1)] = sec_moments_off_diag(n, alpha, mode, **params)
    moments[np.diag_indices(n+1)] = sec_moments_diagonal(n, alpha, mode, **params)
    moments[n,n] = sec_moments_corner(n, mode, **params)
    return moments

def sec_moments_off_diag(n, alpha, mode, **params):
    # 2 <= k < k' <= n
    sign_kpk = diagonal_signs(n)
    sign_ji = diagonal_signs(n)
    B = binom(np.arange(n+1), 2)

    prefactor_ji = sign_ji * (B[:,None] - B[None,:]) / B[:,None]
    prefactor_ji[:,:2] = 0
    prefactor_ji[np.triu_indices(n+1)] = 0

    # G_ji = laguerre_double_integral(lambda_inv_eq3(mode, **params), B)
    if mode == 'constant':
        G_ji = 1 / np.outer(B,B)
    else:
        G_ji = laguerre_double_integral(lambda_inv_eq3(mode, **params), B)
    #G_ji = laguerre_double_integral(partial(lambda_inv_eq3, g=g), B)
    I_ji = prefactor_ji * G_ji
    I_ji[:,:2] = 0
    I_ji[np.triu_indices(n+1)] = 0

    A_kpjik = (alpha[n].T)[:,:,None,None] * alpha[:,None,:,:]
    ret = sign_kpk * np.tensordot(A_kpjik, I_ji, axes=([1,2],[0,1]))
    ret[:,:2] = 0
    return ret[np.tril_indices(n+1)]

def sec_moments_diagonal(n,alpha, mode, **params):
    # k' = k < n
    r = np.arange(n+1)
    j_vec = r[:,None]
    k_vec = r[None,:]
    prefactor_jk = (-1)**(j_vec+k_vec+1) * binom(k_vec+1, 2) / binom(j_vec, 2)
    prefactor_jk[:,:2] = 0
    prefactor_jk[np.triu_indices(n+1)] = 0
    alpha_jk = np.roll(alpha[n], -1, axis=1)

    B = binom(r,2)
    # G_jk = laguerre_double_integral(lambda_inv_eq2(mode, **params), B)
    # print(G_jk)
    # print(G_jk)
    if mode == 'constant':
        G_jk = (2 / (B**2))[None,:]
    else:
        G_jk = laguerre_double_integral(lambda_inv_eq2(mode, **params), B)
    # G_jk = laguerre_double_integral(partial(lambda_inv_eq2, g=g), binom(r,2))
    return np.nansum(G_jk*alpha_jk*prefactor_jk, axis=0)

def sec_moments_corner(n, mode, **params):
    # k' = k = n
    return laguerre_integral(lambda_inv_sq(mode, **params), binom(n,2))
    # return laguerre_integral(partial(lambda_inv_sq, g=g), binom(n,2))

def sigma_ij1(n, Ett_kpk):
    pstar_kkpij = pstar(n)

    r = np.arange(n+1)
    k_vec = r[None,:]
    kp_vec = r[:,None]
    Sigma_ij1 = np.zeros((n+1,n+1))
    for u in [1,2]:
        prefactor_kpk = (-1)**u * k_vec*(k_vec-1)/binom(kp_vec-1, k_vec-u)
        prefactor_kpk[:,:2] = 0
        prefactor_kpk[np.triu_indices(n+1)] = 0
        A_kpk = prefactor_kpk * Ett_kpk
        B_kkpij = np.roll(pstar_kkpij, u-2, axis=0)
        Sigma_ij1 += np.tensordot(A_kpk, B_kkpij, axes=([1,0],[0,1]))
    Sigma_ij1[np.triu_indices(n+1)] = 0
    return Sigma_ij1

def sigma_ij2(n, Ett_kpk):
    term1 = sigma_ij2a(n, Ett_kpk)
    term2 = sigma_ij2b(n, Ett_kpk)
    return term1 + term2

def sigma_ij2a(n, Ett_kpk):
    p_ki = marginal_leaf_prob(n)[n]
    p_kij = np.zeros((n+1,n+1,n+1))
    for j in range(1, n):
        p_kij[:,1:,j] = np.pad(p_ki[:-1,j:],
                                ((1,0),(0,j)),
                                mode='constant')[:,1:]
    k = np.arange(n+1)
    Et2_k = Ett_kpk[np.diag_indices(n+1)]
    A_k = Et2_k * k * (k-1)**2 / (n-k+1)
    ret = np.tensordot(A_k[3:], p_kij[3:], axes=(0,0))
    ret[np.triu_indices(n+1)] = 0
    return ret

def sigma_ij2b(n, Ett_kpk):
    pstar_kkpij = pstar(n)

    r = np.arange(n+1)
    kp = r[:,None]
    k  = r[None,:]
    A_kpk = Ett_kpk * k * (k-1) / binom(kp-1, k-1)
    A_kpk[np.isnan(A_kpk)] = 0
    A_kpk[np.triu_indices(n+1)] = 0

    B_kkpij = np.zeros((n+1,n+1,n+1,n+1))
    for u in [1,2]:
        for v in [1,2]:
            pstar_shift = np.roll(pstar_kkpij, u-v, axis=0)
            for j in range(1,n):
                B_kkpij[:,:,:(n+1)-j,j] += np.diagonal(pstar_shift,
                                                offset=-j, axis1=2, axis2=3) \
                                        + pstar_shift[:,:,j:,j]
    ret = np.tensordot(A_kpk, B_kkpij, ([1,0],[0,1]))
    ret[np.triu_indices(n+1)] = 0
    return ret

def sigma_ij3(n, Ett_kpk):
    kp = np.arange(n+1)
    p_kpi = marginal_leaf_prob(n)[n]
    pfold_kpi = p_kpi + p_kpi[:,::-1]
    val_i = np.zeros(n+1)
    val_i[(n//2)+1:n] = (2*Ett_kpk[2,2])/(n-1) \
                      + ((2*Ett_kpk[3:,2]/(kp[3:]-1)) \
                      @ pfold_kpi[3:,(n//2)+1:n])
    ret = np.diag(val_i)
    return np.fliplr(ret)



def main():
    n = 5

    print("Constant")
    mode = "constant"
    params = dict()
    Sigma_i = sigma_i(n, mode, **params)
    Sigma_ij = sigma_ij(n, mode, **params)
    print("SFS:\n", Sigma_i, '\n')
    print("Sigma_ij:\n", Sigma_ij, '\n')

    print("Exp. growth, g=1")
    mode = "exponential"
    params = {"g":1}
    Sigma_i = sigma_i(n, mode, **params)
    Sigma_ij = sigma_ij(n, mode, **params)
    print("SFS:\n", Sigma_i, '\n')
    print("Sigma_ij:\n", Sigma_ij, '\n')

    print("Two-epoch, tau=1, f=0.5")
    mode = "two-epoch"
    params = {"tau":1, "f":0.5}
    Sigma_i = sigma_i(n, mode, **params)
    Sigma_ij = sigma_ij(n, mode, **params)
    print("SFS:\n", Sigma_i, '\n')
    print("Sigma_ij:\n", Sigma_ij, '\n')

if __name__ == '__main__':
    main()
