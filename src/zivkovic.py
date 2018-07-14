import numpy as np
from scipy.special import binom, loggamma
from numpy.polynomial.laguerre import laggauss
from math import factorial
from functools import partial

n = 4
g = 0

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

def lambda_inv(x, g=0):
    if g == 0:
        return x
    else:
        return np.log(1+x*g)/g

def lambda_inv_eq2(x,y,g=0):
    return (lambda_inv(y+x,g) - lambda_inv(x,g))**2

def lambda_inv_eq3(x,y,g=0):
    return lambda_inv(x,g) * (lambda_inv(y+x,g) - lambda_inv(x,g))

def lambda_inv_sq(x, g=0):
    return lambda_inv(x, g=g)**2

def laguerre_integral(f, x_scale, degree=10):
    lag_x, lag_w = laggauss(degree)
    if np.isscalar(x_scale):
        scaled_x = lag_x / x_scale
    elif len(x_scale.shape) == 1:
        scaled_x = lag_x[:,None] / x_scale[None,:]
    return lag_w @ f(scaled_x)

def laguerre_double_integral(f, x_scale,  degree=10):
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

alpha = zivkovic_alpha(n)
p_nki = marginal_leaf_prob(n)
pstar_kkpij = pstar(n)

K = np.arange(n+1)
B = binom(K,2)
### <SFS> ###

I1_j = laguerre_integral(partial(lambda_inv, g=g), B)
I1_j[np.isinf(I1_j)] = 0
I1_j[np.isnan(I1_j)] = 0
sign_jk = diagonal_signs(n)
alpha_jk = alpha[n]
T_k = I1_j @ (sign_jk*alpha_jk)
T_k[:2] = 0
p_ki = p_nki[n]
sfs_i = (K*T_k) @ p_ki

print(sfs_i)

### <T_k' T_k> ###

# 2 <= k < k' <= n
sign_kpk = diagonal_signs(n)
sign_ji = diagonal_signs(n)

prefactor_ji = sign_ji * (B[:,None] - B[None,:]) / B[:,None]
prefactor_ji[:,:2] = 0
prefactor_ji[np.triu_indices(n+1)] = 0
G_ji = laguerre_double_integral(partial(lambda_inv_eq3, g=g), B)
I_ji = prefactor_ji * G_ji
I_ji[:,:2] = 0
I_ji[np.triu_indices(n+1)] = 0

A_kpjik = (alpha[n].T)[:,:,None,None] * alpha[:,None,:,:]
Ett_kpk = sign_kpk * np.tensordot(A_kpjik, I_ji, axes=([1,2],[0,1]))
Ett_kpk[:,:2] = 0
Ett_kpk[np.triu_indices(n+1)] = 0

# k' = k < n
j_vec = K[:,None]
k_vec = K[None,:]
prefactor_jk = (-1)**(j_vec+k_vec+1) * binom(k_vec+1, 2) / binom(j_vec, 2)
prefactor_jk[:,:2] = 0
prefactor_jk[np.triu_indices(n+1)] = 0
alpha_jk = np.roll(alpha[n], -1, axis=1)
G_jk = laguerre_double_integral(partial(lambda_inv_eq2, g=g), B)
Ett_kpk[np.diag_indices(n+1)] = np.nansum(G_jk*alpha_jk*prefactor_jk, axis=0)

# k' = k = n
Ett_kpk[n,n] = laguerre_integral(partial(lambda_inv_sq, g=g), B[n])

### SIGMA_ij1 ###

k_vec = K[None,:]
kp_vec = K[:,None]
Sigma_ij1 = np.zeros((n+1,n+1))
for u in [1,2]:
    prefactor_kpk = (-1)**u * k_vec*(k_vec-1)/binom(kp_vec-1, k_vec-u)
    prefactor_kpk[:,:2] = 0
    prefactor_kpk[np.triu_indices(n+1)] = 0
    A_kpk = prefactor_kpk * Ett_kpk
    B_kkpij = np.roll(pstar_kkpij, u-2, axis=0)
    Sigma_ij1 += np.tensordot(A_kpk, B_kkpij, axes=([1,0],[0,1]))
Sigma_ij1[np.triu_indices(n+1)] = 0
print(Sigma_ij1)

### SIGMA_ij2 ###

Sigma_ij2 = np.zeros((n+1,n+1))

### SIGMA_ij3 ###

Sigma_ij3 = np.zeros((n+1,n+1))
