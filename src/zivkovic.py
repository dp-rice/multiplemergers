import numpy as np
from scipy.special import binom, loggamma
from numpy.polynomial.laguerre import laggauss
from math import factorial
from functools import partial

n = 4
g = 0

k = np.arange(n+1)

b = binom(k,2)
b_diff = b[:,None] - b[None,:]
b_diff[b_diff<0] = 0
print(b_diff)

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

alpha = zivkovic_alpha(n)
p_nki = marginal_leaf_prob(n)


I1_j = laguerre_integral(partial(lambda_inv, g=g), b)
I1_j[np.isinf(I1_j)] = 0
sign_jk = (-1)**(k[:,None]+k[None,:])
alpha_jk = alpha[n]
T_k = I1_j @ (sign_jk*alpha_jk)
T_k[:2] = 0
p_ki = p_nki[n]
sfs_i = (k*T_k) @ p_ki

print(sfs_i)

# print(laguerre_double_integral(lambda_inv_eq2, b))
# print(2/b**2)
#
# print(laguerre_double_integral(lambda_inv_eq3, b))
# print(1/(b[None,:]*b[:,None]))

prefactor = sign_jk * (b[:,None] - b[None,:]) / b[:,None]
prefactor[:,:2] = 0
prefactor[np.triu_indices(n+1)] = 0
I_ji = prefactor * laguerre_double_integral(partial(lambda_inv_eq3, g=g), b)
I_ji[:,:2] = 0
I_ji[np.triu_indices(n+1)] = 0

A_kpjk = I_ji @ alpha
# TODO: vectorize
Ett_kpk = np.zeros((n+1,n+1))
for kp in k:
    Ett_kpk[kp] = sign_jk[kp,:] * (alpha[n,:,kp] @ A_kpjk[kp])
Ett_kpk[:,:2] = 0
Ett_kpk[np.triu_indices(n+1)] = 0

print(Ett_kpk)
