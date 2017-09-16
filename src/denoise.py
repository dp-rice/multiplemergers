import numpy as np
import pywt

def anscombe_transform(X):
        return 2.0*np.sqrt(X)

def anscombe_transform_inverse(X):
        return (X/2.0)**2

def f_thresh(x, t):
    return np.sign(x)*np.maximum(0, np.abs(x) - t)

def soft_threshold(wavelet_coeffs, t=1.0):
    return [wavelet_coeffs[0]] + [f_thresh(wc,t) for wc in wavelet_coeffs[1:]]

def denoise(X, wavelet='db1', extension_mode='constant', t=1.0, anscombe=True):
    if anscombe:
        X = anscombe_transform(X)
    wc = pywt.wavedec(X, wavelet, extension_mode)
    wc_denoised = soft_threshold(wc, t=t)
    inv = pywt.waverec(wc_denoised, wavelet, extension_mode)
    if len(inv) > len(X):
        inv = inv[:-(len(inv)-len(X))]
    if anscombe:
        inv = anscombe_transform_inverse(inv)
    return inv
