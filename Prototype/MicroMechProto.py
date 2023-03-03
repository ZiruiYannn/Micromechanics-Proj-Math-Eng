# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:38:44 2023

@author: Andreas
"""

import numpy as np



# Basic scheme
def strain(E, mat, c, tol):
    dims = np.shape(mat)
    # TODO: determine period
    
    eps = np.zeros((dims[0], dims[1], dims[2], 6))
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                eps[i, j, k, :, :] = E
    sig = np.zeros((6, 1))
    tau = np.zeros((6, 1))
    eps_freq = np.zeros((3, 3))
    sig_freq = np.zeros((3, 3))
    tau_freq = np.zeros((3, 3))
    c0 = refMat(c)
    e = np.inf
    sig_freq_0 = 0
    
    while e > tol:
        e = 0
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    sig = stress(eps[i, j, k, :], c[mat[i, j, k], 0], c[mat[i, j, k], 1])
                    tau = sig - stress(stress(eps[i, j, k, :], c0[0], c0[1]))
                    sig_freq = FFT(sig)
                    tau_freq = FFT(tau)
                    if i == 0 and j == 0 and k == 0:
                        eps_freq = vec2ten(E)
                        sig_freq_0 = np.linalg.norm(sig_freq)**2 
                        # TODO: compute local frequency norm and initialize error
                    else:
                        eps_freq = np.tensordot(-greenOp(c0[0], c0[1], [i, j, k], dims), tau_freq)
                        # TODO: compute local frequency norm and add to error
                    eps[i, j, k, :] = IFFT(eps_freq)
        e /= dims[0]*dims[1]*dims[2]
        e = np.sqrt(e) 
        e /= sig_freq_0
                        
    return eps



# Determine reference material and corresponding Lamé coefficients
def refMat(c):
    c0 = np.zeros((1, 2))
    c0[0] = (max(c[:, 0]) + min(c[:, 0])) / 2
    c0[1] = (max(c[:, 1]) + min(c[:, 1])) / 2
    
    return c0



# Determine the stress for a given strain and material
# in 2nd-order tensor notation: sig = 2*mu*eps + lambda*trace(eps)*I_3
def stress(strain, lambd, mu):
    I = np.zeros((6, 1))
    I[1:3] = 1
    
    return 2*mu*strain + lambd*sum(strain[1:3])*I



# Determine Fourier transform
def FFT(real):
    # TODO: implement this
    
    return



# Determine inverse Fourier transform
def IFFT(freq):
    # TODO: implement this
    
    return



# Determine the 2nd-order tensor equivalent for a given vector
def vec2ten(v):
    T = np.diag(v[1:3])
    T[2, 3] = T[3, 2] = v[4]
    T[1, 3] = T[3, 1] = v[5]
    T[1, 2] = T[2, 1] = v[6]
    
    return T



# Determine the vector equivalent for a given 2nd-order tensor
def ten2vec(T):
    v = np.zeros((6, 1))
    v[1:3] = np.diag(T)
    v[4] = T[2, 3]
    v[5] = T[1, 3]
    v[6] = T[1, 2]
    
    return v



# Determine the Greenoperator
def greenOp(lambd0, mu0, inds, dims, prds):
    xi = waveVec(inds, dims, prds)
    gam = np.zeros((3,3,3,3))
    for k in range(3):
        for h in range(3):
            for i in range(3):
                for j in range(3):
                    if k == i:
                        gam[k, h, i, j] += xi[h]*xi[j]
                    if h == i:
                        gam[k, h, i, j] += xi[k]*xi[j]
                    if k == j:
                        gam[k, h, i, j] += xi[h]*xi[i]
                    if h == j:
                        gam[k, h, i, j] += xi[k]*xi[i]
                    gam[k, h, i, j] /= 4*mu0*np.linalg.norm(xi)**2 # TODO: check if this is the correct interpretation
                    gam[k, h, i, j] -= (lambd0+mu0)/(mu0*(lambd0+2*mu0)) * xi[i]*xi[j]*xi[k]*xi[h]/np.linalg.norm(xi)**4 # TODO: check if this is the correct interpretation 
                    
    return gam



# Determine the wave vector
def waveVec(inds, dims, prds):
    xi = np.zeros(3)
    for i in range(3):
        if dims[i]%2 == 0:
            xi[i] = (-dims[i]/2 + inds[i]+1) / prds[i]
        else:
            xi[i] = (-(dims[i]-1)/2 + inds[i]) / prds[i]
            
    return xi










