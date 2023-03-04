# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:38:44 2023

@author: Andreas
"""

import numpy as np



# Basic scheme
def strain(E, mat, c, prds, tol):
    dims = np.shape(mat)
    dims_freq = dims
    dims_freq[2] = dims_freq[2]//2 + 1 # (i)rfftn reduces the last transformation axis
    
    eps = np.zeros((dims[0], dims[1], dims[2], 6))
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                eps[i, j, k, :, :] = E
    sig = np.zeros((dims[0], dims[1], dims[2], 6))
    tau = np.zeros((dims[0], dims[1], dims[2], 6))
    eps_freq = np.zeros((dims_freq[0], dims_freq[1], dims_freq[2], 6)) # not necessary for FFTW
    tau_freq = np.zeros((dims_freq[0], dims_freq[1], dims_freq[2], 6)) # not necessary for FFTW
    gam = np.zeros((3, 3, 3, 3))
    c0 = refMat(c)
    e = np.inf
    
    while e > tol:
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    sig[i, j, k, :] = stress(eps[i, j, k, :], c[mat[i, j, k], 0], c[mat[i, j, k], 1])
                    tau[i, j, k, :] = sig - stress(stress(eps[i, j, k, :], c0[0], c0[1]))
        
        e = error(sig, dims_freq, prds)
        tau_freq = np.fft.rfftn(tau, axes=(0, 1, 2)) # FFTW does this in-place
        
        for i in range(dims_freq[0]):
            for j in range(dims_freq[1]):
                for k in range(dims_freq[2]):
                    gam = greenOp(c0[0], c0[1], [i, j, k], dims_freq, prds)
                    eps_freq = ten2vec(np.tensordot(-gam, vec2ten(tau_freq[i, j, k, :]))) # FFT conserves symmetry            
        eps_freq[0, 0, 0, :] = E
        eps = np.fft.irfftn(eps_freq, s=dims, axes=(0, 1, 2)) # FFTW does this in-place
                        
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
    I = np.zeros(6)
    I[1:3] = 1
    
    return 2*mu*strain + lambd*sum(strain[0:2])*I



# Determine the error
def error(sig, dims_freq, prds):
    sig_freq = np.fft.rfftn(sig, axes=(0, 1, 2))  # FFTW does this in-place
    e = 0
    for i in range(dims_freq[0]):
        for j in range(dims_freq[1]):
            for k in range(dims_freq[2]):
                xi = waveVec([i, j, k], dims_freq, prds)
                e += np.linalg.norm(np.matmul(vec2ten(sig_freq[i, j, k, :]), xi))**2 # FFT conserves symmetry
    e /= dims_freq[0]*dims_freq[1]*dims_freq[2]
    e = np.sqrt(e)
    e /= np.linalg.norm(vec2ten(sig_freq[i, j, k, :]))**2
    
    return



# Determine the 2nd-order tensor equivalent for a given vector
def vec2ten(v):
    T = np.diag(v[0:2])
    T[2, 3] = T[3, 2] = v[4]
    T[1, 3] = T[3, 1] = v[5]
    T[1, 2] = T[2, 1] = v[6]
    
    return T



# Determine the vector equivalent for a given 2nd-order tensor
def ten2vec(T):
    v = np.zeros(6)
    v[0:2] = np.diag(T)
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
                    gam[k, h, i, j] /= 4*mu0*np.linalg.norm(xi)**2 
                    gam[k, h, i, j] -= (lambd0+mu0)/(mu0*(lambd0+2*mu0)) * xi[i]*xi[j]*xi[k]*xi[h]/np.linalg.norm(xi)**4  
                    
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





def main():
    # coefficients for steel
    nu = 0.28 # Poisson ratio
    mod = 210 # Elasticity modulus
    lambd = mod*nu / ((1+nu)*(1-2*nu))
    mu = mod / (2*(1+nu))
    
    E0 = 0.01
    E = np.array([E0, -nu*E0, -nu*E0, 0.0, 0.0, 0.0]) 
    mat = np.zeros(5, 5, 5)
    c = np.array([[lambd, mu]]) 
    prds = np.array([5, 5, 5]) 
    tol = 1e-4
    
    eps = strain(E, mat, c, prds, tol) # sig = E0*mod
    print(eps)
    
    return





main()







