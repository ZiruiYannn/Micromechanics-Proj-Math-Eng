# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:38:44 2023

@author: Andreas
"""

import numpy as np



# Basic scheme
def strain(E, mat, c, prds, tol=1e-4, maxit=50):
    dims = np.shape(mat)
    dims_freq = np.array(dims)
    dims_freq[2] = dims_freq[2]//2 + 1 # rfftn reduces the last transformation axis dimension
    
    eps = np.zeros((dims[0], dims[1], dims[2], 6))
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                eps[i, j, k, :] = E
    sig = np.zeros((dims[0], dims[1], dims[2], 6))
    tau = np.zeros((dims[0], dims[1], dims[2], 6))
    eps_freq = np.zeros((dims_freq[0], dims_freq[1], dims_freq[2], 6), dtype=complex) # not necessary for FFTW
    tau_freq = np.zeros((dims_freq[0], dims_freq[1], dims_freq[2], 6), dtype=complex) # not necessary for FFTW
    xi = np.zeros(3)
    gam = np.zeros((3, 3, 3, 3))
    c0 = refMat(c)
    e = np.inf
    
    it = 0
    while e > tol and it < maxit:
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    sig[i, j, k, :] = stress(eps[i, j, k, :], c[mat[i, j, k], 0], c[mat[i, j, k], 1])
                    tau[i, j, k, :] = sig[i, j, k, :] - stress(eps[i, j, k, :], c0[0], c0[1])
        
        e = error(sig, dims_freq, prds)
        print(e)
        tau_freq = np.fft.rfftn(tau, axes=(0, 1, 2)) # FFTW does this in-place
        
        for i in range(dims_freq[0]):
            for j in range(dims_freq[1]):
                for k in range(dims_freq[2]):
                    xi = waveVec([i, j, k], dims_freq, prds)
                    if np.array_equal(xi, [0.0, 0.0, 0.0]):
                        eps_freq[i, j, k, :] = dims[0]*dims[1]*dims[2]*E
                    else:
                        gam = greenOp(c0[0], c0[1], xi)
                        eps_freq[i, j, k, :] = ten2vec(np.tensordot(-gam, vec2ten(tau_freq[i, j, k, :]))) # FFT conserves symmetry            
        eps = np.fft.irfftn(eps_freq, s=dims, axes=(0, 1, 2)) # FFTW does this in-place
        it += 1
        
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                sig[i, j, k, :] = stress(eps[i, j, k, :], c[mat[i, j, k], 0], c[mat[i, j, k], 1])
                        
    return eps, sig, e



# Determine reference material and corresponding Lamé coefficients
def refMat(c):
    c0 = np.zeros(2)
    c0[0] = (max(c[:, 0]) + min(c[:, 0])) / 2
    c0[1] = (max(c[:, 1]) + min(c[:, 1])) / 2
    
    return c0



# Determine the stress for a given strain and material
# in 2nd-order tensor notation: sig = 2*mu*eps + lambda*trace(eps)*I_3
def stress(strain, lambd, mu):
    I = np.zeros(6)
    I[0:3] = 1
    
    return 2*mu*strain + lambd*sum(strain[0:3])*I



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
    e /= np.linalg.norm(vec2ten(sig_freq[0, 0, 0, :]))**2
    
    return e



# Determine the wave vector
def waveVec(inds, dims, prds):
    xi = np.zeros(3)
    for i in range(3):
        if dims[i]%2 == 0:
            bound = dims[i]//2-1
        else:
            bound = dims[i]//2
            
        if inds[i] <= bound:
            xi[i] = inds[i] / prds[i]
        else:
            xi[i] = (-dims[i] + inds[i]) / prds[i]
            
    return xi



# Determine the Greenoperator
def greenOp(lambd0, mu0, xi):
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



# Determine the 2nd-order tensor equivalent for a given vector
def vec2ten(v):
    T = np.diag(v[0:3])
    T[1, 2] = T[2, 1] = v[3]
    T[0, 2] = T[2, 0] = v[4]
    T[0, 1] = T[1, 0] = v[5]
    
    return T



# Determine the vector equivalent for a given 2nd-order tensor
def ten2vec(T):
    v = np.zeros(6, dtype=T.dtype)
    v[0:3] = np.diag(T)
    v[3] = T[1, 2]
    v[4] = T[0, 2]
    v[5] = T[0, 1]
    
    return v





def main():
    # coefficients for steel
    nu = 0.28 # Poisson ratio
    mod1 = 210e9 # Elasticity modulus
    lambd1 = mod1*nu / ((1+nu)*(1-2*nu))
    mu1 = mod1 / (2*(1+nu))
    
    # altered coefficients 2
    mod2 = 0.4*mod1
    lambd2 = mod2*nu / ((1+nu)*(1-2*nu))
    mu2 = mod2 / (2*(1+nu))
    
    # altered coefficients 3
    mod3 = 0.6*mod1
    lambd3 = mod3*nu / ((1+nu)*(1-2*nu))
    mu3 = mod3 / (2*(1+nu))
    
    E0 = 0.01
    E = np.array([-nu*E0, -nu*E0, E0, 0.0, 0.0, 0.0]) 
    mat = np.zeros((10, 10, 60), dtype=int)
    for i in range(20, 40):
        mat[:, :, i] = 1
    for i in range(40, 60):
        mat[:, :, i] = 2
    c = np.array([[lambd1, mu1], [lambd2, mu2], [lambd3, mu3]]) 
    prds = np.array([10, 10, 50]) 
    
    eps, sig, e = strain(E, mat, c, prds, -np.inf, 30) 

    # # equivalent coefficients
    # mod_eq = mod1*mod2 / (mod1+mod2)
    # lambd_eq = mod_eq*nu / ((1+nu)*(1-2*nu))
    # mu_eq = mod_eq / (2*(1+nu))

    # mat_eq = np.zeros((3, 3, 10), dtype=int)
    # c_eq = np.array([[lambd_eq, mu_eq]]) 
    
    # eps_eq, sig_eq, e_eq = strain(E, mat_eq, c_eq, prds)
    
    # avg1 = np.average(sig[:, :, 0:20, 0])
    # avg2 = np.average(sig[:, :, 20:40, 0])
    # avg3 = np.average(sig[:, :, 40:, 0])
    # print(avg3/avg1)
    
    avg1 = np.average(eps[:, :, 0:20, 2])
    avg2 = np.average(eps[:, :, 20:40, 2])
    avg3 = np.average(eps[:, :, 40:, 2])
    print(avg1/avg2)
    
    return





main()







