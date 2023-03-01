# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:38:44 2023

@author: Andreas
"""

import numpy as np



# Basic scheme
def strain(E, mat, c, tol):
    dims = np.shape(mat)
    # determine period
    
    eps = np.zeros((dims[0], dims[1], dims[2], 3, 3))
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                eps[i, j, k, :, :] = E
    sig = np.zeros((3, 3))
    tau = np.zeros((3, 3))
    c0, lambd0, mu0 = refMat(c)
    e = np.inf
    s = 0
    
    while e > tol:
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    sig = np.tensordot(c[mat[i, j, k]], eps[i, j, k, :, :])
                    tau = sig - np.tensordot(c0, eps[i, j, k, :, :])
                    # TODO: perform FFT on sig 'in place'
                    # TODO: perform FFT on tau 'in place'
                    if i == 0 and j == 0 and k == 0:
                        eps[i, j, k, :, :] = E
                        s = np.linalg.norm(sig)**2
                        # TODO: compute local frequency norm and initialize error
                    else:
                        eps[i, j, k, :, :] = np.tensordot(-greenOp(lambd0, mu0, [i, j, k], dims), tau)
                        # TODO: compute local frequency norm and add to error
                    # TODO: perform inverse FFT of eps[i, j, k] 'in place'
        e /= dims[0]*dims[1]*dims[2]
        e = np.sqrt(e) 
        e /= s
                        
    return eps



# Determine reference material and corresponding Lamé coefficients
def refMat(c):
    # TODO: implement
    return



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
                    gam[k, h, i, j] /= 4*mu0*np.linalg.norm(xi)**2 # TODO: correct interpretation?
                    gam[k, h, i, j] -= (lambd0+mu0)/(mu0*(lambd0+2*mu0)) * xi[i]*xi[j]*xi[k]*xi[h]/np.linalg.norm(xi)**4 # TODO: correct interpretation?
                    
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










