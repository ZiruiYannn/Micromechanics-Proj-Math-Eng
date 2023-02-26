# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:38:44 2023

@author: Andreas
"""

import numpy as np



def strain(E, mat, c, tol):
    dims = np.shape(mat)
    l = dims[0]; m = dims[1]; n = dims[2]
    
    eps = np.zeros((l, m, n, 3, 3))
    for i in range(l):
        for j in range(m):
            for k in range(n):
                eps[i, j, k, :, :] = E
    sig = np.zeros((l, m, n, 3, 3))
    tau = np.zeros((3, 3))
    c0, lambd0, mu0 = refMat(c)
    
    converged = False
    while not converged:
        # main loop
        return
    
    return eps



# Determine reference material and corresponding Lamé coefficients
def refMat(c):
    return



# Perform the convergence test
def convTest(sig, i, j, k, dims, tol):
    return



# Determine the Greenoperator
def greenOp(lambd0, mu0, i, j, k, dims):
    return



# Determine the wave vector
def waveVec(i, j, k, dims):
    return










