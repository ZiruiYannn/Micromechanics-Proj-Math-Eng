# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 14:57:08 2023

@author: Zirui
"""

import numpy as np
import pyfftw
from functions import fft44dtensor, ifft44dtensor


#----------------
'''
This way of storing the data is not sufficient enough. It will be better if I use the last axes
to indicate ij of strain_ij 
'''
#----------------

class micromechanics:
    
    def __init__(self, E, cuboid_dim, material_notation, array_lamda_mu, period, tol, maxit):
        #-------------------
        '''
            E is the initial strain array
            cuboid_dim is a tuple
            material_notation is an array which store the number of material for each point
            array_lamda_mu is a 2 by n array. n^th columne store the lamda and mu of the n^th material. 
            tol tolerance for convergence test
            maxit maximun iteration
        '''
        #--------------------
        self.dim = cuboid_dim
        self.strain = self.init_strain(E, self.dim)
        self.stress = np.zeros((6,)+self.dim, dtype='complex128')
        self.polarization = np.empty((6,)+self.dim, dtype='complex128') #tau
        self.polarization_freq = np.empty((6,)+self.dim, dtype='complex128')
        self.stress_freq = np.zeros((6,)+self.dim, dtype='complex128')
        self.strain_freq = np.empty((6,)+self.dim, dtype='complex128')
        
        
        self.material = material_notation
        self.lamda = array_lamda_mu[0,:]
        self.mu = array_lamda_mu[1,:]
        '''
        haven't checked self.lamda self.mu are of the same size
        '''
        self.num_material = np.size(self.lamda)
        self.lamda_ref = (max(self.lamda) + min(self.lamda))/2
        self.mu_ref = (max(self.mu) + min(self.mu))/2
        
        self.c=np.empty(self.num_material, 6, 6)
        for i in range(self.num_material):
            self.c[i,:,:]=self.lame2coeff(self.lamda[i],self.mu[i])
        
        self.c_ref = self.lame2coeff(self.lamda_ref, self.mu_ref)
        
        self.period = period
        self.tol = tol
        self.maxit = maxit
        
    def init_strain(_E, _dim):
        #------------------------------
        '''
            E is the initial strain array
            dim is the shape of cuboid
            return a 4d array. the length of first dimension is 6, which indicates that we use a 6-entry 
            array to represent the strain tensor of each point. The other three dimension gives the position
            of each point.
        '''
        #------------------------------
        _strain = np.empty((6,)+_dim, dtype='complex128')
        for i in range(6):
            _strain[i,:,:,:] = _E[i]*np.ones(_dim) + 1j*np.zeros(_dim)
        
        return _strain
        
        
    def lame2coeff(_lamda, _mu):
        c_xx = np.array([2*_mu+_lamda, _lamda, _lamda, 0, 0, 0])
        c_yy = np.array([_lamda, _lamda, 2*_mu+_lamda, 0, 0, 0])
        c_zz = np.array([_lamda, _lamda, 2*_mu+_lamda, 0, 0, 0])
        c_xy = np.array([0, 0, 0, 2*_mu, 0, 0 ])
        c_xz = np.array([0, 0, 0, 0, 2*_mu, 0 ])
        c_yz = np.array([0, 0, 0, 0, 0, 2*_mu ])
        
        _c = np.array([c_xx, c_yy, c_zz, c_xy, c_xz, c_yz])
        
        return _c
        
    def c_colon_strain(self, _strain):
        #-----------------------------
        '''
        compute stress at each position 
        
        #np.dot(a,b) if a is N-d and b is 1-d, it is a sum product over the last axis of a and b
        '''
        #-----------------------------
        _stress = np.empty((6,)+self.dim, dtype='complex128')
        a = self.dim[0]
        b = self.dim[1]
        c = self.dim[2]
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    n = self.material[b*c*i + c*j + k] #the number of the material at (i,j,k) 
                    _c = self.c[n,:,:]
                    temp = np.dot(_c, _strain[:,i,j,k])
                    _stress[:,i,j,k] = temp
                    
        return _stress
    
    def c_colon_strain_ref(self, _strain):
        #-----------------------------
        '''
        compute reference stress at each position 
        
        
        '''
        #-----------------------------
        _stress = np.empty((6,)+self.dim, dtype='complex128')
        a = self.dim[0]
        b = self.dim[1]
        c = self.dim[2]
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    temp = np.dot(self.c_ref, _strain[:,i,j,k])
                    _stress[:,i,j,k] = temp
                    
        return _stress
        
    def error(self):
        a = self.dim[0]
        b = self.dim[1]
        c = self.dim[2]
        _err = 0
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    _wave_vec = self.compute_wave_vec([i,j,k])
                    _err += np.linalg.norm(np.dot(vec2ten(self.stress_freq[:,i,j,k]), _wave_vec))**2  
        _err /= a*b*c
        _err = np.sqrt(e)
        _err /= np.linalg.norm(self.stress_freq[:, a//2, b//2, c//2])**2
        
        return _err
    
    def compute_wave_vec(self, _index):
        #-----------------------------
        '''
        compute the wave vector
        
        the index in py is from 0 to n-1, but the paper use 1 to n        
        '''
        #-----------------------------
        _wave_vec = np.empty(3)
        for i in range(3):
            n = self.dim[i]
            _wave_vec[i] = (-(n//2) + _index[i])/self.period[i]
            
        return _wave_vec
    
    def vec2ten(v):
        #-----------------------------
        '''
        6-entry vector to 3 by 3tensor     
        '''
        #-----------------------------
        ten = np.diag(v[0:2])
        ten[0,1] = ten[1,0] = v[3]
        ten[0,2] = ten[2,0] = v[4]
        ten[2,1] = ten[1,2] = v[5]
        
    def green_oper(self, _index):
        #-----------------------------
        '''
        return a 6 by 6 array (Green operator)    
        '''
        #-----------------------------
        _green = np.empty((6,6))
        
        _wave_vec = self.compute_wave_vec(_index)
        norm_wave_vec = np.linalg.norm(_wave_vec)
        coef_a = 1/4/self.mu_ref
        coef_a /= norm_wave_vec**2
        coef_b = (self.lamda_ref + self.mu_ref)/(self.mu_ref * (self.lamda_ref + 2 * self.mu_ref))
        coef_b /= norm_wave_vec**4
        
        mat = np.dot(_index.T, _index)
        
        _green[0, :] = np.array([coef_a * mat[0,0] * 4 -coef_b * mat[0,0]**2,
                                 0,
                                 0,
                                 coef_a * mat[0,1] * 2,
                                 coef_a * mat[0,2] * 2,
                                 0])
        _green[1, :] = np.array([0,
                                 coef_a * mat[1,1] * 4 -coef_b * mat[1,1]**2,
                                 0,
                                 coef_a * mat[0,1] * 2,
                                 0,
                                 coef_a * mat[1,2] * 2])
        _green[2, :] = np.array([0,
                                 0,
                                 coef_a * mat[2,2] * 4 -coef_b * mat[2,2]**2,
                                 0,
                                 coef_a * mat[0,2] * 2,
                                 coef_a * mat[1,2] * 2])
        _green[3, :] = np.array([coef_a * mat[0,1] * 2,
                                 coef_a * mat[0,1] * 2,
                                 0,
                                 coef_a * (mat[0,0] + mat[1,1]) -coef_b * mat[0,1]**2,
                                 coef_a * mat[1,2],
                                 coef_a * mat[0,2]])
        _green[4, :] = np.array([coef_a * mat[0,2] * 2,
                                 0,
                                 coef_a * mat[0,2] * 2,
                                 coef_a * mat[1,2],
                                 coef_a * (mat[0,0] + mat[2,2]) -coef_b * mat[0,2]**2,
                                 coef_a * mat[0,1]])
        _green[5, :] = np.array([0,
                                 coef_a * mat[1,2] * 2,
                                 coef_a * mat[1,2] * 2,
                                 coef_a * mat[0,2],
                                 coef_a * mat[0,1],
                                 coef_a * (mat[1,1] + mat[2,2]) -coef_b * mat[1,2]**2])
        
        return _green
        
    def step_d(self):
        #-----------------------------
        '''
        the whole step d in the algorithm
        modify self.strain_freq
        '''
        #-----------------------------
        self.strain_freq = np.zeros((6,)+self.dim, dtype='complex128')
        a = self.dim[0]
        b = self.dim[1]
        c = self.dim[2]
        
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    green = self.green_oper([i,j,k])
                    temp = - np.dot(green, self.polarization_freq[:,i,j,k])
                    self.strain_freq[:,i,j,k] = temp
        
        
    def iteration(self):
        self.stress = self.c_colon_strain(self.strain)
        self.stress_freq = fft44dtensor(self.stress)
        
        err = self.error
        iter_count = 0
        
        while err < self.tol and iter_count < self.maxit:
            self.polarization = self.stress - self.c_colon_strain(self.strain)            
            self.polarization_freq = fft44dtensor (self.polarization)
            self.step_d()
            self.strain = ifft44dtensor(self.strain_freq)
            self.stress = self.c_colon_strain(self.strain)
            self.stress_freq = fft44dtensor(self.stress)
            
            err = self.error
            iter_count +=1
            
            
            
        
    
        
        
        
        
        
        
        
        
        
        
        