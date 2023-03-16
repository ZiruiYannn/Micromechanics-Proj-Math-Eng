# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 01:34:55 2023

@author: Hirvio
"""

import numpy as np
import pyfftw

def fft44dtensor(_tensor): 
    #-------------------------------
    '''
    _tensor is a 4dimensional tensor and the length of the first axes of it is 6
    _dim here is different from cuboid_dim. _dim = ((6,)+cuboid_dim) = np.shape(_tensor) 
    '''
    #--------------------------------
    _dim = np.shape(_tensor)
    a = pyfftw.empty_aligned(_dim, dtype='complex128')
    b = pyfftw.empty_aligned(_dim, dtype='complex128')
     
    a = _tensor
    fft_object_a = pyfftw.FFTW(a, b, axis=(1,2,3))
    fft_object_a()
    
    return b

def ifft44dtensor(_tensor): 
    #-------------------------------
    '''
    _tensor is a 4dimensional tensor and the length of the first axes of it is 6
    _dim here is different from cuboid_dim. _dim = ((6,)+cuboid_dim) = np.shape(_tensor) 
    '''
    #--------------------------------
    _dim = np.shape(_tensor)
    a = pyfftw.empty_aligned(_dim, dtype='complex128')
    b = pyfftw.empty_aligned(_dim, dtype='complex128')
     
    a = _tensor
    ifft_object_a = pyfftw.FFTW(a, b, direction='FFTW_BACKWARD', axis=(1,2,3))
    ifft_object_a()
    
    return b



