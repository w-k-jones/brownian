# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 13:28:20 2016
@author: William Jones
triindgen - generates two arrays of triangular indices for comparison (i > j)
History:
    22/11/2016: WJ - Created function to return a series of triangular indices 
       (i.e. i > j)
"""

import numpy as np

def tri_indgen(n):
    a = np.arange(n)
    b = a.T
    i = a[a>b]
    j = b[a>b]
    return i,j