# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 16:11:47 2017

@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt

def mb_2d(p_in,T,m):
    return ((-2*T/m) * np.log(1 - (m/T)**2 * p_in))**0.5
    
p = np.arange(0,1,0.01)

v = mb_2d(p,1,1)

plt.plot(p,v)
plt.show