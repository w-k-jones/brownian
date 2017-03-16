# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:05:49 2017

@author: pc
"""

import numpy as np

n=10
ang = np.pi/12
length = 20.
ln = length/n

l1 = ln * np.cos(ang)**2
l2 = ln * np.sin(ang)**2 #not actually needed

h = l1 * np.tan(ang)

h_arr = np.arange(2)*h
h_arr = np.tile(h_arr,n)

l_arr = np.arange(n)*ln
l_arr = np.tile(l_arr,[2,1])
l_arr[1] += l1
l_arr = np.reshape(l_arr.T,np.size(l_arr))

out = np.tile(l_arr,[2,1])
out[1] = h_arr
out = out.T

print out