# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 10:40:16 2016
Library of common functions for Brownian code
@author: William Jones and Luc Moseley
History:
    23/11/2016: WJ - created, imported functions from animate.py
"""


import numpy as np


"""
Define sub-routines used throghout
"""

"""
NaN array creation routine just to save typing it out
"""
def nanarr(size):
    arr = np.full(size, np.nan)
    return arr


"""
Generates list of i and j indices for lower half matrix indexing (i.e. i > j)
    For upper half triangle switch i and j. Useful for speeding up routines by
    avoiding loops
"""
def tri_indgen(n):
    a = np.arange(n)
    a = np.tile(a,[n,1])
    b = a.T
    i = a[a>b]
    j = b[a>b]
    return i, j


"""
Returns array of times to collision. 
"""
def t_collide(vel, pos, sz_arr):
    n = np.shape(vel)[0]
    [j_arr, i_arr] = tri_indgen(n) #reversed i,j to use upper half triangle
    temp = np.full([np.size(i_arr), 2], np.nan)
    
    a = np.sum((vel[i_arr]-vel[j_arr])**2, axis=1)
    b = 2*np.sum((vel[i_arr]-vel[j_arr])*(pos[i_arr]-pos[j_arr]), axis=1)
    c = np.sum((pos[i_arr]-pos[j_arr])**2, axis=1)-(sz_arr[i_arr]+sz_arr[j_arr])**2
    
    chk = b**2 - 4*a*c
    wh = chk >= 0
    
    temp[wh,0] = (-b[wh]+chk[wh]**0.5)/(2*a[wh])
    temp[wh,1] = (-b[wh]-chk[wh]**0.5)/(2*a[wh])
    temp = np.nanmin(temp, axis=1)
    temp[temp < 1E-10] = np.nan

    return temp


"""
Returns array of distances to collisions
"""
def get_dist(pos,sz_arr):
    n = np.shape(pos)[0]
    [j_arr, i_arr] = tri_indgen(n) #reversed i,j to use upper half triangle
    
    temp = np.sum((pos[i_arr]-pos[j_arr])**2, axis=1)**0.5 - (sz_arr[i_arr]+sz_arr[j_arr])
    return temp
