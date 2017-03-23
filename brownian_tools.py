# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 10:40:16 2016
Library of common functions for Brownian code
@author: William Jones and Luc Moseley
History:
    23/11/2016: WJ - created, imported functions from animate.py
    23/11/2016: WJ - started work on t_wall; new general wall collision code
    24/11/2016: WJ - Added in new collision routine for walls, incomplete, 
        needs vectorising
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
    
def rand_mb2d(T,m):
    p_in = np.random.uniform()
    return ((-2*T/m) * np.log(1 - (m/T)**2 * p_in))**0.5


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
def t_collide(pos, vel, sz_arr):
    n = np.shape(vel)[0]
    nd = np.shape(vel)[1]
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
Returns array of collision times to wall
"""
def t_wall(pos, vel, sz_arr, coord):
    wall_vec = coord[1]-coord[0]
    # Calculate normalised normal vector to wall
    wall_n = np.dot(np.array([[0,-1],[1,0]]), wall_vec)
    wall_n = wall_n/np.sum(wall_n**2)**0.5
    
    #Position vector to wall
    pos_n = np.dot(coord[0]-pos, wall_n)
    vel_n = np.dot(vel, wall_n)
    
    #minimum time to collision
    temp = np.nanmin([(pos_n-sz_arr)/vel_n, (pos_n+sz_arr)/vel_n])
    if temp < 1E-10:
        temp = np.nan
        
    #TODO: check for end of walls, important for obtuse shapes
    
    return temp
    
    

"""
Returns array of distances to collisions
"""
def get_dist(pos,sz_arr):
    n = np.shape(pos)[0]
    [j_arr, i_arr] = tri_indgen(n) #reversed i,j to use upper half triangle
    
    temp = np.sum((pos[i_arr]-pos[j_arr])**2, axis=1)**0.5 - (sz_arr[i_arr]+sz_arr[j_arr])
    return temp

"""
Finds vector normal for sawtooth (to create trigger wall on correct side)
"""

def normal(x_v2, y_v2, spike_side):
    x_v3, y_v3 = y_v2, x_v2
    if spike_side == 1:
        y_v3 = y_v3 * -1
    elif spike_side == -1:
        x_v3 = x_v3 * -1
    return np.array([x_v3[1]-x_v3[0], y_v3[1]-y_v3[0]])
