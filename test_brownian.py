# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 21:02:33 2016

@author: William Jones



History:
    22/10/2016: WJ - Created file, replicating matlab script collisiontrial.m

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Define NaN array creation routine to save time
def nanarr(size):
    arr = np.full(size, np.nan)
    return arr


# Defaults for number of balls, ball size and ball mass
n_balls = [20]
size = [0.01]
m_balls = [1]

tot_balls = np.sum(n_balls)
# TODO: generalise, array form with different numbers of different balls
for i in range(0,len(n_balls)-1):
    if i == 0:
        size_arr = np.full(n_balls[0], size[0])
        m_arr = np.full(n_balls[0], m_balls[0])
    else:
        size_arr = np.concatenate((size_arr, np.full(n_balls[i], size[i])))
        m_arr = np.concatenate((m_arr, np.full(n_balls[i], m_balls[i])))
# NOTE: this is a pretty inefficient way of doing things, should inialise an
#   array length of total(n_balls) first and then assign values.

        
# Initialise position and velocity of balls
p = np.random.uniform(low=size[0], high=1-size[0], size=[tot_balls,2])
v = np.random.uniform(low=0, high=1, size=[tot_balls,2])

# Define number of steps (for plotting) and step length, find max time
n_steps = 100
l_step = 0.01
t_max = l_step*n_steps

# Initialise time and sim time variables
t = 0
t_sim = 0

# Initialise indices for last collided balls, for plotting
ii = 0
jj = 0

# Find initial momnetum and energy
momentum = [np.sum((np.sum(v**2, axis=1)**0.5))]
energy = [np.sum((np.sum(v**2, axis=1))/2)]

while t < t_max:
    # Define NaN arrays for wall and ball collision times
    t_wall_x = nanarr([tot_balls])
    t_wall_y = nanarr([tot_balls])
    t_col = nanarr([tot_balls,tot_balls])
    
    for i in range(0,tot_balls-1):
        
        # Find collision times for vertical walls
        temp = np.nanmax([-(p[i,0]-size)/v[i,0],-(p[i,0]-1+size)/v[i,0]])
        if temp >= 0:
            t_wall_x[i] = temp
    
        # Find collision times for horizontal walls
        temp = np.nanmax([-(p[i,1]-size)/v[i,1],-(p[i,1]-1+size)/v[i,1]])
        if temp >= 0:
            t_wall_y[i] = temp
        
        # Find collision times between balls
        if i < tot_balls-1:
            # Loop over all higher index balls to avoid double checking
            for j in range(i+1,tot_balls-1):
                # Calculate quadratic coefficients
                a = np.sum((v[i]-v[j])**2)
                b = 2*np.sum((v[i]-v[j])*(p[i]-p[j]))
                c = np.sum((p[i]-p[j])**2-(size[0]+size[0])**2)
                
                chk = b**2 - 4*a*c
                if chk >= 0:
                    temp = np.nanmin([(-b+chk**0.5)/(2*a),(-b-chk**0.5)/(2*a)])
                    if temp > 1E-10:
                        t_col[i,j] = temp

    # Find minimum times to collision from each
    t_x_min = np.nanmin(t_wall_x)
    t_y_min = np.nanmin(t_wall_y)
    t_col_min = np.nanmin(t_col)

    # Find overall minimum time to collision 
    t_min = np.nanmin([t_x_min,t_y_min,t_col_min])

    # Check if no collision, if so stop
    if np.isnan(t_min):
        sys.exit("No collisions detected")

    # Plot every step up to collison time
    while t+l_step < t_sim+t_min:
        t = t+l_step
        t_diff = t-t_sim
    
        p_temp = p + (v*t_diff)    
        plt.plot([0,1,1,0,0],[0,0,1,1,0],'b-',p_temp[:,0],p_temp[:,1],'bo',p_temp[ii,0],p_temp[ii,1],'ro',p_temp[jj,0],p_temp[jj,1],'ro')
        plt.axis([-0.1,1.1,-0.1,1.1])
        plt.show()
        
        
    #Update sim time
    t_sim = t_sim + t_min
    #Move balls to collision positions
    p = p + v*t_min
    
    if np.sum(t_wall_x == t_min) > 0:
        print("collision x")
        v[t_wall_x == t_x_min,0] = -v[t_wall_x == t_x_min,0]
        wh = np.where(t_wall_x == t_min)
        ii = wh[0][0]

    if np.sum(t_wall_y == t_min) > 0:
        print("collision y")
        v[t_wall_y == t_y_min,1] = -v[t_wall_y == t_y_min,1]
        wh = np.where(t_wall_y == t_min)
        ii = wh[0][0]

    if np.sum(t_col == t_min) > 0:
        print("collision balls")
        wh = np.where(t_col == t_min)
        num_wh = wh[0].size
        for k in range(0,num_wh-1):
            ii = wh[0][i]
            jj = wh[1][i]

            dij = p[ii]-p[jj]
            rij = dij/(np.sum(dij**2)**0.5)
            
            dv = np.sum((v[ii]-v[jj])*rij)*rij
            
            v[ii] = v[ii] + dv
            v[jj] = v[jj] - dv

    momentum = np.concatenate((momentum, [np.sum((np.sum(v**2, axis=1)**0.5))]))
    energy = np.concatenate((energy, [np.sum((np.sum(v**2, axis=1))/2)]))
    
    print("collision",ii,jj)

print(momentum)
print(energy)
x = np.arange(0,momentum.size)
plt.plot(x,momentum,'b-',x,energy,'r-')
plt.show
    

    