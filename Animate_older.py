"""
Created on Sat Oct 22 21:02:33 2016
@authors: William Jones & Luc Moseley  
History:
    22/10/2016: WJ - Created file, replicating matlab script collisiontrial.m
    26/10/2016: LM - fixed array-making bugs, animated properly
    15/11/2016: WJ - Checked animation code working, added momentum based 
                     collisions
    17/11/2016: WJ - Testing generic wall format (defined as a line)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



# Define NaN array creation routine to save time
# This hasn't been used yet, so may not need, but leave for now 
def nanarr(size):
    arr = np.full(size, np.nan)
    return arr

#Test diagnonal wall
test_wall = np.array([[0,0],
                      [1,1]])

#Define a, b coefficients
#test_co = nanarr(1,2)
#Should calculate but just test for now
test_co = np.array([1,0])
    
    
# Defaults for number of balls, ball size and ball mass
n_balls = np.array([15])
sizes = np.array([0.05])
m_balls = np.array([1])

tot_balls = np.sum(n_balls)
# TODO: generalise, array form with different numbers of different balls

for i in range(len(n_balls)):
    if i == 0:
        size_arr = np.full(n_balls[0], sizes[0])
        m_arr = np.full(n_balls[0], m_balls[0])
    else:
        size_arr = np.concatenate((size_arr, np.full(n_balls[i], sizes[i])))
        m_arr = np.concatenate((m_arr, np.full(n_balls[i], m_balls[i])))
# NOTE: this is a pretty inefficient way of doing things, should inialise an
#   array length of total(n_balls) first and then assign values.

#you say this is inefficient, but the way you suggest (below) actually takes
#   longer -- I timed it with %timeit and the original is roughly 33% faster

"""
size_arr = nanarr(tot_balls)
m_arr = nanarr(tot_balls)
for i in range(len(n_balls)):
    if i == 0:
        size_arr[: n_balls[i]] = sizes[i]
        m_arr[: n_balls[i]] = m_balls[i]
    else:
        size_arr[sum(n_balls[:i]): sum(n_balls[:i+1])] = sizes[i]
        m_arr[sum(n_balls[:i]): sum(n_balls[:i+1])] = m_balls[i]
"""
        
# Initialise position and velocity of balls
p = np.random.uniform(low=sizes[0], high=1-sizes[0], size=[tot_balls,2])
v = np.random.uniform(low=-1, high=1, size=[tot_balls,2]) #low changed to -1

# Define number of steps (for plotting) and step length, find max time
n_steps = 100
l_step = 0.01
t_max = l_step*n_steps

# Initialise time
t = 0

# Initialise indices for last collided balls, for plotting
ii = 0
jj = 0

# Find initial momnetum and energy
momentum = [np.sum((np.sum(v**2, axis=1)**0.5))]
energy = [np.sum((np.sum(v**2, axis=1))/2)]

def step():
    global p,t,t_sim
    # Define NaN arrays for wall and ball collision times
    t_wall_x = nanarr([tot_balls])
    t_wall_y = nanarr([tot_balls])
    t_col = nanarr([tot_balls,tot_balls])
    #t_wall_test = nanarr([tot_balls])
    
    #while t < t_max:
        
    for i in range(tot_balls):
        
        # Find collision times for vertical walls
        temp = np.nanmax([-(p[i,0]-sizes[0])/v[i,0],-(p[i,0]-1+sizes[0])/v[i,0]]) 
        #is sizes[0] right?
        #walls set as 0 & 1 here
        if temp >= 0:
            t_wall_x[i] = temp
    
        # Find collision times for horizontal walls
        temp = np.nanmax([-(p[i,1]-sizes[0])/v[i,1],-(p[i,1]-1+sizes[0])/v[i,1]])
        if temp >= 0:
            t_wall_y[i] = temp
        
        # Find collision times between balls
        # Loop over all higher index balls to avoid double checking
        for j in range(i+1,tot_balls):
            # Calculate quadratic coefficients
            a = np.sum((v[i]-v[j])**2)
            b = 2*np.sum((v[i]-v[j])*(p[i]-p[j]))
            c = np.sum((p[i]-p[j])**2-(size_arr[i]+size_arr[j])**2)
            
            chk = b**2 - 4*a*c
            if chk >= 0:
                temp = np.nanmin([(-b+chk**0.5)/(2*a),(-b-chk**0.5)/(2*a)])
                if temp > 1E-10:
                    t_col[i,j] = temp

        
        #Test wall collisions
        #temp = np.nanmin([])
        
        
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
    if l_step < t_min:
        t = t+l_step
        p += v*l_step
            
    else:
        #update time
        t += t_min

        #Move balls to collision positions
        p += v*t_min
            
           
        if np.sum(t_wall_x == t_min) > 0:
            v[t_wall_x == t_x_min,0] = -v[t_wall_x == t_x_min,0]
            ii = np.where(t_wall_x == t_min)
            

        if np.sum(t_wall_y == t_min) > 0:
            v[t_wall_y == t_y_min,1] = -v[t_wall_y == t_y_min,1]
            jj = np.where(t_wall_y == t_min)

        
        if np.sum(t_col == t_min) > 0:
            kk,ll = np.where(t_col == t_min)
            for m in range(0,kk.size):

                dij = p[ll[m]]-p[kk[m]]
                rij = dij/(np.sum(dij**2)**0.5)
                dm = m_arr[kk[m]]/m_arr[ll[m]]
            
                ua = np.sum((v[kk[m]]-v[ll[m]])*rij)
                
                va = ((dm-1)*ua)/(1+dm)
                dva = va - ua
            
                v[kk[m]] = v[kk[m]] + dva*rij
                v[ll[m]] = v[ll[m]] - dm*dva*rij
        
    return p, 
        #momentum = np.concatenate((momentum, [np.sum((np.sum(v**2, axis=1)**0.5))]))
        #energy = np.concatenate((energy, [np.sum((np.sum(v**2, axis=1))/2)]))
    
fig = plt.figure()
# Limits manually entered as 0 & 1 for now, should be whatever box limits we set
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-0.5, 1.5), ylim=(-0.2, 1.2))

# Enter the locations of the particles (size manually =5. should be proportional to size)
particles, = ax.plot([], [], 'bo', ms=5.)
time_text = ax.text(0.02, 0.93, '', transform=ax.transAxes)

#Make box boundary (manually set for now, should automate later)
rect = plt.Rectangle([0,0], 1, 1, lw=2, fc='none')

ax.add_patch(rect)

# Initialize animation
def init():
    global rect
    particles.set_data([], [])
    time_text.set_text('')
    return particles, time_text, rect

# Perform animation step
def animate(i):
    global ax, fig, rect
    
    #step function forward
    step()

    #set marker size based on ball size
    #this size is only appropraite if width of figure = 2 (needs work)
    ms = int(fig.dpi * fig.get_figwidth() * sizes[0])
    particles.set_markersize(ms)
    
    # set position at each step
    particles.set_data(p[:, 0], p[:, 1])
    time_text.set_text('time = %.1f s' % t)
    return particles, time_text,rect

ani = animation.FuncAnimation(fig, animate, frames=100, interval=10, blit=True, init_func=init)
                              
plt.show()

#print(momentum)
#print(energy)
#x = np.arange(0,momentum.size)
#plt.plot(x,momentum,'b-',x,energy,'r-')
#plt.show