"""
Created on Wed Jan 25 09:54:15 2017
@authors: William Jones & Luc Moseley  
History:
    25/01/2017: LM - Set up separate file with vars imported from System
    08/02/2017: LM - Fixed bugs & fully assimilated with new code
Variable naming convention:
    w_... Wall related variable (should probably make class)
    n_... # of ...
    p_... position of ...
    v_... velocity of ...
    ..._n normal vector of ...
    ..._p parallel vector of ...
"""

"""
You cannot run this twice in same shell window, because the original 
step function updates global variables - this file simply takes global 
variables as they were at end of previous simulation when run again.
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from brownian_tools import *
import brownian_classes as brw
import brownian_system as brs

#set up Bernoulli tube system
in_co = np.array([[0,0],
                  [0.2,0],
                  [0.4,0.2],
                  [0.6,0.2],
                  [0.8,0],
                  [1,0],
                  [1,1],
                  [0.8,1],
                  [0.6,0.8],
                  [0.4,0.8],
                  [0.2,1],
                  [0,1]])

#co-ords for opposing sawtooth
in_co = np.array([[0,0],
                  [0.2,0.2],
                  [0.2,0],
                  [0.4,0.2],
                  [0.4,0],
                  [0.6,0.2],
                  [0.6,0],
                  [0.8,0.2],
                  [0.8,0],
                  [1,0.2],
                  [1,1],
                  [0.8,0.8],
                  [0.8,1],
                  [0.6,0.8],
                  [0.6,1],
                  [0.4,0.8],
                  [0.4,1],
                  [0.2,0.8],
                  [0.2,1],
                  [0,0.8]])

wal = brw.wall_shape(in_co)
bal = brw.balls([100,10],[0.01,0.05],[1,25],2,wal)
inst = brs.system(wal,bal)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
xlim=(wal.xlim[0]-1, wal.xlim[1]+1), ylim=(wal.ylim[0]-1, wal.ylim[1]+1))

# Enter the locations of the particles
particles, = ax.plot([], [], 'bo', ms=5.)
if len(bal.r2) > 1:
    particles1, = ax.plot([], [], 'bo', ms=10.)
if len(bal.r2) > 2:
    particles2, = ax.plot([], [], 'bo', ms=15.)
time_text = ax.text(0.025, 0.93, '', transform=ax.transAxes)
energy_text = ax.text(0.025, 0.88, '', transform=ax.transAxes)
pressure_text = ax.text(0.025, 0.83, '', transform=ax.transAxes)
avg_press_text = ax.text(0.025, 0.78, '', transform=ax.transAxes)

#Make container boundary
ax.plot(wal.co_plt[:,0],wal.co_plt[:,1])

#ax.add_patch(rect)

# Initialize animation
def init():
    #global rect
    particles.set_data([], [])
    if len(bal.r2) > 1:
        particles1.set_data([], [])
    if len(bal.r2) > 2:
        particles2.set_data([], [])
    
    time_text.set_text('')
    energy_text.set_text('')
    pressure_text.set_text('')
    avg_press_text.set_text('')
    
    if len(bal.r2) == 1:
        return particles, time_text, energy_text, pressure_text, avg_press_text #rect
    elif len(bal.r2) == 2:
        return particles, particles1, time_text, energy_text, pressure_text, avg_press_text #rect
    else:
        return particles, particles1, particles2, time_text, energy_text, pressure_text, avg_press_text #rect

# Perform animation step
def animate(i):
    global ax, fig #rect
    
    #step function forward
    #p, energy, pressure, avg_press = inst.step()
    p, t = inst.step()
    energy = 0
    pressure = 0
    avg_press = 0

    #set marker size based on ball size
    ms = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*bal.r2[0]
    particles.set_markersize(ms)
    particles.set_data(p[0:bal.n_balls[0], 0], p[0:bal.n_balls[0], 1])
    if len(bal.r2) > 1:
        ms1 = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*bal.r2[1] #this abs could be the problem
        particles1.set_markersize(ms1)
        particles1.set_markerfacecolor('r')
        particles1.set_data(p[bal.n_balls[0]:sum(bal.n_balls[0:2]), 0], p[bal.n_balls[0]:sum(bal.n_balls[0:2]), 1])
    if len(bal.r2) > 2:
        ms2 = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*bal.r2[2]
        particles2.set_markersize(ms2)
        particles2.set_data(p[sum(bal.n_balls[0:2]):sum(bal.n_balls[0:3]), 0], p[sum(bal.n_balls[0:2]):sum(bal.n_balls[0:3]), 1])
        
    time_text.set_text('Time = %.1f s' % t)
    energy_text.set_text('Energy = %.2f J' % energy)
    pressure_text.set_text('Pressure = %.2f mPa' % (pressure*1000.))
    avg_press_text.set_text('Average Pressure = %.1f mPa' % (avg_press*1000.))
    
    if len(bal.r2) == 1:
        return particles, time_text, energy_text, pressure_text, avg_press_text #rect
    elif len(bal.r2) == 2:
        return particles, particles1, time_text, energy_text, pressure_text, avg_press_text #rect
    else:
        return particles, particles1, particles2, time_text, energy_text, pressure_text, avg_press_text #rect

ani = animation.FuncAnimation(fig, animate, frames=600, interval=10, blit=True, init_func=init)
                              
plt.show()