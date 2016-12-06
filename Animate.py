"""
Created on Sat Oct 22 21:02:33 2016
@authors: William Jones & Luc Moseley  
History:
    22/10/2016: WJ - Created file, replicating matlab script collisiontrial.m
    26/10/2016: LM - Fixed array-making bugs, animated properly
    1/11/2016: WJ - Added collision code back in, currently not working properly
    3/11/2016: WJ - Fixed looping code - would collide regardless after each 
        loop. Added momentum to collisions so different masses work. Added 
        size_arr, collisions working but not animation.
    8/11/2016: LM - Made marker size vary according to ball size Balls are 
        overlapping walls and sometimes each other
    9/11/2016: LM - Made box class and walls use this class, added energy check
    10/11/2016: LM - Put in pressure read for animation
    22/11/2016: WJ - Updated collision code with vectorised t_collide 
        subroutine. Still won't run properly for me (plot opens and immediately 
        closes). No error message or debugging.
    23/11/2016: WJ - Fixed new collision code, vectorised wall code. Separated
        functions to brownian_tools module.
    24/11/2016: WJ - Added in new collision routine for walls, incomplete.
    6/12/2016: WJ - Added new wall class with initialisation routines.
    6/12/2016: LM - Added sawtooth function with trigger wall on tooth side.
        Needs assimilating to rest of new code.
    
Variable naming convention:
    w_... Wall related variable (should probably make class)
    n_... # of ...
    p_... position of ...
    v_... velocity of ...
    ..._n normal vector of ...
    ..._p parallel vector of ...
"""

import sys
import numpy as np
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from brownian_tools import *

"""
Define # of dimensions
"""
nd = 2

"""
Define wall class
[move to different file?]
"""
class Wall:
    """
    Initialise with coordinate 0,0,...
    """
    def __init__(self, nd):
        #no. of dimensions
        self.nd = nd
        #corner co-ords
        self.co = np.full([1,self.nd],0)
        #Limits of box
        self.x_v = np.array([0,0])
        self.y_v = np.array([0,0])
        #no. of walls (# of co-ords -1)
        self.n = 0
        """TODO: add check for corresponding periodic boundary"""
        #Periodic boundary flag
        self.pb_flag = 0
        #Index of corresponding opposite boundary
        self.pb_ind = 0
        #Roughness value (randomly rotates wall vector for collision)
        self.rb = 0
        #Friction value (proportionally reduces parallel v)
        self.fb = 0

    """Get normalised wall vector for each wall"""
    def get_v(self):
        self.v = np.full([self.n,self.nd],np.nan)
        for i in np.arange(self.n):
            self.v[i] = self.co[i+1] - self.co[i]
        
        self.v = self.v / (np.sum(self.v**2, axis=1))**0.5

    """Get perpendicular normal vector for each wall (2d)"""
    def get_n(self):
        if nd == 2:
            self.n = np.concatenate((-self.v[:,1],self.v[:,0]), axis=1)
    
    """Set wall coords and find vectors for 2d case"""
    def set_2d(self):
        if np.sum(self.co[-1] != self.co[0]) >= 1:
            self.co = np.concatenate((self.co, self.co[-1]), axis=0)
        
        self.x_v = np.array([np.nanmin(self.co[:,0]),np.nanmax(self.co[:,0])])
        self.y_v = np.array([np.nanmin(self.co[:,1]),np.nanmax(self.co[:,1])])
        self.n = np.size(self.co)[0] - 1
        self.get_v()
        self.get_n()
        
    def rect(self,x0,y0,x1,y1):
        self.co = np.array([[x0,y0],
                            [x1,y0],
                            [x1,y1],
                            [x0,y1]])
        self.set_2d()
    
    
    
    #def check_inside(pos,vel):
        
 

"""
Set # of balls, radii and masses
Define as arrays for multiple different balls
"""
n_balls = np.array([100,5,1])
radii = np.array([0.025,0.1,0.2])
m_balls = np.array([5,100,400])

"""
Define number of steps (for plotting) and step length, find max time
"""
n_steps = 100
l_step = 0.01

t_max = l_step*n_steps

"""
Define wall coordinates
"""
wall = Wall(nd)

#Simple box
wall.co = np.array([[-1,-1],
                    [1,-1],
                    [1,1],
                    [-1,1]])

wall.set_2d()

"""
#Check if last corner equal to first, if not duplicate, then get number of 
#walls
if (sum(w_coord[0] == w_coord[-1]) != np.shape(w_coord)[1]):
    n_wall = np.shape(w_coord)[0]
    w_coord = np.concatenate((w_coord, w_coord[-1]), axis=0)
    
else:
    n_wall = np.shape(w_coord)[0] - 1

#Now calculate wall unit vectors normal vectors
w_vec = nanarr([n_wall,nd])
w_norm = nanarr([n_wall,nd])
for i in np.arange(n_wall):
    w_vec = w_coord[i+1]-w_coord[i]
    w_vec = w_vec/np.sum(w_vec**2)**0.5
    #Calculate normalised normal vector to wall
    w_norm = np.dot(np.array([[0,-1],[1,0]]), w_vec)
    w_norm = w_norm/np.sum(w_norm**2)**0.5
"""

# TODO: convert all routines to use new wall class

"""
Creates arrays of ball sizes and masses for input values
"""
tot_balls = np.sum(n_balls)
for i in range(len(n_balls)):
    if i == 0:
        size_arr = np.full(n_balls[0], radii[0])
        m_arr = np.full(n_balls[0], m_balls[0])
    else:
        size_arr = np.concatenate((size_arr, np.full(n_balls[i], radii[i])))
        m_arr = np.concatenate((m_arr, np.full(n_balls[i], m_balls[i])))

class Container:
    def __init__(self):
        self.x_v = 0
        self.y_v = 0
    
    def rect(self, vertices_array):
        """input two vertices corresponding to opposite corners of rectangle as
        list [[Xo,Yo],[X1,Y1]]"""
        
        self.x_v = np.asarray(vertices_array)[:,0]
        self.y_v = np.asarray(vertices_array)[:,1]
        self.x_v.sort()
        self.y_v.sort()
    
    def sawtooth(self, vertices_array, num_teeth, spike_side = 1, spike_dir = 1):
        """input two vertices to create wall between as list [[Xo,Yo],[X1,Y1]],
        no. of triangular teeth as integer,
        side want spikes on: 1 = right, -1 = left
        direction want spikes in: -1 reverses
        """
        #for horizontal lines, spike_side = 1 will put spikes on bottom, -1 on top
        #spikes will point up or right if spike_dir = 1, down or left for -1
        
        self.vertices = np.asarray(vertices_array)
        self.vdist = pdist(self.vertices)
        self.width = self.vdist / num_teeth
        self.height = self.width
        self.side = spike_side
        self.dir = spike_dir
        self.x_v2 = np.asarray(vertices_array)[:,0]
        self.y_v2 = np.asarray(vertices_array)[:,1]
        self.theta = 0
        
        self.norm = normal(self.x_v2, self.y_v2, self.side)
        self.norm = self.norm * length(self.norm)**(-1)
        self.trig_wall_x_v = self.x_v2 + self.norm[0] * (self.height + 1.1 * max(radii))
        self.trig_wall_y_v = self.y_v2 + self.norm[1] * (self.height + 1.1 * max(radii))
        #self.trig_wall = np.array
        
        #sets self.angle_sign = +- 1 depending on positive or negative slope
        self.angle_sign = np.sign((self.y_v2[1]-self.y_v2[0])/(self.x_v2[1]-self.x_v2[0]))
        if self.angle_sign == 0:
            self.angle_sign = 1
        
        #determines whether vertices are in ascending or descending order of y-value as multiplier
        c = 1
        if self.y_v2[1] < self.y_v2[0]:
            if self.x_v2[1] != self.x_v2[0]:
                c = -1
        
        #generates vertical sawtooth to be rotated and moved into place
        self.y_v = np.linspace(0, self.vdist, self.vdist * 1000)
        self.x_v = self.side*0.5*self.height*(signal.sawtooth(2 * np.pi * self.width**(-1) * self.y_v) + 1) 
        if self.dir == -1:
            self.x_v = self.x_v[::-1]
        
        #rotation if necessary
        if self.x_v2[0] != self.x_v2[1]:
            self.theta = findangle(np.array([self.x_v2[1]-self.x_v2[0], self.y_v2[1]-self.y_v2[0]]))
            self.x_prime = self.x_v * math.cos(self.theta) + c*self.angle_sign * self.y_v * math.sin(self.theta)
            self.y_prime = (-1) * c*self.angle_sign * self.x_v * math.sin(self.theta) + self.y_v * math.cos(self.theta)
        elif self.x_v2[0] == self.x_v2[1]:
            self.x_prime = self.x_v
            self.y_prime = self.y_v
        
        #shifts rotated sawtooth to correct vertex 
        if self.angle_sign == 1:
            self.x_fin = c * self.x_prime + min(self.x_v2)
            self.y_fin = c * self.y_prime + min(self.y_v2)
        else:
            self.x_fin = c * self.x_prime + max(self.x_v2)
            self.y_fin = c * self.y_prime + min(self.y_v2) 
        
#walls = Container()
#walls.rect([[-1,-1],[1,1]])
        
"""
Initialise position and velocity of balls
"""
p1 = np.random.uniform(low=wall.x_v[0]+max(radii), 
                       high=wall.x_v[1]-max(radii), size=[tot_balls,1])
p2 = np.random.uniform(low=wall.y_v[0]+max(radii), 
                       high=wall.y_v[1]-max(radii), size=[tot_balls,1])
p = np.concatenate((p1,p2), axis=1)
v = np.random.uniform(low=-0.5, high=0.5, size=[tot_balls,2])

#TODO check if balls are inside/outside of box

# Initialise time
t = 0

# Initialise indices for last collided balls, for plotting
ii = 0
jj = 0

# Find initial momentum and energy
momentum = [np.sum((np.sum(v**2, axis=1)**0.5))]
energy = np.sum((np.sum(v**2, axis=1)*m_arr)/2)

#start pressure counter and pressure at 0
pressure = 0
press_temp = 0
avg_press = 0
press1 = []

"""
Step function, step to next collision or step time limit
"""
def step():
    global p, t, energy, pressure, avg_press, press_temp, press1
    
    n = np.shape(v)[0]
    [j_arr,i_arr] = tri_indgen(n) #reversed i,j to use upper half triangle
    
    # Find collision times between balls
    t_col = t_collide(p,v,size_arr)
    
    # Find collision times for vertical walls
    t_2wall = nanarr([tot_balls,wall.n])
    for j in range(wall.n):
        w_ends = wall.co[j:j+1]
        for i in range(tot_balls):
            t_2wall[i,j] = t_wall(p[i],v[i],size_arr[i],w_ends)
    
    
    """
    t_wall_x = np.nanmax([(walls.x_v[0]-p[:,0]+size_arr)/v[:,0], 
                          (walls.x_v[1]-p[:,0]-size_arr)/v[:,0]], axis = 0)
    t_wall_x[t_wall_x < 0] = np.nan
    
    t_wall_y = np.nanmax([(walls.x_v[0]-p[:,1]+size_arr)/v[:,1], 
                          (walls.x_v[1]-p[:,1]-size_arr)/v[:,1]], axis = 0)
    t_wall_y[t_wall_y < 0] = np.nan
    """
    
    """Old loop code - slow
    for i in range(tot_balls):
        
        # Find collision times for vertical walls
        temp = np.nanmax([(walls.x_v[0]-p[i,0]+size_arr[i])/v[i,0], 
                          (walls.x_v[1]-p[i,0]-size_arr[i])/v[i,0]]) 
        #walls set as 0 & 1 here
        if temp >= 0:
            t_wall_x[i] = temp
        
        # Find collision times for horizontal walls
        temp = np.nanmax([-(p[i,1]-walls.y_v[0]-size_arr[i])/v[i,1],
                          -(p[i,1]-walls.y_v[1]+size_arr[i])/v[i,1]])
        if temp >= 0:
            t_wall_y[i] = temp
            
    
        
        
        Old iterative code, much slower. See t_collide function for replacement
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
        
    """
            
    # Find minimum times to collision from each
    t_wall_min = np.nanmin(t_2wall)
    t_col_min = np.nanmin(t_col)

    # Find overall minimum time to collision 
    t_min = np.nanmin([t_wall_min,t_col_min])

    # Check if no collision, if so stop
    if np.isnan(t_min):
        sys.exit("No collisions detected")

    # Plot every step up to collison time
    if l_step < t_min:
        t = t+l_step
        p += v*l_step
        
        #calculate pressure over last 100 l_steps if 100 have been taken
        press1.append(press_temp)
        if len(press1) < 100:
            pressure = np.sum(press1) / len(press1)
        else:
            pressure = np.sum(press1[len(press1):len(press1)-101:-1]) / 100.
        avg_press = np.sum(press1) / len(press1)
        press_temp = 0
        
        energy = np.sum((np.sum(v**2, axis=1)*m_arr)/2)
        
        return p, energy, pressure, avg_press
            
    else:
        #update time
        t += t_min

        #Move balls to collision positions
        p += v*t_min
        
        """Check if balls are overlapping each other (tolerance of 1E-10 for 
        float uncertainty)"""
        """#uncomment to enable
        temp = get_dist(p,size_arr)
        if np.size(temp[temp<-1E-10]) > 0:
            print temp[temp<0]
        """
        
        #Temporary wall code
        for j in range(wall.n):
            if np.sum(t_2wall == t_min) > 0:
                ii = np.where(t_2wall[:,j] == t_min)
                #This is only for box case
                v[ii,(j+1)%2] = -v[ii,(j+1)%2]
                #TODO make v_wall routine
                
                
        """
        if np.sum(t_wall_x == t_min) > 0:
            v[t_wall_x == t_x_min,0] = -v[t_wall_x == t_x_min,0]
            ii = np.where(t_wall_x == t_min)
            #find momentum add
            for x in ii:
                press_temp += sum(m_arr[x]*(sum(v[x]**2)**0.5)) * l_step / (walls.y_v[1]-walls.y_v[0])

        if np.sum(t_wall_y == t_min) > 0:
            v[t_wall_y == t_y_min,1] = -v[t_wall_y == t_y_min,1]
            jj = np.where(t_wall_y == t_min)
            for x in jj:
                press_temp += sum(m_arr[x]*(sum(v[x]**2)**0.5)) * l_step / (walls.x_v[1]-walls.x_v[0])
        """
        
        if np.sum(t_col == t_min) > 0:
            wh = t_col == t_min
            kk = i_arr[wh]
            ll = j_arr[wh]
            for m in range(kk.size):

                dij = p[ll[m]]-p[kk[m]]
                rij = dij/(np.sum(dij**2)**0.5)
                dm = m_arr[kk[m]]/m_arr[ll[m]]
            
                ua = np.sum((v[kk[m]]-v[ll[m]])*rij)
                
                va = ((dm-1)*ua)/(1+dm)
                dva = va - ua
            
                v[kk[m]] = v[kk[m]] + dva*rij
                v[ll[m]] = v[ll[m]] - dm*dva*rij
    
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
xlim=(min(wall.x_v)-1, max(wall.x_v)+1), ylim=(min(wall.x_v)-1, max(wall.y_v)+1))

# Enter the locations of the particles
particles, = ax.plot([], [], 'bo', ms=5.)
if len(radii) > 1:
    particles1, = ax.plot([], [], 'bo', ms=10.)
if len(radii) > 2:
    particles2, = ax.plot([], [], 'bo', ms=15.)
time_text = ax.text(0.025, 0.93, '', transform=ax.transAxes)
energy_text = ax.text(0.025, 0.88, '', transform=ax.transAxes)
pressure_text = ax.text(0.025, 0.83, '', transform=ax.transAxes)
avg_press_text = ax.text(0.025, 0.78, '', transform=ax.transAxes)

#Make box boundary (manually set for now, should automate later)
rect = plt.Rectangle([wall.x_v[0],wall.y_v[0]], wall.x_v[1] - wall.x_v[0], wall.y_v[1] - wall.y_v[0], lw=2, fc='none')

ax.add_patch(rect)

# Initialize animation
def init():
    global rect
    particles.set_data([], [])
    if len(radii) > 1:
        particles1.set_data([], [])
    if len(radii) > 2:
        particles2.set_data([], [])
    
    time_text.set_text('')
    energy_text.set_text('')
    pressure_text.set_text('')
    avg_press_text.set_text('')
    
    if len(radii) == 1:
        return particles, time_text, energy_text, pressure_text, avg_press_text, rect
    elif len(radii) == 2:
        return particles, particles1, time_text, energy_text, pressure_text, avg_press_text, rect
    else:
        return particles, particles1, particles2, time_text, energy_text, pressure_text, avg_press_text, rect

# Perform animation step
def animate(i):
    global ax, fig, rect
    
    #step function forward
    step()

    #set marker size based on ball size
    ms = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*radii[0]
    particles.set_markersize(ms)
    particles.set_data(p[0:n_balls[0], 0], p[0:n_balls[0], 1])
    if len(radii) > 1:
        ms1 = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*radii[1] #this abs could be the problem
        particles1.set_markersize(ms1)
        particles1.set_markerfacecolor('r')
        particles1.set_data(p[n_balls[0]:sum(n_balls[0:2]), 0], p[n_balls[0]:sum(n_balls[0:2]), 1])
    if len(radii) > 2:
        ms2 = fig.dpi * fig.get_figwidth()/(ax.get_xlim()[1] - ax.get_xlim()[0]) * 2*radii[2]
        particles2.set_markersize(ms2)
        particles2.set_data(p[sum(n_balls[0:2]):sum(n_balls[0:3]), 0], p[sum(n_balls[0:2]):sum(n_balls[0:3]), 1])
        
    time_text.set_text('Time = %.1f s' % t)
    energy_text.set_text('Energy = %.2f J' % energy)
    pressure_text.set_text('Pressure = %.2f mPa' % (pressure*1000.))
    avg_press_text.set_text('Average Pressure = %.1f mPa' % (avg_press*1000.))
    
    if len(radii) == 1:
        return particles, time_text, energy_text, pressure_text, avg_press_text, rect
    elif len(radii) == 2:
        return particles, particles1, time_text, energy_text, pressure_text, avg_press_text, rect
    else:
        return particles, particles1, particles2, time_text, energy_text, pressure_text, avg_press_text, rect

ani = animation.FuncAnimation(fig, animate, frames=600, interval=10, blit=True, init_func=init)
                              
plt.show()
