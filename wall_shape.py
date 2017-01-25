# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 10:26:32 2016

@author: pc
"""

import numpy as np
import math as math

class wall_shape:
    """Get normalised wall vector for each wall"""
    def get_vec(self):
        #self.vec = np.full([self.n,self.nd],np.nan)
        self.vec = self.co[1:] - self.co[0:-1]
        self.vlen = np.sum(self.vec,axis=1)
        
        #vectorised form fix
        self.vec = self.vec / np.tile(np.reshape(self.vlen,[1,self.n]),[1,self.nd])
        
        return self.vec

    """Get perpendicular normal vector for each wall (2d)"""
    def get_norm(self, side=1):
        """side gives wall side of normal"""
        #Maths for collisions is based on normal pointing in the same direction
        #to the wall each time (90degrees CCW). Should be invariant though
        
        self.norm = np.concatenate((-self.vec[:,1],self.vec[:,0]), axis=1)
        #Same as multiplying all by rotation matrix [0,-1]
        #                                           [1, 0]
        
        return self.norm
    
    #Set wall coords and find vectors for 2d case
    def set_2d(self):
        if np.sum(self.co[-1] != self.co[0]) >= 1:
            self.co = np.concatenate((self.co, self.co[-1]), axis=0)
        
        self.x_lim = np.array([np.nanmin(self.co[:,0]),np.nanmax(self.co[:,0])])
        self.y_lim = np.array([np.nanmin(self.co[:,1]),np.nanmax(self.co[:,1])])
        self.n = np.size(self.co)[0] - 1
        self.vec = self.get_vec()
        self.norm = self.get_norm()
        
    """
    Initialise with coordinate array
    """
    def __init__(self, coordinate_array):
        #corner co-ords
        self.co = coordinate_array
        
        #Check if first and last coordinates are the same (should be), if not 
        #stick first coordinate on the bottom of the coordinate array
        if np.sum(self.co[-1] != self.co[0]) >= 1:
            self.co = np.concatenate((self.co, self.co[0]), axis=0)
            
        #no. of dimensions
        self.nd = np.size(self.co)[1]
        
        #no. of walls = #coords in closed shape
        #But if giving finalcoord same as first this is wrong - need to correct
        #for this
        self.n = np.size(self.co)[0]-1

        #TODO: add check for corresponding periodic boundary - subroutine needed
        #Periodic boundary flag
        self.pb_flag = np.full(0,self.n)
        #Index of corresponding opposite boundary
        self.pb_ind = np.full(0,self.n)
        #Roughness value (randomly rotates wall vector for collision)
        self.rb = np.full(0,self.n)
        #Friction value (proportionally reduces parallel v)
        self.fb = np.full(0,self.n)
        #x and y limits of system (smallest and largest x and y values)
        self.x_lim = np.full(0,[1,self.nd])
        self.y_lim = np.full(0,[1,self.nd])


    def t_2wall(self, pos, vel, size):
        #Input dimensions (number of walls and number of particles)
        di = np.size(self.co)[0]-1
        dj = np.size(pos)[0]
        
        t = np.full([dj,di], np.nan)
        
        for i in np.arange(di):
            for j in np.arange(dj):
                
                #position from first corner of wall
                dp0 = pos[i]-self.co[j]
                #normal distance to wall
                p_norm = np.sum(dp0*self.norm[j])
                #normal velocity to wall
                v_norm = np.sum(vel[i]*self.norm[j])
                #time to wall
                dt = np.nanmin([-(p_norm-size[i])/v_norm, 
                                -(p_norm+size[i])/v_norm])
                
                if dt > 1E-10:
                    #Parallel postion to wall
                    p_par = np.sum(dp0*self.vec[j])
                    #Velocity parallel to wall
                    v_par = np.sum(vel[i]*self.vec[j])
                    #position along wall when colliding
                    p_wall_col = p_par +dt*v_par
                    
                    if p_wall_col >= 0 or p_wall_col <= self.vlen[j]:
                        t[j,i] = dt
                    
        return t
        
    def dv_wall(self, vel, jj):
        vel_n = vel/np.sum(vel**2)**0.5
        
        dv = np.sum(vel*self.norm[jj])
        
        new_vel = vel -2*dv*vel_n
        
        return new_vel
        
                    
        
    def check_inside(self, pos,vel,size):
        #Find times t wall collisions
        t = self.t_2wall(pos,vel,size)
        
        #Find no. of times negative (equals nan) for each particle.
        chk = sum(t == np.nan, axis=1)
        #Find where even/odd (odd = 1, even =0). Equivalent to logical array
        chk = chk%2
        
        return chk
    
 
class wall_modifier(wall_shape):    
    def flat(self, starting_vertex):
        """input wall index"""
        
        self.x_v = self.co[starting_vertex: starting_vertex+2, 0]
        self.y_v = self.co[starting_vertex: starting_vertex+2, 1]
        #self.x_v.sort()   #do we need these?
        #self.y_v.sort()
    
    def sawtooth(self, starting_vertex, num_teeth, spike_side = 1, spike_dir = 1):
        """input wall index,
        no. of triangular teeth as integer,
        side want spikes on: 1 = right, -1 = left
        direction want spikes in: -1 reverses
        """
        #for horizontal lines, spike_side = 1 will put spikes on bottom, -1 on top
        #spikes will point up or right if spike_dir = 1, down or left for -1
        #only works in 2d rn 
        
        self.vertices = self.co[starting_vertex: starting_vertex+2]
        self.vdist = pdist(self.vertices)
        self.width = self.vdist / num_teeth
        self.height = self.width
        self.side = spike_side
        self.dir = spike_dir
        self.x_v2 = self.vertices[:,0]
        self.y_v2 = self.vertices[:,1]
        self.theta = 0
        
        self.normal = self.get_n()[starting_vertex]
        self.trig_wall_x_v = self.x_v2 + self.normal[0] * (self.height + 1.1 * max(radii))
        self.trig_wall_y_v = self.y_v2 + self.normal[1] * (self.height + 1.1 * max(radii))
        
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
            
        self.x_v = self.trig_wall_x_v #use these and compare to x_fin after trigger?
        self.y_v = self.trig_wall_y_v
        
        #problem: trigger wall doesn't cover all balls if walls not at right angles