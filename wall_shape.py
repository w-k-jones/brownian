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
        self.vec = self.vec / np.reshape(self.vlen,[self.n,1])
        
        return self.vec

    """Get perpendicular normal vector for each wall (2d)"""
    def get_norm(self, side=1):
        """side gives wall side of normal"""
        #Maths for collisions is based on normal pointing in the same direction
        #to the wall each time (90degrees CCW). Should be invariant though
        
        self.norm = np.concatenate((-np.reshape(self.vec[:,1],[self.n,1]),np.reshape(self.vec[:,0],[self.n,1])), axis=1)
        #Same as multiplying all by rotation matrix [0,-1]
        #                                           [1, 0]
        
        return self.norm
    
    #Set wall coords and find vectors for 2d case
    def set_2d(self):
        #no. of dimensions
        self.nd = np.shape(self.co)[1]
        
        #Check if first and last coordinates are the same (should be), if not 
        #stick first coordinate on the bottom of the coordinate array
        if np.sum(self.co[-1] != self.co[0]) >= 1:
            self.co = np.concatenate((self.co, np.reshape(self.co[0],[1,self.nd])), axis=0)
        
        #no. of walls = #coords in closed shape
        self.n = np.shape(self.co)[0]-1

        #TODO: add check for corresponding periodic boundary - subroutine needed
        #Periodic boundary flag
        self.pb_flag = np.full(0,self.n)
        #Index of corresponding opposite boundary
        self.pb_ind = np.full(0,self.n)
        #Roughness value (randomly rotates wall vector for collision)
        self.rb = np.full(0.,self.n)
        #Friction value (proportionally reduces parallel v)
        self.fb = np.full(0.,self.n)
        
        if np.sum(self.co[-1] != self.co[0]) >= 1:
            self.co = np.concatenate((self.co, self.co[-1]), axis=0)
        
        self.x_lim = np.array([np.nanmin(self.co[:,0]),np.nanmax(self.co[:,0])])
        self.y_lim = np.array([np.nanmin(self.co[:,1]),np.nanmax(self.co[:,1])])
        self.n = np.shape(self.co)[0] - 1
        self.vec = self.get_vec()
        self.norm = self.get_norm()
        
    """
    Initialise with coordinate array
    """
    def __init__(self, coordinate_array):
        #corner co-ords
        self.co = coordinate_array
        
        #now set all values
        self.set_2d()


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
    
