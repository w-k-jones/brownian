# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 09:50:23 2017

@author: pc
"""
#
import numpy as np
import matplotlib.pyplot as plt
#Define system class to hold wall,ball classes and system procedures

class system:
    def __init__(self,wall_in,ball_in):
        self.wall = wall_in
        self.ball = ball_in
        self.t = 0
        self.t_step = 0.1
        self.t_max = 100
        
    def step(self):
        t_2col = self.ball.t2col()[0]
        t_2wall = self.wall.t2wall(self.ball)[0]
        
        if (t_2col < t_2wall) & (t_2col < self.t_step):
            print 'ball collision'
            self.t +=t_2col
            self.ball.p += self.ball.v*t_2col
            self.ball.dv_col()
            
        else:
            if (t_2col > t_2wall) & (t_2wall < self.t_step):
                'wall collision'
                self.t +=t_2wall
                self.ball.p += self.ball.v*t_2wall
                self.wall.dv_wall(self.ball)
            
            else:
                'no collision'
                self.t +=t_step
                self.ball.p += self.ball.v*self.t_step
            
    def plot(self):
        plt.scatter(self.ball.p[:,0],self.ball.p[:,1])
        plt.plot(self.wall.co_plt[:,0],self.wall.co_plt[:,1])
        plt.show
    
    def run(self,n_step):
        mo = np.full([n_step+1,2],np.nan)
        #ang = np.full(n_step+1,np.nan)
        
        mo[0] = np.sum(self.ball.v * np.reshape(self.ball.r,[self.ball.n,1]),axis=0)
        
        for i in range(n_step):
            self.step()
            mo[i+1] = np.sum(self.ball.v * np.reshape(self.ball.r,[self.ball.n,1]),axis=0)
        
        return mo
            
    
        
        