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
        t_2col = self.ball.t2col()
        t_2wall = self.wall.t2wall(self.ball)
        t_2corn = self.wall.t2corn(self.ball)
        
        self.min_t = np.nanmin([t_2col,t_2wall,t_2corn])
        
        if self.min_t <= self.t_step:
            self.t +=self.min_t
            self.ball.p += self.ball.v*self.min_t
            if t_2col == self.min_t:
                self.ball.dv_col()
                print 'b'
                
            if t_2wall == self.min_t:
                self.wall.dv_wall(self.ball)
                print 'w'
                
            if t_2corn == self.min_t:
                self.wall.dv_corn(self.ball)
                print 'c'
            
        else:
            self.t +=self.t_step
            self.ball.p += self.ball.v*self.t_step
            
    def plt(self):
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                             xlim=(self.wall.xlim[0]-0.1,self.wall.xlim[1]+0.1),
                             ylim=(self.wall.ylim[0]-0.1,self.wall.ylim[1]+0.1))
        wal, = ax.plot(self.wall.co_plt[:,0],self.wall.co_plt[:,1])
        bal, = ax.plot(self.ball.p[:,0],self.ball.p[:,1],'bo',
                       ms=fig.dpi
                       *fig.get_figwidth()/(ax.get_xlim()[1]-ax.get_xlim()[0])
                       *2*self.ball.r[0])
        fig.show
    
    def run(self,n_step):
        mo = np.full([n_step+1,2],np.nan)
        #ang = np.full(n_step+1,np.nan)
        
        mo[0] = np.sum(self.ball.v * np.reshape(self.ball.r,[self.ball.n,1]),axis=0)
        
        for i in range(n_step):
            self.step()
            mo[i+1] = np.sum(self.ball.v * np.reshape(self.ball.r,[self.ball.n,1]),axis=0)
        
        return mo
        
    def run_plt(self,n_step):
        for i in range(n_step):
            self.step()
            if i % 50 == 0:
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
                                     xlim=(self.wall.xlim[0]-0.1,self.wall.xlim[1]+0.1),
                                     ylim=(self.wall.ylim[0]-0.1,self.wall.ylim[1]+0.1))
                wal, = ax.plot(self.wall.co_plt[:,0],self.wall.co_plt[:,1])
                bal, = ax.plot(self.ball.p[:,0],self.ball.p[:,1],'bo',
                           ms=fig.dpi
                              *fig.get_figwidth()/(ax.get_xlim()[1]-ax.get_xlim()[0])
                              *2*self.ball.r[0])
                fig.show
            
        
    """    
    def ani_init(self):
        def init():
            outline = wall.co
            
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

            
    
        
        in_co3=np.array([[0,0],[0.2,0],[0.4,0.2],[0.4,0],[0.6,0.2],[0.6,0],[0.8,0.2],[0.8,0],[1,0],[1,1],[0.8,1],[0.6,0.8],[0.6,1],[0.4,0.8],[0.4,1],[0.2,0.8],[0.2,1],[0,1]])
        """