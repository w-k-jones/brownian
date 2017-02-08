# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 09:50:23 2017

@author: pc
"""
#

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
        t_2wall = self.wall.t2wall()
        
        if t_2col < t_2wall & t_2col < self.t_step:
            self.t +=t_2col
            self.ball.p += self.ball.v*t_2col
            self.ball.dv_col()
            
        else:
            if t_2col > t_2wall & t_2wall < self.t_step:
                self.st +=t_2wall
                self.ball.p += self.ball.v*t_2wall
                self.wall.dv_wall(self.ball)
            
            else:
                self.ball.p += self.ball.v*self.t_step
            
    