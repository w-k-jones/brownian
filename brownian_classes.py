# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 10:26:32 2016
@author: wj
Module containing classes 'wall_shape' and 'balls' for brownian system. 
Wall_shape is initialised using a 2d array of co-ordinates of the wall 
vertices. Balls is initialised using arrays/scalars of the ball #s, radii and 
masses and the wall object of the system it is contained within.
History:
    10/02/2017, wj: Commenting, adding corner code.
    16/02/2017, lm; Small adjustments so can run properly with animation file
"""

import numpy as np

class wall_shape:
    #Wall_shape class, contains info on wall properties and the functions for 
    #calculating times to collisions with walls and the resultant velocity 
    #change.
    
    """Get normalised wall vector for each wall"""
    def get_vec(self):
        #Caculates unit vector of each wall and the length.
        self.vec = self.co_plt[1:] - self.co_plt[0:-1]
        self.vlen = np.sum(self.vec**2,axis=1)**0.5
        
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
        
    """
    Calculate area contained within wall
    """
    def get_A(self):
        tmp = (-self.vec[:,0]*self.co[:,1]*self.vlen
               +self.vec[:,1]*self.co[:,0]*self.vlen)
    
        self.A = 0.5*np.sum(tmp)
        
        #print self.A
    
    #Set wall coords and find vectors for 2d case
    def set_2d(self):
        #no. of dimensions
        self.nd = np.shape(self.co)[1]

        print 'Dimensionality = ',self.nd
        
        #Check if first and last coordinates are the same. If they are remove 
        #last co-ordinate for list of vertices (.co)
        if np.sum(self.co[-1] != self.co[0]) == 0:
            self.co = self.co[0:-1]
        
        #list of vertices for plotting, needs first co-ordinate appended on the end
        self.co_plt = np.concatenate((self.co, np.reshape(self.co[0],[1,self.nd])), axis=0)
        
        #note, done this way to avoid a long if expression
        
        #no. of walls = #coords in closed shape
        self.n = np.shape(self.co)[0]
        print 'Vertices = ',self.n

        #TODO: add check for corresponding periodic boundary - subroutine needed
        #Periodic boundary flag
        self.pb_flag = np.full(self.n,0)
        #Index of corresponding opposite boundary
        self.pb_ind = np.full(self.n,0)
        #Roughness value (randomly rotates wall vector for collision)
        #angle of 1 s.d. note, angle cannot be more than pi/4
        self.rb = np.full(self.n,0.)
        #Friction value (proportionally reduces parallel v)
        self.fb = np.full(self.n,0.)
        #wall temperture
        self.T = np.full(self.n,np.nan)
        #wall pressure
        self.P = np.full(self.n,0.)
        
        
        #x and y limits - smallest and largest co-ords in each direction
        self.xlim = np.array([np.nanmin(self.co[:,0]),np.nanmax(self.co[:,0])])
        self.ylim = np.array([np.nanmin(self.co[:,1]),np.nanmax(self.co[:,1])])
        print 'setting system limits'
        
        #get wall vectors and normals.
        print 'finding wall vectors'
        self.vec = self.get_vec()
        print 'finding wall normal vectors'
        self.norm = self.get_norm()
        self.get_A()
        print 'System area = ',self.A
        
   
        
    """
    Initialise with coordinate array
    """
    def __init__(self, coordinate_array):
        #corner co-ords
        self.co = coordinate_array
        
        print '>>>Initialising wall properties'
        
        #now set all values
        self.set_2d()
        
        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    
        
        

    
    """calculates time to collision with corner in case of reflex angles"""
    def t2corn(self, ball):
        #check for different dimensionality (e.g. cylinder)
        if ball.nd > self.nd:
            p_b = ball.p[0:self.nd]
            p_b[-1] = (np.sum((ball.p[self.nd-1:])**2,axis=1))**0.5
            
            v_b = ball.v[0:self.nd]
            v_b[-1] = (np.sum((ball.v[self.nd-1:])**2,axis=1))**0.5
            
        else:
            p_b = ball.p
            v_b = ball.v
        
        #Define temporry array for times
        temp = np.full([ball.n,self.n,2],np.nan)

        #Calculate quadratic a,b,c coefficients (similar to 2 ball collision)
        a = np.reshape(np.sum((v_b)**2, axis=1),[ball.n,1])
        b = 2*np.sum(np.reshape(v_b,[ball.n,1,self.nd])* 
                     (np.reshape(p_b,[ball.n,1,self.nd]) 
                      - np.reshape(self.co,[1,self.n,self.nd])), axis=2)
        c = (np.sum((np.reshape(p_b,[ball.n,1,self.nd])
                    - np.reshape(self.co,[1,self.n,self.nd]))**2, axis=2)
                    - (np.reshape(ball.r**2,[ball.n,1])))
    
        #Calculate discriminant
        chk = b**2 - 4*a*c
        wh = chk >= 0
        ti,tj = np.where(wh)
        
        #Find both solutions and take smallest (first collision)
        temp[wh,0] = (-b[wh]+chk[wh]**0.5)/(2*np.tile(a,[1,self.n])[wh])
        temp[wh,1] = (-b[wh]-chk[wh]**0.5)/(2*np.tile(a,[1,self.n])[wh])
        temp = np.nanmin(temp, axis=2)
        #If negative solution is in the past, so set as NaN to ignore.
        temp[temp<=0] = np.nan
        
        #Find smallest positive solution
        t_min = np.nanmin(temp)
        
        #Set values for ball and corner indexes to use in changing velocities.
        self.i_ball_c = np.where(temp == t_min)[0]
        self.i_corn = np.where(temp == t_min)[1]

        #Return time
        return t_min#,self.i_arr[wh],self.j_arr[wh]
        
    """Change velocity of ball after corner collision"""
    def dv_corn(self,ball):
        #Calculate normal vector from ball to vertex
        t_norm = self.co[self.i_corn] - ball.p[self.i_ball_c]
        #Normalise
        t_norm = t_norm/(np.sum(t_norm**2)**0.5)
        
        #adjust for roughness value
        if self.rb[self.i_corn] > 0:
            ang = np.pi
            while ang > np.pi/4.:
                ang = np.random.normal(loc=0.,scale=self.rb[self.i_corn])
                
            t_norm[0] = t_norm[0]*np.cos(ang) -t_norm[1]*np.sin(ang)
            t_norm[1] = t_norm[0]*np.sin(ang) +t_norm[1]*np.cos(ang)
        
        #Dot product of ball velocity with normal
        dv = np.sum(ball.v[self.i_ball_c]*t_norm)
        
        #adjust for wall temperature
        if np.isfinite(self.T[self.i_corn]):
            sig = (self.T[self.i_corn]/ball.m[self.i_ball_c])**2
            dv_T = abs(np.random.normal(loc=0.,scale=sig))
            #while dv_T/dv >= 0.5:
                #dv_T = np.random.normal(sig)
            
            if dv <=0:
                dv -=dv_T
            else:
                dv +=dv_T
        else: dv = 2*dv
        
        #print ball.v[self.i_ball_c]
        #print t_norm
        #print dv *t_norm
        #change ball velocity by 2x dv in direction of -normal
        ball.v[self.i_ball_c] -= dv*t_norm
        #print ball.v[self.i_ball_c]
        
        return dv
        
    
    """Calculate time to next collision between ball and wall"""
    #new, vectorised and object code
    def t2wall(self, ball):
        #check for different dimensionality (e.g. cylinder)
        if ball.nd > self.nd:
            p_b = ball.p[0:self.nd]
            p_b[-1] = (np.sum((ball.p[self.nd-1:])**2,axis=1))**0.5
            
            v_b = ball.v[0:self.nd]
            v_b[-1] = (np.sum((ball.v[self.nd-1:])**2,axis=1))**0.5
            
        else:
            p_b = ball.p
            v_b = ball.v
            
        #Calculate relative normal positions and velocities
        p_rel = np.sum((np.reshape(self.co,[1,self.n,self.nd]) 
                        - np.reshape(p_b,[ball.n,1,self.nd]))
                        * np.reshape(self.norm,[1,self.n,self.nd]),axis=2)
                    
        v_rel = np.sum((np.reshape(v_b,[ball.n,1,self.nd]))
                        * np.reshape(self.norm,[1,self.n,self.nd]),axis=2)
        
        #Time to collision with infinite length wall
        #Subtract time taken for ball to travel one radius to account for size 
        #of ball
        t = p_rel/v_rel - np.abs(np.reshape(ball.r,[ball.n,1])/v_rel)
        
        #Now calculate position along length of wall at collision
        
        #Parallel position and velocity components
        p_par = np.sum((-np.reshape(self.co,[1,self.n,self.nd]) 
                        + np.reshape(p_b,[ball.n,1,self.nd]))
                        * np.reshape(self.vec,[1,self.n,self.nd]),axis=2)
        
        v_par = np.sum((np.reshape(v_b,[ball.n,1,self.nd]))
                        * np.reshape(self.vec,[1,self.n,self.nd]),axis=2)
        
        #Position at collision divided by wall length
        p_col = (p_par + v_par*t)/np.reshape(self.vlen,[1,self.n])
        
        #If >1 or <0 then collision is not along the length of the wall, so
        #ignore
        t[p_col < 0] = np.nan
        t[p_col > 1] = np.nan
                                 
        #Negative time is in the past so ignore
        t[t<=0] = np.nan
                                 
        #Calculate time to next collision
        t_min = np.nanmin(t)
        
        #Find indices of ball and wall section and save for dv calculation
        [i_ball,i_wall] = np.where(t == t_min)
        self.i_ball = i_ball
        self.i_wall = i_wall
        
        #return next time
        return t_min#,i_ball,i_wall
        
    """Calculate if a ball is inside the system"""
    #TODO: debug ball appearing half way through line on way into system.
    def isinside(self,p_in,v_in,r_in):
        #repeats t2wall but for set p,v
        if p_in.size > self.nd:
            #print self.nd
            p_b = np.copy(p_in[0:self.nd])
            #print p_b
            p_b[-1] = (np.sum((p_in[(self.nd-1):])**2))**0.5
            
            v_b = np.copy(v_in[0:self.nd])
            v_b[-1] = (np.sum((v_in[(self.nd-1):])**2))**0.5
            #print p_b
            
        else:
            p_b = p_in
            v_b = v_in
            
        print p_in
        print p_b
        p_rel = np.sum((self.co - np.reshape(p_b,[1,self.nd]))
                       * np.reshape(self.norm,[self.n,self.nd]),axis=1)
        
        v_rel = np.sum(np.reshape(v_b,[1,self.nd])
                       * np.reshape(self.norm,[self.n,self.nd]),axis=1)
                       
        #print p_rel,v_rel
        
        t = p_rel/v_rel
        #print t
        
        r_a = np.absolute(r_in/v_rel)
        
        t -= r_a
        
        p_par = np.sum((-self.co + np.reshape(p_b,[1,self.nd]))
                        * np.reshape(self.vec,[self.n,self.nd]),axis=1)
        
        v_par = np.sum(np.reshape(v_b,[1,self.nd])
                       * np.reshape(self.vec,[self.n,self.nd]),axis=1)
        
        p_col = (p_par + v_par*t)/self.vlen

        t[p_col < 0] = np.nan
        t[p_col > 1] = np.nan
                  
        t[t<0] = np.nan
        
        print t

        chk = np.sum(np.isfinite(t)) % 2
        #1 if inside, 0 if outside

        if chk == 1:
            #double check backwards
            t2 = p_rel/(-v_rel) - r_a
            
            p_col2 = (p_par -v_par*t2)/self.vlen
            t2[p_col2 < 0] = np.nan
            t2[p_col2 > 1] = np.nan
            t2[t2<0] = np.nan
            print t2
            chk = np.sum(np.isfinite(t2)) % 2    
        """
            #check for corner collision
            if chk == 1:
                #Calculate quadratic a,b,c coefficients (similar to 2 ball collision)
                a = np.sum(v_b**2)
                b = 2*np.sum((v_b*p_b)- self.co, axis=1)
                c = np.sum((p_b-self.co)**2, axis=1) - r_in
    
                #Calculate discriminant
                chk_c = b**2 - 4*a*c
                #print chk_c
                if np.sum(chk_c >= 0) > 0:
                    chk = 0
        """
        return chk
        
    """Change ball velocity after collision with wall"""
    def dv_wall(self, ball):
        #Account for different dimensions
        if ball.nd > self.nd:
            v_b = ball.v[self.i_ball,0:self.nd]
            v_b[-1] = (np.sum((ball.v[self.i_ball,self.nd-1:])**2,axis=1))**0.5
            
        else:
            v_b = ball.v[self.i_ball]

        #set temporary wall normal
        t_norm = self.norm[self.i_wall]

        #adjust for roughness value
        if self.rb[self.i_wall] > 0:
            ang = np.pi
            while ang > np.pi/4.:
                ang = np.random.normal(loc=0.,scale=self.rb[self.i_wall])
                
            t_norm[0] = t_norm[0]*np.cos(ang) -t_norm[1]*np.sin(ang)
            t_norm[1] = t_norm[0]*np.sin(ang) +t_norm[1]*np.cos(ang)

        #Calculate normal dv
        
        chk = 0
        #adjust for wall temperature
        if np.isfinite(self.T[self.i_wall]):
            dv = np.sum(v_b*t_norm)
            
            while chk == 0:
                if self.T[self.i_wall] > 0:
                    sig = (self.T[self.i_wall]/ball.m[self.i_ball])**2
                    u = dv-np.random.normal(loc=0.,scale=sig)
                else: u = dv.copy
            
                w = ((ball.m[self.i_ball]-1)*u)/(ball.m[self.i_ball]+1)
                dv_tot = w-u
                
                ball.v[self.i_ball] = ball.v[self.i_ball] + dv_tot*t_norm
                v_b = ball.v[self.i_ball]
                self.P[self.i_wall] += ball.m[self.i_ball]*abs(dv_tot)
                if dv >= 0:
                    if np.sum(v_b*t_norm) < 0:
                        chk = 1
                else:
                    if np.sum(v_b*t_norm) > 0:
                        chk = 1
                    
        else: 
            dv = 2*np.sum(v_b*t_norm)
        
            ball.v[self.i_ball] = ball.v[self.i_ball] - dv*t_norm

            self.P[self.i_wall] += ball.m[self.i_ball]*abs(dv)

        return dv
        
        """ball collision code for wall with temperature
        u = np.sum((self.v[self.i_ball]-self.v[self.j_ball])*rij)
        
        w = ((dm-1)*u)/(1+dm)
        dv = w - u
        
        self.v[self.i_ball] = self.v[self.i_ball] + dv*rij
        """
    
class balls:
    #generate balls
    def __init__(self,n_in,r_in,m_in,nd,wall):
        print '>>>Initialising particle properties'
        self.nd = nd
        print 'Dimensionality = ',self.nd
        self.n = np.sum(n_in)
        print 'Number of particles = ',self.n
        
        print 'setting particle sizes and masses'
        if np.size(n_in) == 1:
            self.r = np.full(n_in, r_in)
            self.m = np.full(n_in, m_in)
            self.r2 = np.array([r_in])
            self.n_balls = np.array([n_in])
        
        else:
            self.r2 = r_in
            self.n_balls = n_in
            for i in range(np.size(n_in)):
                if i == 0:
                    self.r = np.full(n_in[0], r_in[0])
                    self.m = np.full(n_in[0], m_in[0])
                else:
                    self.r = np.concatenate((self.r, np.full(n_in[i], r_in[i])))
                    self.m = np.concatenate((self.m, np.full(n_in[i], m_in[i])))
        
        print 'Initialising position and velocity'
        self.p = np.full([self.n,self.nd],np.nan)
        self.v = np.full([self.n,self.nd],np.nan)  
        
        j = 0
        cnt = 0
        while j < self.n:
            
            cnt +=1
            
            if self.nd == 2:
                new_p = np.array([np.random.uniform(low=wall.xlim[0]+self.r[j], 
                               high=wall.xlim[1]-self.r[j]),
                              np.random.uniform(low=wall.ylim[0]+self.r[j], 
                               high=wall.ylim[1]-self.r[j])])
                new_v = np.random.normal(loc=0., scale=1.6, size=[1,2])
                #new_p = np.reshape(new_p,[1,2])
            elif self.nd == 3:
                new_p = np.array([np.random.uniform(low=wall.xlim[0]+self.r[j], 
                               high=wall.xlim[1]-self.r[j]),
                              np.random.uniform(low=wall.ylim[0]+self.r[j], 
                               high=wall.ylim[1]-self.r[j]),
                              np.random.uniform(low=wall.ylim[0]+self.r[j], 
                               high=wall.ylim[1]-self.r[j])])
#<<<<<<< HEAD
                new_p = np.reshape(new_p,[1,3])
                new_v = np.random.normal(loc=0., scale=1.6, size=[1,3])
#=======
                #new_p = np.reshape(new_p,[1,3])
            new_v = np.random.normal(loc=0., scale=1, size=self.nd)
#>>>>>>> origin/new_classes
            
            if wall.isinside(new_p,new_v,self.r[j]) == 1:
                """if (np.nanmin(np.sum((self.p - np.reshape(new_p,[1,self.nd]))**2) 
                           - (self.r + self.r[j]))) < 0:
                    
                    self.p[j] = new_p
                    self.v[j] = new_v
                    j+=1
                    print j
                
                else: print "ball overlap" """
                self.p[j] = new_p
                self.v[j] = new_v
                j+=1
                print j
            #else: print "outside"
#<<<<<<< HEAD
            
        print 'Particles initialised: ',j
        print 'Attempts: ',cnt
        #print self.p
#=======
#>>>>>>> origin/new_classes
                    
        a = np.arange(self.n)
        a = np.tile(a,[self.n,1])
        b = a.T
        self.i_arr = a[a>b]
        self.j_arr = b[a>b]

        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                
            
    def t2col(self):
        temp = np.full([self.i_arr.size,2],np.nan)

        a = np.sum((self.v[self.i_arr] - self.v[self.j_arr])**2, axis=1)
        b = 2*np.sum((self.v[self.i_arr] - self.v[self.j_arr])
                     *(self.p[self.i_arr] - self.p[self.j_arr]), axis=1)
        c = (np.sum((self.p[self.i_arr]-self.p[self.j_arr])**2, axis=1)
             - (self.r[self.i_arr]+self.r[self.j_arr])**2)
    
        chk = b**2 - 4*a*c
        wh = chk >= 0
    
        temp[wh,0] = (-b[wh]+chk[wh]**0.5)/(2*a[wh])
        temp[wh,1] = (-b[wh]-chk[wh]**0.5)/(2*a[wh])
        temp = np.nanmin(temp, axis=1)
        temp[temp < 1E-10] = np.nan
        
        t_min = np.nanmin(temp)
        
        wh = np.where(temp == t_min)[0]

        self.i_ball = self.i_arr[wh]
        self.j_ball = self.j_arr[wh]

        return t_min#,self.i_arr[wh],self.j_arr[wh]

    def dv_col(self):
        dij = self.p[self.j_ball] - self.p[self.i_ball]
        rij = dij/(np.sum(dij**2)**0.5)
        dm = self.m[self.i_ball]/self.m[self.j_ball]

        u = np.sum((self.v[self.i_ball]-self.v[self.j_ball])*rij)
        
        w = ((dm-1)*u)/(1+dm)
        dv = w - u
        
        self.v[self.i_ball] = self.v[self.i_ball] + dv*rij
        self.v[self.j_ball] = self.v[self.j_ball] - dm*dv*rij
