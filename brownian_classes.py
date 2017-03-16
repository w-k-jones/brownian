# -*- coding: utf-8 -*-
"""
Created on Wed Dec 07 10:26:32 2016
@author: wj
Module containing classes 'wall_shape' and 'balls' for brownian system. 
Wall_shape is initialised using a 2d array of co-ordinates of the wall 
vertices. Balls is initialised using arrays/scalars of the ball #s, radii and 
masses and the wall object of the system it is contained within.
History:
    10/02/2017, WJ: Commenting, adding corner code.
    16/02/2017, LM; Small adjustments so can run properly with animation file
    06/03/2017, WJ: Tidying, commenting and cleaning up code
    07/03/2017, WJ: Added periodic boundaries and more tidying
"""

import numpy as np

"""
Wall_shape class, contains wall property variables and the functions for 
calculating times to collisions with walls and the resultant velocity 
change.
"""
class wall_shape:
    
    """
    Get normalised wall vector for each wall
        Finds unit vector along each wall and the length of each wall
    """
    def get_vec(self, co):
        #create array of co-ordinates with repeat of first co at end
        co_2 = np.concatenate((co, np.reshape(co[0],[1,self.nd])), axis=0)
        #subtract end of each wall from the start
        vec = co_2[1:] - co_2[0:-1]
        #calculate the magnitde of each vector = wall length
        vlen = np.sum(vec**2,axis=1)**0.5
        #normalise to unit vector
        vec = vec / np.reshape(vlen,[np.shape(co)[0],1])
        return vec, vlen

    """
    Get perpendicular normal vector for each wall (2d):
        rotates each wall unit vector by 90 degrees CCW to find normal vectors.
    """
    def get_norm(self, vec):
        n = np.shape(vec)[0]
        norm = np.concatenate((-np.reshape(vec[:,1],[n,1]),
                               np.reshape(vec[:,0],[n,1])), axis=1)
        #Same as multiplying all by rotation matrix [0,-1]
        #                                           [1, 0]
        return norm
        
    """
    Calculate area contained within wall:
        Uses wall co-ordinates, vectors and lengths to perform line integral of
        field with constant curl around the boundary, which - by Stokes theorem
        - is equal to the area times a constant
        In this case the field is [-y,x] and the constant is 0.5
    """
    def get_A(self, co):
        #get wall unti vectors and lengths
        vec, vlen = self.get_vec(co)
        #cross product vectors with [-y,x]
        tmp = (-vec[:,0]*co[:,1]*vlen
               +vec[:,1]*co[:,0]*vlen)
        #sum (discrete integral) over all walls and multiply by 0.5 giving area
        A = 0.5*np.sum(tmp)
        A = np.abs(A)
        return A
    
    """
    Set initial wall properties:
        Uses input co-ordinates to set wall properties of:
            nd, n, coordinates, plotting coordinates, unit vectors, wall 
            lengths, normal vectors, x and y limits and area
        Also initialises blank arrays for a number of properties
    """
    def __init__(self, coordinate_array):
        print '>>>Initialising wall properties'
        self.co = coordinate_array
        #no. of dimensions
        self.nd = np.shape(self.co)[1]
        print 'Dimensionality = ',self.nd
        #Check if first and last coordinates are the same. If they are remove 
        #last co-ordinate for list of vertices (.co)
        if np.sum(self.co[-1] != self.co[0]) == 0:
            self.co = self.co[0:-1]
        #list of vertices for plotting, needs first co-ordinate appended on the end
        self.co_plt = np.concatenate(
                                     (
                                      self.co,
                                      np.reshape(self.co[0], [1, self.nd])
                                      ),
                                     axis=0
                                     )
        #note, done this way to avoid a long if expression        
        #no. of walls = #coords in closed shape
        self.n = np.shape(self.co)[0]
        print 'Vertices = ',self.n
        #Periodic boundary flag
        self.pb_flag = np.full(self.n, 0)
        #Index of corresponding opposite boundary
        self.pb_ind = np.full(self.n, -1)
        #Roughness value (randomly rotates wall vector for collision)
        #angle of 1 s.d. note, angle cannot be more than pi/4
        self.rb = np.full(self.n, 0.)
        #Friction value (proportionally reduces parallel v)
        self.fb = np.full(self.n, 0.)
        #wall temperature
        self.T = np.full(self.n, np.nan)
        #heat bath boundary flag
        self.hb_flag = np.full(self.n, 0)
        #heat bath flow speed
        self.hb_v = np.full(self.n, 0.)
        #expected rate /t
        self.pb_t = np.full(self.n, 0.)
        #x and y limits - smallest and largest co-ords in each direction
        self.xlim = np.array(
                             [
                              np.nanmin(self.co[:, 0]),
                              np.nanmax(self.co[:, 0])
                              ]
                             )
        self.ylim = np.array(
                             [
                              np.nanmin(self.co[:, 1]),
                              np.nanmax(self.co[:, 1])
                              ]
                             )
        print 'setting system limits'        
        #get wall vectors and normals.
        print 'finding wall vectors'
        self.vec,self.vlen = self.get_vec(self.co)
        print 'finding wall normal vectors'
        self.norm = self.get_norm(self.vec)
        self.A = self.get_A(self.co)
        print 'System area = ',self.A
        #end initialisation
        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'        
    
    """
    Add extra set of co-ordinates for internal walls
    def add_wall(self,new_co):
        #get length of array
        n = np.shape(new_co)[0]
        new_vec,new_vlen = self.get_vec(new_co)
    """
        
    """
    Calculate time to corner collision
    """
    def t2corn(self, ball):
        p_b = ball.p
        v_b = ball.v
        #Define temporry array for times
        temp = np.full([ball.n,self.n,2],np.nan)
        #Calculate quadratic a,b,c coefficients (similar to 2 ball collision)
        a = np.sum((v_b)**2, axis=1)
        a = np.reshape(a, [ball.n, 1])
        b = np.sum(
                   (np.reshape(p_b, [ball.n, 1, self.nd]) 
                       -np.reshape(self.co, [1, self.n, self.nd])) 
                       *np.reshape(v_b, [ball.n, 1, self.nd]),
                   axis=2
                   )
        b = 2*b
        c = np.sum(
                   (np.reshape(p_b, [ball.n, 1, self.nd])
                       -np.reshape(self.co, [1, self.n, self.nd]))**2
                   ,axis=2
                   )
        c -= np.reshape(ball.r**2,[ball.n,1])
        #Calculate discriminant
        chk = b**2 - 4*a*c
        #Find real solutions
        wh = chk >= 0
        #Find both solutions and take smallest (first collision)
        temp[wh, 0] = (-b[wh]+chk[wh]**0.5)/(2*np.tile(a, [1, self.n])[wh])
        temp[wh, 1] = (-b[wh]-chk[wh]**0.5)/(2*np.tile(a, [1, self.n])[wh])
        temp = np.nanmin(temp, axis=2)
        #If negative solution is in the past, so set as NaN to ignore.
        temp[temp<=0] = np.nan
        #Find smallest positive solution
        t_min = np.nanmin(temp)
        #Set values for ball and corner indexes to use in changing velocities.
        self.i_ball_c = np.where(temp == t_min)[0]
        self.i_corn = np.where(temp == t_min)[1]
        #Return time
        return t_min
        
    """Change velocity of ball after corner collision"""
    def dv_corn(self,ball):
        #Calculate normal vector from ball to vertex
        t_norm = self.co[self.i_corn] -ball.p[self.i_ball_c]
        #Normalise
        t_norm = t_norm/(np.sum(t_norm**2)**0.5)
        """
        THIS CODE NOT YET IMPLEMENTED
        #adjust for roughness value
        if self.rb[self.i_corn] > 0:
            ang = np.pi
            while ang > np.pi/4.:
                ang = np.random.normal(loc=0.,scale=self.rb[self.i_corn])
            t_norm[0] = t_norm[0]*np.cos(ang) -t_norm[1]*np.sin(ang)
            t_norm[1] = t_norm[0]*np.sin(ang) +t_norm[1]*np.cos(ang)
        """
        #Dot product of ball velocity with normal
        dv = np.sum(ball.v[self.i_ball_c] * t_norm)
        #adjust for wall temperature
        if np.isfinite(self.T[self.i_corn]):
            sig = (self.T[self.i_corn]/ball.m[self.i_ball_c])**0.5
            dv_T = np.abs(np.random.normal(loc=0.,scale=sig))
            #Set change in same direction as dv
            if dv <=0:
                dv -=dv_T
            else:
                dv +=dv_T
        #If no T then elastic collision
        else: dv = 2*dv
        #Change ball velcity
        ball.v[self.i_ball_c] -= dv*t_norm
        #Store last dv change
        self.dv = dv
        #Return dv
        return dv
        
    
    """Calculate time to next collision between ball and wall"""
    #new, vectorised and object code
    def t2wall(self, ball):
        p_b = ball.p
        v_b = ball.v
        #Calculate relative normal positions and velocities
        p_rel = np.sum(
                       (np.reshape(self.co, [1,self.n, self.nd]) 
                           -np.reshape(p_b, [ball.n, 1,self.nd]))
                           *np.reshape(self.norm, [1,self.n, self.nd]),
                       axis=2
                       )
        v_rel = np.sum(
                       np.reshape(v_b, [ball.n, 1, self.nd])
                           *np.reshape(self.norm, [1, self.n, self.nd]),
                       axis=2
                       )
        #Time to collision with infinite length wall
        #Subtract time taken for ball to travel one radius to account for size 
        #of ball
        t = p_rel/v_rel -np.abs(np.reshape(ball.r, [ball.n, 1]) /v_rel)
        #Now calculate position along length of wall at collision
        #Parallel position and velocity components
        p_par = np.sum(
                       (np.reshape(p_b, [ball.n, 1, self.nd])
                           -np.reshape(self.co, [1, self.n, self.nd]))
                           *np.reshape(self.vec, [1, self.n, self.nd]),
                       axis=2
                       )
        v_par = np.sum(
                       np.reshape(v_b, [ball.n, 1, self.nd])
                           *np.reshape(self.vec, [1, self.n, self.nd]),
                       axis=2
                       )
        #Position at collision divided by wall length
        p_col = (p_par + v_par*t) /np.reshape(self.vlen, [1, self.n])
        #If >1 or <0 then collision is not along the length of the wall, so
        #ignore
        t[p_col < 0] = np.nan
        t[p_col > 1] = np.nan
        #Negative time is wrong direction so also ignore
        t[t <= 0] = np.nan
        #Calculate time to next collision
        t_min = np.nanmin(t)
        #Find indices of ball and wall section and save for dv calculation
        [i_ball, i_wall] = np.where(t == t_min)
        self.i_ball = i_ball
        self.i_wall = i_wall
        #return next time
        return t_min
        
    """
    Check if a ball is inside the system:
        Used for initialisation of ball object. Works on the principle that a 
        particle inside the system will pass through an odd number of walls on 
        the way out, whereas a ball outside will pass through an odd number
    """
    def isinside(self,p_in,v_in,r_in):
        #change dimension to 2d if ball is greater than 2d
        if p_in.size > self.nd:
            #position
            p_b = np.copy(p_in[0:self.nd])
            p_b[-1] = (np.sum((p_in[(self.nd-1):])**2))**0.5
            #velocity
            v_b = np.copy(v_in[0:self.nd])
            v_b[-1] = (np.sum((v_in[(self.nd-1):])**2))**0.5
        else:
            p_b = p_in
            v_b = v_in
        #Calculate relative position and velocity along normal
        p_rel = np.sum(
                       (self.co - np.reshape(p_b, [1, self.nd]))
                           *np.reshape(self.norm, [self.n, self.nd]),
                       axis=1
                       )
        v_rel = np.sum(
                       np.reshape(v_b, [1, self.nd])
                           *np.reshape(self.norm, [self.n, self.nd]),
                       axis=1
                       )     
        #Calculate times to collsion from relative p/v minus the time to travel
        #one ball radius
        t = p_rel/v_rel -np.abs(r_in/v_rel)
        #Calculate parallel velcoity and position to calculate distance along
        #wall at collision
        p_par = np.sum(
                       (np.reshape(p_b, [1, self.nd]) - self.co)
                           *np.reshape(self.vec, [self.n, self.nd]),
                       axis=1
                       )
        v_par = np.sum(
                       np.reshape(v_b,[1,self.nd])
                           *np.reshape(self.vec, [self.n, self.nd]),
                       axis=1
                       )
        #Relative position along wall
        p_col = (p_par + v_par*t) /self.vlen
        #If <0 or >1 is out of bounds of wall
        t[p_col < 0] = np.nan
        t[p_col > 1] = np.nan
        #If -ve behind ball so ignore
        t[t < 0] = np.nan
        #Calculate check of number of times to collision
        chk = np.sum(np.isfinite(t)) % 2
        #1 if inside, 0 if outside
        if chk == 1:
            #second check in opposite direction (-v_rel)
            t2 = p_rel /(-v_rel) -np.abs(r_in/v_rel)
            p_col2 = (p_par-v_par*t2) /self.vlen
            t2[p_col2 < 0] = np.nan
            t2[p_col2 > 1] = np.nan
            t2[t2 < 0] = np.nan
            chk = np.sum(np.isfinite(t2)) % 2
        return chk
        
    """Change ball velocity after collision with wall"""
    def dv_wall(self, ball):
        v_b = ball.v[self.i_ball]
        #set temporary wall normal
        t_norm = self.norm[self.i_wall]
        if self.pb_ind[self.i_wall] >= 0:
            i_pb = self.pb_ind[self.i_wall]
            i_pb = i_pb.astype('int')
            p_b = ball.p[self.i_ball]
            p_par = np.sum((p_b-self.co[self.i_wall]) *self.vec[self.i_wall])
            p_par = p_par/self.vlen[self.i_wall]
            new_par = 1-p_par
            new_p = self.co[i_pb] + self.vec[i_pb]*new_par*self.vlen[i_pb]
            p_rel = np.sum((p_b-self.co[self.i_wall]) *self.norm[self.i_wall])
            new_p += p_rel*self.norm[i_pb]
            ball.p[self.i_ball] = new_p
            self.dv = 0
            return 0
        """
        #adjust for roughness value
        if self.rb[self.i_wall] > 0:
            ang = np.pi
            while ang > np.pi/4.:
                ang = np.random.normal(loc=0.,scale=self.rb[self.i_wall])
            t_norm[0] = t_norm[0]*np.cos(ang) -t_norm[1]*np.sin(ang)
            t_norm[1] = t_norm[0]*np.sin(ang) +t_norm[1]*np.cos(ang)
        """
        dv = np.sum(v_b*t_norm)
        #adjust for wall temperature
        if np.isfinite(self.T[self.i_wall]):
            #WHY 2 TIMES?
            dv_T = np.random.normal(
                                    loc=0.,
                                    scale=(2*self.T[self.i_wall]
                                           /ball.m[self.i_ball])**0.5
                                    )
            dv_T = np.abs(dv_T)
            if dv >= 0:
                dv += dv_T
            else:
                dv -= dv_T
        else: 
            dv = 2*dv
        self.dv = dv
        ball.v[self.i_ball] = ball.v[self.i_ball] - dv*t_norm
        #self.P[self.i_wall] += ball.m[self.i_ball]*abs(dv)
        return dv
        
class balls:
    #generate balls
    def __init__(self,n_in,r_in,m_in,nd,T_in,v_in,wall):
        print '>>>Initialising particle properties'
        self.nd = nd
        print 'Dimensionality = ',self.nd
        self.n = np.sum(n_in)
        print 'Number of particles = ',self.n
        self.T = T_in
        print 'Initial temperature = ',self.T*120
        self.v_tot = v_in
        print 'Initial bulk velocity = ',self.v_tot
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
        self.p = np.full([self.n, self.nd], np.nan)
        self.v = np.full([self.n, self.nd], np.nan)  
        j = 0
        cnt = 0
        while j < self.n:
            cnt +=1
            if self.nd == 2:
                new_p = np.array(
                                 [
                                  np.random.uniform(
                                      low=wall.xlim[0] +self.r[j],
                                      high=wall.xlim[1] -self.r[j]),
                                  np.random.uniform(
                                      low=wall.ylim[0] +self.r[j],
                                      high=wall.ylim[1] -self.r[j])
                                  ]
                                 )
                new_v = np.random.normal(
                                         loc=self.v_tot,
                                         scale=(self.T/self.m[j])**0.5,
                                         size=[1,2]
                                         )
                #new_p = np.reshape(new_p,[1,2])
            elif self.nd == 3:
                new_p = np.array(
                                 [
                                  np.random.uniform(
                                      low=wall.xlim[0] +self.r[j], 
                                      high=wall.xlim[1] -self.r[j]),
                                  np.random.uniform(
                                      low=wall.ylim[0] +self.r[j], 
                                      high=wall.ylim[1] -self.r[j]),
                                  np.random.uniform(
                                      low=wall.ylim[0] +self.r[j], 
                                      high=wall.ylim[1]-self.r[j])
                                  ]
                                 )
                new_p = np.reshape(new_p, [1, 3])
                new_v = np.random.normal(
                                         loc=self.v_tot,
                                         scale=(self.T/self.m[j])**0.5,
                                         size=[1,3]
                                         )
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
        print 'Particles initialised: ',j
        print 'Attempts: ',cnt                   
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

    def get_v_tot(self):
        self.v_tot = np.sum(self.v,axis=0)
        
    def get_mv_tot(self):
        self.mv_tot = np.median(self.v*np.reshape(self.m,[self.n,1]),axis=0)
        return self.mv_tot

    def get_T(self):
        self.T = np.mean(
                         np.std(
                                self.v *np.reshape(self.m**0.5, [self.n, 1]),
                                axis=0
                                )**2
                         )
        return self.T
