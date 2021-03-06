"""
Created on Wed Jan 25 09:50:23 2017

@author: wj
History:
    16/02/2017, LM: Small adjustments so can run properly with animation file
    07/03/2017, WJ: General tidying
"""
#
import numpy as np
import matplotlib.pyplot as plt
#Define system class to hold wall,ball classes and system procedures

class system:
    def __init__(self,wall_in,ball_in):
        print '>>>Initialising System'
        print 'Inputing wall properties'
        self.wall = wall_in
        print 'Inputing particle properties'
        self.ball = ball_in
        print 'Setting system time'
        self.t = 0
        self.t_step = 0.1
        self.t_max = 100
        
        #Incident and reflected velocity tracker
        self.vel = np.full(2, np.nan).reshape([2,-1])
        
        self.b = 0
        self.c = 0
        self.w = 0
        self.type = 0
        
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
                self.b +=1
                self.type = 1

            if t_2wall == self.min_t:
                self.wall.dv_wall(self.ball)
                self.w +=1
                self.type = 2
                self.vel = np.concatenate((self.vel,self.wall.velw), axis=1)

            if t_2corn == self.min_t:
                self.wall.dv_corn(self.ball)
                self.c +=1
                self.type = 3
                self.vel = np.concatenate((self.vel,self.wall.velc), axis=1)

        else:
            self.t +=self.t_step
            self.ball.p += self.ball.v*self.t_step
            self.type = 0
        self.ball.get_T()
            
        return self.ball.p, self.t, self.ball.T
            
    def plt_sys(self):
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                             xlim=(self.wall.xlim[0]-0.1,self.wall.xlim[1]+0.1),
                             ylim=(self.wall.ylim[0]-0.1,self.wall.ylim[1]+0.1))
        wal, = ax.plot(self.wall.co_plt[:,0],self.wall.co_plt[:,1])
        bal, = ax.plot(self.ball.p[:,0],self.ball.p[:,1],'bo',
                       ms=fig.dpi
                       *fig.get_figwidth()/(ax.get_xlim()[1]-ax.get_xlim()[0])
                       *2*self.ball.r[0]
                       )
        fig.show
    
    def run(self,n_step):
        print 'Running system for ',n_step,' steps'
        #Initialise arrays for recording system properties
        # elpased time
        self.t_el =  np.full(n_step+1,0)
        # ball index
        self.i_b = np.full([n_step+1,2],-1)
        #wall index
        self.i_w = np.full([n_step+1,1],-1)
        #collision type
        self.typ = np.full(n_step+1,0)
        # collision velocities
        self.dv = np.full([n_step+1,2],0)
        # mean velocity
        self.v_av = np.full([n_step+1,2],0)
        self.v_av[0] = self.ball.get_v_av()
        # total momentum
        self.mv = np.full([n_step+1,2],0)
        self.mv[0] = self.ball.get_mv_tot()
        # system temperature
        self.T = np.full(n_step+1,0)
        self.T[0] = self.ball.get_T()
        # total energy
        self.E = np.full(n_step+1,0)
        self.E[0] = self.ball.get_E_tot()
        # kinetic energy
        self.KE = np.full(n_step+1,0)
        self.KE[0] = self.ball.get_E_KE()
        # heat
        self.Q = np.full(n_step+1,0)
        self.Q[0] = self.ball.get_E_Q()
        # angular momentum
        self.mv_ang = np.full(n_step+1,0)
        self.mv_ang[0] = self.ball.get_mv_ang()
        
        for i in range(n_step):
            # print progress
            if i > 0 and i % (n_step/10) == 0:
                print i*100./n_step,' % complete'
            # step forwards
            self.step()
            # record system properties
            self.t_el[i+1] = self.t
            self.typ[i+1] = self.type
            if self.type == 1:
                self.i_b[i+1,0] = self.ball.i_ball
                self.i_b[i+1,1] = self.ball.j_ball
                self.dv[i+1] = [self.ball.dv1, self.ball.dv2]
            elif self.type == 2:
                self.i_w[i+1] = self.wall.i_wall
                self.i_b[i+1,0] = self.wall.i_ball
                self.dv[i+1] = [self.wall.dv_in, self.wall.dv_out]
            elif self.type == 3:
                self.i_w[i+1] = self.wall.i_corn
                self.i_b[i+1,0] = self.wall.i_ball_c
                self.dv[i+1] = [self.wall.dv_in, self.wall.dv_out]
            #update particle properties
            self.v_av[i+1] = self.ball.get_v_av()
            self.mv[i+1] = self.ball.get_mv_tot()
            self.T[i+1] = self.ball.get_T()
            self.E[i+1] = self.ball.get_E_tot()
            self.KE[i+1] = self.ball.get_E_KE()
            self.Q[i+1] = self.ball.get_E_Q()
            self.mv_ang[i+1] = self.ball.get_mv_ang()
            #print np.shape(np.sum(self.ball.v**2,axis=1)**0.5)
            #print self.v_av[i+1]
            
        # collision histogram code
        #self.vel = np.delete(self.vel,0,1)
        #plt.hist(self.vel[0],bins=30)
        #plt.show()
        #plt.hist(self.vel[1],bins=30)
        #plt.show()
        print 'System time elapsed: ',self.t
        print 'Particle collisions: ',self.b
        print 'Wall collisions: ',self.w
        print 'Corner collisions: ',self.c
        v_av_all = np.mean(self.v_av)
        print 'Average particle velocity = ', v_av_all
        v_err_all = np.std(self.v_av)/n_step**0.5
        print 'v error = ', v_err_all
        
        #Calculate mean free path and collision time
        t_col_tot = 0
        t_col_n = 0
        i_b_col = self.i_b[self.typ == 1]
        t_col = self.t_el[self.typ == 1]
        for j in np.arange(self.ball.n):
            n_t = np.sum(i_b_col == j)
            if n_t >= 2:
                i = np.where(i_b_col == j)[0]
                i = np.sort(i)
                t_i = t_col[i]
                t_col_tot += t_i[-1]-t_i[0]
                t_col_n += (n_t-1)
        self.t_col = t_col_tot/t_col_n
        self.t_col_err = self.t_col/t_col_n**0.5
        print 'Average collision time = ',self.t_col
        print 'Error = ',self.t_col_err
        self.d_col = self.t_col*v_av_all
        self.d_col_err = self.d_col/t_col_n**0.5
        print 'Mean free path = ', self.d_col
        print 'Error = ',self.d_col_err
        #print 'Wall collision velocities = ', self.vel_distrib
        return self.t_el, \
               self.t_col, \
               self.t_col_err, \
               self.d_col, \
               self.d_col_err, \
               self.v_av, \
               self.mv, \
               self.T, \
               self.E, \
               self.KE, \
               self.Q, \
               self.mv_ang
        
        
        
        
    """
     Old run and plot code, use .run and .plt_... functions instead 
    """ 
    def run_plt(self,n_step):
        print 'Running system for ',n_step,' steps'
        
        print 'Recording momentum, angular momentum and temperature'
        self.T = np.full([n_step+1],np.nan)
        self.T[0] = np.mean(np.std(self.ball.v,axis=0)**2)
        
        self.mv_tot = np.full([n_step+1,self.ball.nd],np.nan)
        self.mv_tot[0] = np.sum(self.ball.v
                                *np.reshape(self.ball.m,[self.ball.n,1])
                                ,axis=0)
        
        v_tot = np.sum(self.ball.v,axis=0)
        #print v_tot
        
        self.KE = np.full([n_step+1,self.ball.nd],np.nan)
        self.KE[0] = np.sum(self.ball.m)*np.sum(v_tot**2)/2
        #print self.KE[0]
        
        self.Q = np.full([n_step+1,self.ball.nd],np.nan)
        self.Q[0] = np.sum(np.reshape(self.ball.m,[self.ball.n,1])
                           *(self.ball.v -np.reshape(v_tot,[1,self.ball.nd]))**2)/2
        #print self.Q[0]
        
        self.mv_ang = np.full([n_step+1],np.nan)
        self.mv_ang[0] = np.sum(self.ball.m
                                *(self.ball.v[:,0]*(self.ball.p[:,1]-0.5)
                                  -self.ball.v[:,1]*(self.ball.p[:,0]-0.5)))
        
        self.P = np.copy(self.wall.P)
        self.t_el =  np.full(n_step+1,0)
        
        self.dv = np.full([n_step+1],np.nan)
        self.dv[0] = 0

        self.dv_i = np.full([n_step+1],self.wall.n+1)
        
        self.i_b = np.full([n_step+1,2],-1)
        
        
        for i in range(n_step):
            self.step()
            self.t_el[i+1] = self.t
            self.ball.get_T()
            self.T[i+1] = self.ball.T

            self.ball.get_mv_tot()
            self.mv_tot[i+1] = self.ball.mv_tot
            
            self.mv_ang[i+1] = np.sum(self.ball.m
                                *(self.ball.v[:,0]*(self.ball.p[:,1]-0.5)
                                  -self.ball.v[:,1]*(self.ball.p[:,0]-0.5)))
            
            self.P = (self.wall.P/self.wall.vlen)/self.t
            

            """
            v_tot = np.sum(self.ball.v,axis=0)
            
            self.KE[i+1] = np.sum(self.ball.m)*np.sum(v_tot**2)/2
            self.Q[i+1] = np.sum(np.reshape(self.ball.m,[self.ball.n,1])
                               *(self.ball.v -np.reshape(v_tot,[1,self.ball.nd]))**2)/2
            """
            
            if self.type == 2:
                self.dv[i+1] = abs(self.wall.dv)
                self.dv_i[i+1] = self.wall.i_wall

            if self.type == 1:
                self.i_b[i+1,0] = self.ball.i_ball
                self.i_b[i+1,1] = self.ball.j_ball
            
            if i % 2000 == 0:
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
                                     xlim=(self.wall.xlim[0]-0.1,self.wall.xlim[1]+0.1),
                                     ylim=(self.wall.ylim[0]-0.1,self.wall.ylim[1]+0.1))
                wal, = ax.plot(self.wall.co_plt[:,0],self.wall.co_plt[:,1])
                bal, = ax.plot(self.ball.p[:,0],self.ball.p[:,1],'bo',
                               ms=fig.dpi
                                   *fig.get_figwidth()
                                   /(ax.get_xlim()[1]-ax.get_xlim()[0])
                                   *2*self.ball.r[0]
                               )
                fig.show
                
        print 'Particle collisions: ',self.b
        print 'Wall collisions: ',self.w
        print 'Corner collisions: ',self.c
        
        print 'Plotting temperature and momentum'
        fig2 = plt.figure(figsize=(6,6))
        ax2 = fig2.add_subplot(111)
        ax2.plot(self.t_el,self.T*120)
        ax2.set_title('System Temperature')
        ax2.set_xlabel('Elapsed time /ps')
        ax2.set_ylabel('Temperature /K')
        fig3 = plt.figure(figsize=(6,6))
        ax3 = fig3.add_subplot(111)
        ax3.plot(self.t_el,self.mv_tot)
        ax3.set_title('Linear Momentum')
        ax3.set_xlabel('Elapsed time /ps')
        ax3.set_ylabel('Momentum /AMU.m/s')
        ax3.legend(['x component','y component'])
        
        fig4 = plt.figure(figsize=(6,6))
        ax4 = fig4.add_subplot(111)
        ax4.plot(self.t_el,self.mv_ang)
        ax4.set_xlabel('Elapsed time /ps')
        ax4.set_ylabel('Angular Momentum /AMU/s')
        """
        ax5 = fig3.add_subplot(224)
        ax5.plot(self.t_el,self.P)
        """
        
        fig2.show
        fig3.show
        fig4.show
        """
        fig5 = plt.figure(figsize=(6,6))
        ax5 = fig5.add_subplot(111)
        plot_dv = self.dv[self.dv_i < 8]
        #print plot_dv
        plot_dv_i = self.dv_i[self.dv_i < 8]
        plot_dv_i = plot_dv_i.astype(int)
        #print plot_dv_i
        plot_dv = plot_dv * self.wall.norm[plot_dv_i,0]
        #print plot_dv
        ax5.hist(plot_dv, 50, normed=1, facecolor='green', alpha=0.75)
        fig5.show
        """
        t_col_tot = 0
        n_tot = 0
        for j in np.arange(self.ball.n):
            n_col = np.sum(self.i_b == j)
            if n_col >= 2:
                i = np.where(self.i_b == j)[0]
                i = np.sort(i)
                t_i = self.t_el[i]
                t_col_tot += np.mean(t_i[-1]-t_i[0])
                n_tot += (n_col-1)
        self.t_col = t_col_tot/n_tot
        print self.t_col
        """
        fig4 = plt.figure(figsize=(6,6))
        ax6 = fig4.add_subplot(111)
        ax6.plot(self.t_el,self.KE)
        ax6.plot(self.t_el,self.Q)
        fig3.show
        """
        print 'System time elapsed: ',self.t
            
        
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
