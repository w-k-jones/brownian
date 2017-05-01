# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 12:17:39 2017

@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt
import brownian_classes as brw
import brownian_system as brs
"""
#co-ords fo benoulli tube shape
in_co = np.array([[0,0],
                  [0.2,0],
                  [0.4,0.3],
                  [0.6,0.3],
                  [0.8,0],
                  [1,0],
                  [1,1],
                  [0.8,1],
                  [0.6,0.7],
                  [0.4,0.7],
                  [0.2,1],
                  [0,1]])*10

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
"""
#box to test
in_co = np.array([[0,0],
                  [0,10],
                  [10,10],
                  [10,0]])


wal = brw.wall_shape(in_co)
#wal.pb_ind = np.array([2,3,0,1])
wal.T[:] = 2.5
#wal.pb_ind[5] = 11

bal = brw.balls(50,0.1,1,2,2.5,0.,wal)
sys = brs.system(wal,bal)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(sys.wall.xlim[0]-0.1,sys.wall.xlim[1]+0.1),
                     ylim=(sys.wall.ylim[0]-0.1,sys.wall.ylim[1]+0.1))
wal, = ax.plot(sys.wall.co_plt[:,0],sys.wall.co_plt[:,1])
bal, = ax.plot(sys.ball.p[:,0],sys.ball.p[:,1],'bo',
               ms=fig.dpi
               *fig.get_figwidth()/(ax.get_xlim()[1]-ax.get_xlim()[0])
               *2*sys.ball.r[0]
               )
out = sys.run(10000)
mv = out[6]
E = out[8]
T = out[7]
t_el = out[0]
fig_1 = plt.figure(figsize=(6,6))
ax_1 = fig_1.add_subplot(111)
fig_2 = plt.figure(figsize=(6,6))
ax_2 = fig_2.add_subplot(111)
fig_3 = plt.figure(figsize=(6,6))
ax_3 = fig_3.add_subplot(111)

ax_1.plot(t_el,mv[:,0],'b-')
ax_2.plot(t_el,E,'g-')
ax_3.plot(t_el,T*120,'r-')
ax_1.set_title('System Total Momentum')
ax_1.set_xlabel('Elapsed time $/ps$')
ax_1.set_ylabel('Momentum $/AMU.nm.ps^{-1}$')


ax_2.set_title('System Total Energy')
ax_2.set_xlabel('Elapsed time /ps')
ax_2.set_ylabel('Energy /AMU.m^2/s^2')


ax_3.set_title('System Temperature')
ax_3.set_xlabel('Elapsed time /ps')
ax_3.set_ylabel('Temperature /K')

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
wall_p = sys.dv[:,0]
wall_p = wall_p[sys.typ == 2]
ax.hist(wall_p,40)
T_mean = np.mean(T)
ax.plot(np.arange(0,6,0.1), sys.w/7*np.arange(0,6,0.1)/T_mean*np.exp(-np.arange(0,6,0.1)**2*0.5/T_mean),'r--',linewidth=2)
ax.set_title('Wall Impulse Distribution')
ax.set_xlabel('Impulse $/AMU.nm.ps^{-1}$')
ax.set_ylabel('Frequency')
"""
t_col = np.full([10],np.nan)
t_col_err = np.full([10],np.nan)
d_col = np.full([10],np.nan)
d_col_err = np.full([10],np.nan)
for i in np.arange(10):
    #for j in np.arange(10):
    bal = brw.balls(25,0.05*(i+1),1,2,2.5,0.,wal)
    sys = brs.system(wal,bal)
    out = sys.run(5000)
    t_col[i] = out[1]
    t_col_err[i] = out[2]
    d_col[i] = out[3]
    d_col_err[i] = out[4]
    

print t_col
print d_col

fig1 = plt.figure(figsize=(6,6))
ax1 = fig1.add_subplot(111)
ax1.plot((np.arange(10)+1)*0.05,t_col,'b-', linewidth=2)
ax1.plot((np.arange(10)+1)*0.05,t_col+t_col_err,'b--', linewidth=2)
ax1.plot((np.arange(10)+1)*0.05,t_col-t_col_err,'b--', linewidth=2)
ax1.set_title('Mean Collision Time vs Particle Radius')
ax1.set_xlabel('Particle Radius $/nm$')
ax1.set_ylabel('Mean Collision Time $/ps$')

fig2 = plt.figure(figsize=(6,6))
ax2 = fig2.add_subplot(111)
ax2.plot((np.arange(10)+1)*0.05,d_col,'g-', linewidth=2)
ax2.plot((np.arange(10)+1)*0.05,d_col+d_col_err,'g--', linewidth=2)
ax2.plot((np.arange(10)+1)*0.05,d_col-d_col_err,'g--', linewidth=2)
ax2.plot((np.arange(10)+1)*0.05,1/(2**0.5*4*0.25*(np.arange(10)+1)*0.05),'k:', linewidth=3)
ax2.set_title('Mean Free Path vs Particle Radius')
ax2.set_xlabel('Particle Radius $/nm$')
ax2.set_ylabel('Mean Free Path $/nm$')
"""