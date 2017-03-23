# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 12:17:39 2017

@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt
import brownian_classes as brw
import brownian_system as brs

#co-ords fo benoulli tube shape
in_co = np.array([[0,0],
                  [0.2,0],
                  [0.4,0.2],
                  [0.6,0.2],
                  [0.8,0],
                  [1,0],
                  [1,1],
                  [0.8,1],
                  [0.6,0.8],
                  [0.4,0.8],
                  [0.2,1],
                  [0,1]])

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

#box to test
in_co = np.array([[0,0],
                  [0,10],
                  [10,10],
                  [10,0]])


wal = brw.wall_shape(in_co)
#wal.T[:] = 2.5
wal.pb_ind = np.array([2,3,0,1])
"""
bal = brw.balls(50,0.1,1,2,2.5,0.,wal)
sys = brs.system(wal,bal)
sys.run_plt(100000)
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