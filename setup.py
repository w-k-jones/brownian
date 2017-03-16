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

t_col = np.full([20],np.nan)
d_col = np.full([20],np.nan)
for i in np.arange(20):
    #for j in np.arange(10):
    bal = brw.balls((i+1)*10,0.1,1,2,2.5,0.,wal)
    sys = brs.system(wal,bal)
    out = sys.run(1000)
    t_col[i] = out[0]
    d_col[i] = out[1]
    

print t_col
print d_col

fig1 = plt.figure(figsize=(6,6))
ax1 = fig1.add_subplot(111)
ax1.plot(np.arange(20)+1,t_col)

fig2 = plt.figure(figsize=(6,6))
ax2 = fig1.add_subplot(111)
ax2.plot(np.arange(20)+1,d_col)