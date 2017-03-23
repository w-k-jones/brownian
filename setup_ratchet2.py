# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 23:17:46 2017

@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt
import brownian_classes as brw
import brownian_system as brs
"""
in_co = np.array([[0,0],
                  [2,1],
                  [2,0],
                  [4,1],
                  [4,0],
                  [6,1],
                  [6,0],
                  [8,1],
                  [8,0],
                  [10,1],
                  [10,0],
                  [15,0],
                  [15,2],
                  [-5,2],
                  [-5,0]])
"""

in_co = np.array([[  0.       ,   0.       ],
                  [  1.8660254,   0.5      ],
                  [  2.       ,   0.       ],
                  [  3.8660254,   0.5      ],
                  [  4.       ,   0.       ],
                  [  5.8660254,   0.5      ],
                  [  6.       ,   0.       ],
                  [  7.8660254,   0.5      ],
                  [  8.       ,   0.       ],
                  [  9.8660254,   0.5      ],
                  [ 10.       ,   0.       ],
                  [ 11.8660254,   0.5      ],
                  [ 12.       ,   0.       ],
                  [ 13.8660254,   0.5      ],
                  [ 14.       ,   0.       ],
                  [ 15.8660254,   0.5      ],
                  [ 16.       ,   0.       ],
                  [ 17.8660254,   0.5      ],
                  [ 18.       ,   0.       ],
                  [ 19.8660254,   0.5      ],
                  [20,0],
                  [21,0],
                  [21,2],
                  [-1,2],
                  [-1,0]]
                 )

wal = brw.wall_shape(in_co)
wal.T[:] = 1
#wal.T[:20] = 100

wal.pb_ind[21] = 23
wal.pb_ind[23] = 21

in_co2 = np.array([[0,0],
                  [10,0],
                  [10,10],
                  [0,10]])
wal = brw.wall_shape(in_co)

#wal.pb_ind = np.array([2,3,0,1])
#bal = brw.balls(107,0.1,1,2,1.,0.,wal)

n_run = 1

fig_1 = plt.figure(figsize=(6,6))
ax_1 = fig_1.add_subplot(111)
fig_2 = plt.figure(figsize=(6,6))
ax_2 = fig_2.add_subplot(111)
fig_3 = plt.figure(figsize=(6,6))
ax_3 = fig_3.add_subplot(111)

for i in np.arange(n_run):
    wal = brw.wall_shape(in_co2)
    #wal.T[:] = 1.
    #wal.T[0:20] = 100

    #wal.pb_ind[21] = 23
    #wal.pb_ind[23] = 21
    wal.pb_ind = np.array([2,3,0,1])
    bal = brw.balls(96,0.1,1,2,2.5,0.,wal)
    sys = brs.system(wal,bal)
    sys.plt_sys
    out = sys.run(1000)
    mv = out[4]
    E = out[6]
    T = out[5]
    t_el = out[0]
    ax_1.plot(t_el,mv[:,0],'b-',alpha=1/n_run**0.75)
    ax_2.plot(t_el,E,'g-',alpha=1/n_run**0.75)
    ax_3.plot(t_el,T*120,'r-',alpha=1/n_run**0.75)
    if i == 0:
        mv_all = out[4]
        E_all = out[6]
        T_all = out[5]
        t_el_all = out[0]
    else:
        mv_all += out[4]
        E_all += out[6]
        T_all += out[5]
        t_el_all += out[0]

ax_1.plot(t_el_all/n_run,mv_all[:,0]/n_run,'b--')
ax_1.set_title('System X Momentum')
ax_1.set_xlabel('Elapsed time /ps')
ax_1.set_ylabel('Momentum /AMU.m/s')

ax_2.plot(t_el_all/n_run,E_all/n_run,'g--')
ax_2.set_title('System Total Energy')
ax_2.set_xlabel('Elapsed time /ps')
ax_2.set_ylabel('Energy /AMU.m^2/s^2')

ax_3.plot(t_el_all/n_run,T_all*120/n_run,'r--')
ax_3.set_title('System Temperature')
ax_3.set_xlabel('Elapsed time /ps')
ax_3.set_ylabel('Temperature /K')

"""
fig_4 = plt.figure(figsize=(6,6))
ax_4 = fig_4.add_subplot(111)
ax_4.hist(sys.ball.v[:,0])
ax_4.hist(sys.ball.v[:,1])
ax_4.hist(np.sum(sys.ball.v**2,axis=1)**0.5)
"""
