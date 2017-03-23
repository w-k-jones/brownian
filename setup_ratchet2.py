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
                  [ 20.       ,   0.       ],
                  [21,0],
                  [21,2],
                  [-1,2],
                  [-1,0]]
                 )

wal = brw.wall_shape(in_co)
wal.T[:] = 1
#wal.T[:20] = 100

#wal.pb_ind[21] = 23
#wal.pb_ind[23] = 21
"""
in_co2 = np.array([[0,0],
                  [10,0],
                  [10,10],
                  [0,10]])
wal = brw.wall_shape(in_co2)
"""
#wal.pb_ind = np.array([2,3,0,1])
#bal = brw.balls(107,0.1,1,2,1.,0.,wal)

n_run = 100
l_run = 1000

fig_1 = plt.figure(figsize=(6,6))
ax_1 = fig_1.add_subplot(111)
fig_2 = plt.figure(figsize=(6,6))
ax_2 = fig_2.add_subplot(111)
fig_3 = plt.figure(figsize=(6,6))
ax_3 = fig_3.add_subplot(111)

mv_arr = np.full([n_run,l_run+1],0.)

for i in np.arange(n_run):
    wal = brw.wall_shape(in_co)
    #wal.T[:] = 1.
    #wal.T[0:20] = 100

    wal.pb_ind[21] = 23
    wal.pb_ind[23] = 21
    #wal.pb_ind = np.array([2,3,0,1])
    bal = brw.balls(20,0.1,1,2,1.,[0.84,0.],wal)
    sys = brs.system(wal,bal)
    """
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
    """
    out = sys.run(l_run)
    mv = out[6]
    mv_arr[i] = mv[:,0]
    E = out[8]
    T = out[7]
    t_el = out[0]
    ax_1.plot(t_el,mv[:,0],'b-',alpha=1/n_run**0.75)
    ax_2.plot(t_el,E,'g-',alpha=1/n_run**0.75)
    ax_3.plot(t_el,T*120,'r-',alpha=1/n_run**0.75)
    if i == 0:
        mv_all = out[6]
        E_all = out[8]
        T_all = out[7]
        t_el_all = out[0]
    else:
        mv_all += out[6]
        E_all += out[8]
        T_all += out[7]
        t_el_all += out[0]

mv_err = np.std(mv_arr,axis=0)

ax_1.plot(t_el_all/n_run,mv_all[:,0]/n_run,'b--',linewidth=2)
ax_1.plot(t_el_all/n_run,mv_all[:,0]/n_run+mv_err,'b:',linewidth=2)
ax_1.plot(t_el_all/n_run,mv_all[:,0]/n_run-mv_err,'b:',linewidth=2)
ax_1.plot(t_el_all/n_run,np.full(l_run+1,0),'k--',linewidth=2)
ax_1.set_title('System Total X Momentum')
ax_1.set_xlabel('Elapsed time $/ps$')
ax_1.set_ylabel('Momentum $/AMU.nm.ps^{-1}$')

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
