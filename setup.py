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
                  [0,1],
                  [1,1],
                  [1,0]])


wal = brw.wall_shape(in_co)
wal.T[:] = 2.5
wal.pb_ind[1] = 3
wal.pb_ind[3] = 1
bal = brw.balls(50,0.01,1,2,2.5,0.,wal)
sys = brs.system(wal,bal)

sys.run_plt(10000)