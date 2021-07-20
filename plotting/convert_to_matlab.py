# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:02:51 2020

@author: warnuk
"""

import numpy as np
import scipy.io

trona = np.load('3d_arrays/trona.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\trona.mat', dict(x=trona[0], y=trona[1], z=trona[2]))

halite = np.load('3d_arrays/halite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\halite.mat', dict(x=halite[0], y=halite[1], z=halite[2]))


gaylussite = np.load('3d_arrays/gaylussite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\gaylussite.mat', dict(x=gaylussite[0], y=gaylussite[1], z=gaylussite[2]))


nahcolite = np.load('3d_arrays/nahcolite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\nahcolite.mat', dict(x=nahcolite[0], y=nahcolite[1], z=nahcolite[2]))


pirssonite = np.load('3d_arrays/pirssonite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\pirssonite.mat', dict(x=pirssonite[0], y=pirssonite[1], z=pirssonite[2]))


calcite = np.load('3d_arrays/calcite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\calcite.mat', dict(x=calcite[0], y=calcite[1], z=calcite[2]))

dolomite = np.load('3d_arrays/dolomite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\dolomite.mat', dict(x=dolomite[0], y=dolomite[1], z=dolomite[2]))


northupite = np.load('3d_arrays/northupite.npy')
scipy.io.savemat(r'C:\Users\warnuk\Documents\MATLAB\northupite.mat', dict(x=northupite[0], y=northupite[1], z=northupite[2]))