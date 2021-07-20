# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 08:49:55 2019

@author: Will Arnuk
"""

""" 3D Plotting

This should:
    - iterate through output at each co2 (y) (c or o)
    - compile list of concentrations (z) at each temp (x) for which
        a certain mineral precipitates 
    - plot xyz as 3D scatterplot or surface?
    
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def file_detect(directory, filetype):
    
    """ This function iterates through all files in a directory, returning
    a list of files that match a specified ending. """
    
    files = []
    for file in os.listdir(directory):
        if file.endswith(filetype):
            #print(os.path.join("/mydir", file))
            files.append(file)
    return(files)

def compile_xyz(simulation_name, system, mineral, first_last, var='sal'):
    
    base = os.path.join(simulation_name, "closed" if system=="c" else "open")
    
    directories = os.listdir(base)
    
    X, Y, Z = [], [], []
    
    for directory in directories:
        if not directory.startswith("."):
            data_dir = os.path.join(base, directory, "output")
            files = file_detect(data_dir, ".j{}%".format(system))
        
            for file in files:
                filename = file.split(".j{}%".format(system))[0]
                water, system, pco2, temp = filename.split("_")
                temp = float(temp.split("C")[0])
                pco2 = float(pco2.split("ppm")[0])
                
                df = pd.read_csv(os.path.join(data_dir, file))
                df[var] = pd.read_csv(os.path.join(data_dir, filename + ".j{}&".format(system)))[var]
                target = np.nan
                indexer = np.ones(df.shape[0])
                for row in range(0, df.shape[0]):
                    indexer[row] = indexer[row] * (mineral in 
                                                   df.columns.values[
                                                       np.nonzero(
                                                           df.iloc[row,:].values)][1:])
                subset = df.loc[indexer == 1]
                subset.shape[0]
                if subset.shape[0] > 0:
                    if first_last == 'first':
                        target = subset.loc[subset['fc'] == subset['fc'].min()][var].values[0]
                    elif first_last == 'last':
                        target = subset.loc[subset['fc'] == subset['fc'].max()][var].values[0]
                    else:
                        raise ValueError('invalid entry for argument \'first_last\'. enter \'first\' or \'last\'')
                
                X.append(temp)
                Y.append(pco2)
                Z.append(target)
    
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    x, y = np.unique(X), np.unique(Y)
    x.sort(), y.sort()
    x, y = np.meshgrid(x, y)
    z = np.zeros(y.shape)

    for i in range(0, z.shape[0]):
        for k in range(0, z.shape[1]):
            z[i][k] = Z[(X == x[i][k]) & (Y == y[i][k])]
    
    return(np.array([x, y, z]))
            
def explode_xyz(arr):
    
    x = []
    y = []
    z = []
     
    for i in range(0, arr[0].shape[0]):
        
        for k in range(0, arr[0].shape[1]):
            
            x.append(arr[0][i][k])
            y.append(arr[1][i][k])
            z.append(arr[2][i][k])
        
    return(np.array([x, y, z]).T)

def kernel_mean(array):
    smooth = np.ones(array.shape)
    smooth[0:2,:,:] = array[0:2,:,:]
    
    for i in range(0,smooth.shape[1]):
        for k in range(0,smooth.shape[2]):
            low_x, low_y = i-1, k-1
            high_x, high_y = i+2, k+2
            if i == 0:
                low_x = 0
            elif i == smooth.shape[1]-1:
                high_x = smooth.shape[1]+1
            if k == 0:
                low_y = 0
            elif k == smooth.shape[2]-1:
                high_y = smooth.shape[2]+1
            
            subset = array[2,low_x:high_x,low_y:high_y]
            if (subset != np.nan).any():
                smooth[2,i,k] = np.nanmean(subset)
            else:
                smooth[2,i,k] = np.nan
    return(smooth)

# Change these
load_existing = False
simulation_name = "steamboat-no-dol_071921"
system = 'c'

# Set the path to the simulation 3D arrays folder
directory = os.path.join("3d_arrays", simulation_name, system)

# load 3D arrays if True, make 3D arrays if False
if load_existing:
    trona = np.load(os.path.join(directory, "trona.npy"))
    halite = np.load(os.path.join(directory, "halite.npy"))
    calcite = np.load(os.path.join(directory, "calcite.npy"))
    gaylussite = np.load(os.path.join(directory, "gaylussite.npy"))
    pirssonite = np.load(os.path.join(directory, "pirssonite.npy"))
    northupite = np.load(os.path.join(directory, "northupite.npy"))
    nahcolite = np.load(os.path.join(directory, "nahcolite.npy"))
    natron = np.load(os.path.join(directory, "natron.npy"))
    
else:
    if not os.path.exists(directory) and not os.path.isdir(directory):
        os.makedirs(directory)

    trona = compile_xyz(simulation_name, system, "TRONA", "first", "sal")
    np.save(os.path.join(directory, "trona"), trona)

    halite = compile_xyz(simulation_name, system, "HALITE", "first", "sal")
    np.save(os.path.join(directory, "halite"), halite)

    calcite = compile_xyz(simulation_name, system, "CALCITE", "first", "sal")
    np.save(os.path.join(directory, "calcite"), calcite)

    gaylussite = compile_xyz(simulation_name, system, "GAYLUSSITE", "first", "sal")
    np.save(os.path.join(directory, "gaylussite"), gaylussite)

    pirssonite = compile_xyz(simulation_name, system, "PIRSSONITE", "first", "sal")
    np.save(os.path.join(directory, "pirssonite"), pirssonite)
    
    northupite = compile_xyz(simulation_name, system, "NORTHUPITE", "first", "sal")
    np.save(os.path.join(directory, "northupite"), northupite)

    nahcolite = compile_xyz(simulation_name, system, "NAHCOLITE", "first", "sal")
    np.save(os.path.join(directory, "nahcolite"), nahcolite)

    natron = compile_xyz(simulation_name, system, "NATRON", "first", "sal")
    np.save(os.path.join(directory, "natron"), natron)

halite = kernel_mean(halite)
trona = kernel_mean(trona)
gaylussite = kernel_mean(gaylussite)
northupite = kernel_mean(northupite)
pirssonite = kernel_mean(pirssonite)
calcite = kernel_mean(calcite)

# Make plot figure and 3D axis
fig = plt.figure()
ax = plt.axes(projection='3d')

h = ax.plot_surface(halite[0], halite[1], halite[2], color='pink', alpha=0.75,label='Halite')
t = ax.plot_surface(trona[0], trona[1], trona[2], color='black', alpha=0.75,label='Trona')
g = ax.plot_surface(gaylussite[0], gaylussite[1], gaylussite[2], color='yellow', alpha=0.75, label='Gaylussite')
p = ax.plot_surface(pirssonite[0], pirssonite[1], pirssonite[2], color='green', alpha=0.75, label='Pirssonite')
cc = ax.plot_surface(calcite[0], calcite[1], calcite[2], color='blue',alpha=0.75, label='Calcite')
no = ax.plot_surface(northupite[0], northupite[1], northupite[2], color='gray',alpha=0.75, label='Northupite')
#nah = ax.plot_surface(nahcolite[0], nahcolite[1], nahcolite[2], color='red',alpha=0.75, label='Nahcolite')
#nat = ax.plot_surface(natron[0], natron[1], natron[2], color='brown',alpha=0.75, label='Natron')

plt.show()

