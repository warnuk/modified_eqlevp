#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: warnuk
"""
# import the libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import color_mapper
import file_handler

def sequence_plot(infile, outfile, title, colormap, variable):
    # reading the output from the specified trial
    data = pd.read_csv(infile)

    # defining which colors are used to represent each mineral
    if colormap:
        colormap = color_mapper.color_read(colormap)

    cols = data.columns.values.tolist()
    for i in ['fc', 'aH2O', 'ionStr', 'salinity']:
        cols.remove(i)
    y = {}

    for mineral in cols:
        y[mineral] = data[mineral].values

    negligible = []

    for mineral in y:
        if y[mineral].sum() == 0:
            negligible.append(mineral)

    for mineral in negligible:
        y.pop(mineral)

    if variable == 'aH2O':
        x = data[variable].values
        plt.gca().invert_xaxis()
        plt.xlabel('aH2O')
    elif variable == 'fc':
        x = np.log10(data[variable].values)
        plt.xlabel('Log(fc)')
    elif variable == 'ionStr':
        x = np.log10(data[variable].values)
        plt.xlabel('Log(Ionic Strength)')
    elif variable == 'salinity':
        x = data[variable].values
        plt.xlabel('Salinity (g/kg)')
    else:
        raise AttributeError("'{variable}' does not exist".format(variable=variable))

    for mineral in y:
        if colormap and mineral in colormap:
            plt.plot(x, y[mineral], color=colormap[mineral], label=mineral)
        else:
            plt.plot(x, y[mineral], label=mineral)

    plt.yscale('log')
    plt.ylabel('mmol')

    plt.title(title)
    plt.legend()
    plt.savefig(outfile)
    plt.close()

def phase_diagram(input_file, output_dir, variable, title, plot_file, colormap):

    all_files = file_handler.file_detect(output_dir, '&')
    filetypes = []
    for i in all_files:
        filetypes.append('.' + i.split('.')[1])

    # find unique minerals
    all_minerals = []

    for file in all_files:
        data = pd.read_csv(r'{odir}/{f}'.format(odir=output_dir, f=file))

        minerals = pd.Series(data['label']).tolist()
        minerals = np.unique(minerals).tolist()

        if 'No_minerals' in minerals:
            minerals.remove('No_minerals')

        for seq in minerals:
            for mineral in seq.split('_'):
                all_minerals.append(mineral)

    unique_minerals = np.unique(all_minerals).tolist()

    inputs = pd.read_csv(input_file)
    trials = inputs.label
    temperature = inputs.t.values

    first_precip = np.ones((len(unique_minerals), temperature.shape[0]))
    last_precip = np.ones((len(unique_minerals), temperature.shape[0]))

    for key, t in enumerate(temperature):
        trial = trials[key]
        file = trial+filetypes[key]
        df = pd.read_csv(output_dir+'/'+file)

        sequence = df['label'].values

        for idx, mineral in enumerate(unique_minerals):
            indexer = np.ones(df.shape[0])
            for step in np.arange(0, df.shape[0]):
                if mineral in str(sequence[step]):
                    pass
                else:
                    indexer[step] = 0

            subset = df[indexer == 1]
            if subset.shape[0] == 0:
                first_precip[idx, key] = np.nan
                last_precip[idx, key] = np.nan
            else:
                first_precip[idx, key] = subset.iloc[0,:][variable]
                last_precip[idx, key] = subset.iloc[-1,:][variable]

    plt.figure(figsize=(11,8.5))

    if colormap:
        colormap = color_mapper.color_read(colormap)

        for idx, mineral in enumerate(unique_minerals):
            x = temperature
            y = first_precip[idx]
            if mineral in colormap:
                plt.plot(x, y, color=colormap[mineral], linestyle='-', label=mineral)
            else:
                plt.plot(x, y, '-', label=mineral)

    else:
        for idx, mineral in enumerate(unique_minerals):
            x = temperature
            y = first_precip[idx]
            plt.plot(x, y, '-', label=mineral)

    plt.legend()

    if variable == 'aH2O':
        plt.gca().invert_yaxis()
        plt.ylabel('Activity of H2O in Water (aH2O)')
    elif variable == 'fc':
        plt.yscale('log')
        plt.ylabel('Concentration Factor (fc)')
    elif variable == 'ionStr':
        plt.ylabel('Ionic Strength (M)')
    elif variable == 'salinity':
        plt.ylabel('Salinity (g/kg)')

    plt.xlabel('Temperature (°C)')
    plt.title(title)
    plt.savefig(plot_file, dpi=300)
    plt.close()
    
def compile_xyz(mineral, base, co, first_last, var='aH2O'):
    directories = os.listdir(base)
    Y = []
    for i in directories:
        Y.append(int(i.split('ppm')[0]))
    Y.sort()
    Y = np.array(Y)
    X = np.arange(0, 50, 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros(Y.shape)
    for i in range(0, Y.shape[0]):
        for k in range(0, Y.shape[1]):
            pco2 = Y[i,0]
            t = X[i, k] + 1
            
            f = r'{base}/{pco2}ppm_output/{co}/output/trial{t}.j{co0}&'.format(base=base, pco2=pco2, co=co, t=t, co0=co[0])
            df = pd.read_csv(f)
    
            target = np.nan
            
            indexer = np.ones(df.shape[0])
            for row in range(0, df.shape[0]):
                sequence = df['label'].iloc[row]
                indexer[row] = indexer[row] * (mineral in sequence)
            subset = df.loc[indexer == 1]
            subset.shape[0]
            if subset.shape[0] > 0:
                if first_last == 'first':
                    target = subset.loc[subset['fc'] == subset['fc'].min()][var].values[0]
                elif first_last == 'last':
                    target = subset.loc[subset['fc'] == subset['fc'].max()][var].values[0]
                else:
                    raise ValueError('invalid entry for argument \'first_last\'. enter \'first\' or \'last\'')
            Z[i, k] = target
            
    return(np.array([X, Y, Z]))
    
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