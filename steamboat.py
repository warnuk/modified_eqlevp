#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 20:55:23 2021

@author: warnuk
"""

import os
import math
import time
import shutil
import eqlevp
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# set a base water name for the simulation labels and folders
water_name = "steamboat"

# enter water chemistry
ph = 6.55
na = 84.45
k = 3.32
li = None
ca = 2.74
mg = 1.28
cl = 39.49
so4 = 0
no3 = None
b = None
si = 0.
alk = 55.03

# add minerals to database
add = ['calcite']

# remove minerals from database
remove = ['nesquehonite','brucite','magnesite','hydromagnesite','dolomite']

# select open ("o") and/or closed ("c")
systems = ["c"]

# Set the range of temperatures (ºC) to run and the increment
low_T = 1
high_T = 50
inc_T = 1

# Set the range of CO2 concentrations (ppm) to run and the increment
low_pCO2 = 400
high_pCO2 = 2000
inc_pCO2 = 100

# show text output? True or False
verbose = False

# make plots for each simulation? True or False
plot = True


base = os.getcwd()
sim_dir = "{}_{}".format(water_name, datetime.date.today().strftime("%m%d%y"))


for system in systems:
    long_system = "closed" if system == "c" else "open"

    for pco2 in np.arange(low_pCO2, high_pCO2+inc_pCO2, inc_pCO2):
        long_pco2 = "{}ppm".format(pco2)

    
        filepath = os.path.join(base, sim_dir, long_system, long_pco2, "output")
        if not os.path.exists(filepath) and not os.path.isdir(filepath):
            os.makedirs(filepath)
        
        if plot:
            plotpath = os.path.join(base, sim_dir, long_system, long_pco2, "plots")
            if not os.path.exists(plotpath) and not os.path.isdir(plotpath):
                os.makedirs(plotpath)

        for temp in np.arange(low_T, high_T+inc_T, inc_T):

            t1 = time.perf_counter()

            label = "{}_{}_{}ppm_{}C".format(water_name, system, pco2, temp)
            log_pco2 = math.log10(pco2/1e6)

            sim = eqlevp.simulation(label=label, temp=temp, dens=1, ph=ph, na=na, 
                              k=k, li=li, ca=ca, mg=mg, cl=cl, so4=so4, no3=no3, 
                              b=b, si=si, alk=alk, pco2=log_pco2, system=system,
                              units='molal', add=add, remove=remove,
                              verbose=0 if not verbose else 1)
            sim.run_simulation()

            
            min_file = "{}.j{}%".format(label, system)
            chem_file = "{}.j{}&".format(label, system)
            event_file = "{}.j{}@".format(label, system)
            log_file = "{}.log".format(label)
            tra_file = "{}.tra".format(label)

            files = [min_file, chem_file, event_file, log_file]
            for file in files:
                shutil.move(file, os.path.join(filepath, file))
            os.remove(tra_file)

            min_file = os.path.join(filepath, min_file)
            chem_file = os.path.join(filepath, chem_file)

            

            t2 = time.perf_counter()
            dt = t2-t1

            print('%40s' % label, '%15s' % '{} seconds'.format(round(dt, 2)))

            if plot:

                fig, (ax1, ax2) = plt.subplots(2, figsize=(8.5,11))

                chemdata = pd.read_csv(chem_file)
                x = np.log10(chemdata['fc'].values)
                ions = ['na', 'k', 'alk', 'ca', 'mg', 'cl']
                for ion in ions:
                    y = chemdata[ion].values
                    ax1.plot(x, y, label=ion)

                df = pd.read_csv(min_file)

                x = np.log10(df.fc.values)
                for i in range(1, df.shape[1]):
                    y = df.iloc[:,i]
                    l = df.columns[i]
                
                    ax2.plot(x, y, label=l)
                
                ax2.set_xlabel("log(fc)")
                ax1.sharex(ax2)

                ax1.set_ylabel("mmol/kg H2O")
                ax2.set_ylabel("mmol precipitated")

                ax1.set_yscale("log")
                ax2.set_yscale("log")

                ax1.legend()
                ax2.legend()

                fig.suptitle("{}, {}ºC, {}, {} system".format(water_name.capitalize(), temp, long_pco2, long_system))
                plt.savefig(os.path.join(plotpath, "{}.png".format(label)), dpi=300)
                plt.close()
