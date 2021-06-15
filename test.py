#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 08:58:47 2021

@author: warnuk
"""
import eqlevp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


if __name__ == "__main__":
    from time import perf_counter
    
    t1 = perf_counter()
    
    pco2_ppm = 1200
    log_pco2 = np.log10(pco2_ppm / 1e6)
    
    temp = 30
    
    test = eqlevp.simulation(label='test', temp=temp, dens=1, ph=6.55, na=84.5, 
                             k=3.3, li=0, ca=2.7, mg=1.3, cl=39.5, so4=0, 
                             alk=56.2, no3=0, si=0, b=0)
    
    test.run_eql(log_pco2, "o", units="molar", add_minerals=['calcite'],
                 rem_minerals=['dolomite', 'nesquehonite', 'brucite',
                               'magnesite', 'hydromagnesite', 'antarcticite',
                               'aragonite', 'burkeite', 'glaserite'],
                 verbose=True, call_evp=True)

    t2 = perf_counter()
    
    print()
    print("Simulation time: {} seconds".format(round(t2-t1, 2)))

    df = pd.read_csv(test.min_file)

    x = np.log10(df.fc.values)
    for i in range(1, df.shape[1]):
        y = df.iloc[:,i]
        label = df.columns[i]
        
        plt.plot(x, y, label=label)

    plt.title("30ºC, 1000ppm, open system")
    plt.xlabel("log(fc)")
    plt.ylabel("moles precipitated")
    plt.yscale('log')
    plt.legend()
    plt.show()