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
    
    pco2_ppm = 3000
    log_pco2 = np.log10(pco2_ppm / 1e6)
    temp = 30
    system = "c"
    
    plot = True
            

    test = eqlevp.simulation(label='test', temp=temp, dens=1, ph=6.55, na=84.5, 
                              k=3.3, ca=2.7, mg=1.3, cl=39.5, 
                              alk=56.2, pco2=log_pco2, system=system,
                              units='molar', add=['calcite'], 
                              remove=['dolomite','nesquehonite','brucite',
                                      'magnesite','hydromagnesite'],
                              verbose='0')
    test.run_simulation()
    

    t2 = perf_counter()
    
    print()
    print("Simulation time: {} seconds".format(round(t2-t1, 2)))

    if plot:
        min_file = "{}.j{}%".format(test.label, test.system)
        df = pd.read_csv(min_file)

        x = np.log10(df.fc.values)
        for i in range(1, df.shape[1]):
            y = df.iloc[:,i]
            label = df.columns[i]
        
            plt.plot(x, y, label=label)

        plt.title("{}ºC, {}ppm, {} system".format(temp, pco2_ppm, "closed" if test.system == "c" else "open"))
        plt.xlabel("log(fc)")
        plt.ylabel("moles precipitated")
        plt.yscale('log')
        plt.legend()
        plt.show()
