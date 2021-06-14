#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 08:58:47 2021

@author: warnuk
"""
import eqlevp
import os
import shutil
import numpy as np
from time import perf_counter

if __name__ == "__main__":
    for pco2 in range(400, 2001, 100):
        for option in ['c', 'o']:
            base = r'{}ppm/{}/'.format(pco2, option)
            os.makedirs(base)
            
            for temperature in range(10, 41):
                t1 = perf_counter()
                
                log_pco2 = np.log10(pco2/1e6)
            
                sim = eqlevp.simulation(label='{}degrees'.format(temperature),
                                        temp=temperature, dens=1, ph=6.55, na=84.5,
                                        k=3.3, li=0, ca=2.7, mg=1.3, cl=39.5,
                                        so4=0, alk=56.2, no3=0, si=0, b=0)
                                        
                sim.run_eql(log_pco2, option, add_minerals=['calcite'],
                            rem_minerals=['dolomite', 'nesquehonite', 'brucite',
                                            'magnesite', 'hydromagnesite'],
                            verbose=False, call_evp=True)
                
                for file in [sim.log_file, sim.chem_file,
                                sim.event_file, sim.min_file,
                                sim.transfer_file]:
                    shutil.move(file, base+file)
                
                
                t2 = perf_counter()
            
                print("{} ppm ... {} ... "\
                      "{}ºC ... {} s".format(pco2, option, 
                                             temperature, round(t2-t1, 3)))
