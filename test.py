#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 08:58:47 2021

@author: warnuk
"""
import eqlevp
from math import log, exp

if __name__ == "__main__":
    from time import perf_counter
    
    t1 = perf_counter()
    
    pco2_ppm = 1000
    log_pco2 = log(pco2_ppm / 1e6, 10)
    
    temp = 30
    
    test = eqlevp.simulation(label='test', temp=temp, dens=1, ph=6.55, na=84.5, 
                             k=3.3, li=0, ca=2.7, mg=1.3, cl=39.5, so4=0, 
                             alk=56.2, no3=0, si=0, b=0)
    
    test.run_eql(log_pco2, "c", add_minerals=['calcite'],
                 rem_minerals=['dolomite', 'nesquehonite', 'brucite',
                               'magnesite', 'hydromagnesite'],
                 verbose=True, call_evp=True)

    t2 = perf_counter()
    
    print(t2-t1)