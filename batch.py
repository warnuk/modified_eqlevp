#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 20:55:23 2021

@author: warnuk
"""

import os
import math

for pco2 in range(400, 500, 50):
    for s in ['c', 'o']:
        for t in range(1, 11):
            label = "test_{}ppm_{}_{}C".format(pco2, s, t)
            with open("input.dat", "w+") as file:
                file.write(",".join(["label", label]))
                file.write("\n")
                file.write(",".join(["temp", str(t)]))
                file.write("\n")
                file.write(",".join(["dens", str(1)]))
                file.write("\n")
                file.write(",".join(["ph", str(6.55)]))
                file.write("\n")
                file.write(",".join(["na", str(84.5)]))
                file.write("\n")
                file.write(",".join(["k", str(3.3)]))
                file.write("\n")
                file.write(",".join(["ca", str(2.7)]))
                file.write("\n")
                file.write(",".join(["mg", str(1.3)]))
                file.write("\n")
                file.write(",".join(["cl", str(39.5)]))
                file.write("\n")
                file.write(",".join(["alk", str(56.2)]))
                file.write("\n")
                file.write(",".join(["pco2", str(math.log10(pco2/1e6))]))
                file.write("\n")
                file.write(",".join(["system", s]))
                file.write("\n")
                file.write(",".join(["units", "molar"]))
                file.write("\n")
                file.write(",".join(["remove", "brucite,dolomite,nesquehonite,magnesite,hydromagnesite"]))
                file.write("\n")
                file.write(",".join(["output", "0"]))
                file.write("\n")
                file.close()


            os.system("./eql")
            
            