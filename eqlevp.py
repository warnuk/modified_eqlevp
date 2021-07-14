#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 20:55:23 2021

@author: warnuk
"""

import os
import math


class simulation:
    
    def __init__(self, label, temp=None, dens=None, ph=None, na=None, k=None, 
                 li=None, ca=None, mg=None, cl=None, so4=None, no3=None, 
                 b=None, si=None, alk=None, pco2=None, system='c', 
                 units='molar', dil=None, add=None, remove=None, max_sal=None, 
                 pkmol=None, pkeq=None, output_step=None, print_step=None, 
                 output=None, verbose=None, increment=None):
        self.label = label
        self.temp = temp
        self.dens = dens
        self.ph = ph
        self.na = na
        self.k = k
        self.li = li
        self.ca = ca
        self.mg = mg
        self.cl = cl
        self.so4 = so4
        self.no3 = no3
        self.b = b
        self.si = si
        self.alk = alk
        self.pco2 = pco2
        self.system = system
        self.units = units
        self.dil = dil
        self.add = add
        self.remove = remove
        self.max_sal = max_sal
        self.pkmol = pkmol
        self.pkeq = pkeq
        self.output_step = output_step
        self.print_step = print_step
        self.output = output
        self.verbose = verbose
        self.increment = increment
        
        if self.add:
            if type(self.add) == list:
                self.add = ",".join(self.add)
        if self.remove:
            if type(self.remove) == list:
                self.remove = ",".join(self.remove)
    
    def run_simulation(self):
        with open("input.dat", "w+") as file:
            file.write(",".join(["label", self.label]))
            file.write("\n")
            
            if self.temp:
                file.write(",".join(["temp", str(self.temp)]))
                file.write("\n")
            if self.dens:
                file.write(",".join(["dens", str(self.dens)]))
                file.write("\n")
            if self.ph:
                file.write(",".join(["ph", str(self.ph)]))
                file.write("\n")
            if self.na:
                file.write(",".join(["na", str(self.na)]))
                file.write("\n")
            if self.k:
                file.write(",".join(["k", str(self.k)]))
                file.write("\n")
            if self.li:
                file.write(",".join(["li", str(self.li())]))
                file.write("\n")
            if self.ca:
                file.write(",".join(["ca", str(self.ca)]))
                file.write("\n")
            if self.mg:
                file.write(",".join(["mg", str(self.mg)]))
                file.write("\n")
            if self.cl:
                file.write(",".join(["cl", str(self.cl)]))
                file.write("\n")
            if self.so4:
                file.write(",".join(["so4", str(self.so4)]))
                file.write("\n")
            if self.no3:
                file.write(",".join(["no3", str(self.no3)]))
                file.write("\n")
            if self.b:
                file.write(",".join(["b", str(self.b)]))
                file.write("\n")
            if self.si:
                file.write(",".join(["si", str(self.si)]))
                file.write("\n")
            if self.alk:
                file.write(",".join(["alk", str(self.alk)]))
                file.write("\n")
            if self.pco2:
                file.write(",".join(["pco2", str(self.pco2)]))
                file.write("\n")
            if self.system:
                file.write(",".join(["system", self.system]))
                file.write("\n")
            if self.units:
                file.write(",".join(["units", self.units]))
                file.write("\n")
            if self.dil:
                file.write(",".join(["dil", str(self.dil)]))
                file.write("\n")
            if self.add:
                file.write(",".join(["add", self.add]))
                file.write("\n")
            if self.remove:
                file.write(",".join(["remove", self.remove]))
                file.write("\n")
            if self.max_sal:
                file.write(",".join(["max_sal", str(self.max_sal)]))
                file.write("\n")
            if self.pkmol:
                file.write(",".join(["pkmol", str(self.pkmol)]))
                file.write("\n")
            if self.pkeq:
                file.write(",".join(["pkeq", str(self.pkeq)]))
                file.write("\n")
            if self.output_step:
                file.write(",".join(["output_step", str(self.output_step)]))
                file.write("\n")
            if self.print_step:
                file.write(",".join(["print_step", str(self.print_step)]))
                file.write("\n")
            if self.output:
                file.write(",".join(["output", str(self.output)]))
                file.write("\n")
            if self.verbose:
                file.write(",".join(["verbose", str(self.verbose)]))
                file.write("\n")
            if self.increment:
                file.write(",".join(["increment", str(self.increment)]))
                file.write("\n")

            file.close()
            
            os.system("./eql")