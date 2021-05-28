# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:12:10 2021

@author: warnu
"""
"""
EQL
System: Na-K-Li-Ca-Mg-Cl-SO4-NO3-CO2-HCO3-CO3-H4SiO4-B(OH)3-H-OH-H2O
Initial equilibrium for the simulation of evaporation EVP
Ionic interaction (Pitzer)
Solution of the system of equations by Newton-Raphson method
Temperature dependence included for minerals
  and for some Pitzer's coefficients (0-50C)
Options to change PCO2 and to dilute the initial solution
Connected files: MATRICE1-MURTF2-MURTF3-COMPLEX3-COEFFT4-DENSITE
Program chained to EVP through transfer file: *.tra
Last revision: march 1999
author: F. RISACHER
conversion to FORTRAN 90: A. CLEMENT
conversion to Python 3: W. ARNUK
"""
import pandas as pd
import numpy as np
from datetime import datetime

def read_file(file):
    with open(file, 'r+') as file:
        text = list(filter(None, file.read().split('\n')))
        file.close()
        values = [line.split(',') for line in text]
    return(values)

#complex3 = read_file("complex3")



# EQL local variables
repinit_S, fich_S, molemin_S, fichmin_S = ".", ".", ".", "."
nchcat, nchani = 0, 0
n = 25
ntot = 12
n0 = 10
max_cat = 5
max_an = 5
epsilon=1e-8





aqu = pd.read_csv("aqu.dat", header=None)
aq_S = aqu.iloc[:,0].values.astype(str)
atom = aqu.iloc[:,1].values
nch = aqu.iloc[:,2].values

text = read_file("matrice1")
kmat = np.zeros((n+1, n+1))
for i in range(1, n+1):
    for j in range(1, n+1):
        kmat[i, j] = float(text[i-1][j-1])     
kmat[13, 0], kmat[15, 0], kmat[20, 0], kmat[21, 0] = -1, -1, 3, 5

text = read_file("complex3")
complex3 = pd.DataFrame([i[1:] for i in text], dtype="float64", 
                        index=[i[0] for i in text],
                        columns = ['at', 'bt', 'ct', 'dt', 'et'])


text = read_file("murtf2")

(nc, na, nm) = (int(i) for i in text[0])
nt = nc + na

wmin = np.zeros((nm+1, nt+1))
lmin = np.zeros(nm+1)
nwmin = np.zeros(nm+1)

ion_db = pd.DataFrame(text[1:nt+2], columns = ['species', 'at', 'bt', 'ct', 'dt', 'et'])
ion_db['species'] = ion_db['species'].str.lower()
ion_db = ion_db.astype({'species': str,
                        'at': np.dtype("float"),
                        'bt': np.dtype("float"),
                        'ct': np.dtype("float"),
                        'dt': np.dtype("float"),
                        'et': np.dtype("float")})
ion_db = pd.merge(pd.Series(aq_S, name="species"), ion_db, 
                  how="left", on="species").dropna()

mineral_list = text[nt+2:]
mineral_S = [i[0] for i in mineral_list]
mineral_S.insert(0, None)
min_db = np.zeros((nm+1, 5))

for k in range(1, nm+1):
    line = text[nt+1+k]
    mineral_S[k] = line[0]
    ncomp = int(line[1])
    c_ion = np.zeros(ncomp+1)
    nom_ion_S = np.empty(ncomp+1, dtype=np.object)
    for i in range(1, ncomp+1):
        c_ion[i] = float(line[i*2])
        nom_ion_S[i] = str(line[1+i*2])
        for j in range(0, nt+1):
            x_S = nom_ion_S[i].lower()
            if x_S == aq_S[j]:
                wmin[k, j] = c_ion[i]
    (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
    min_db[k] = np.array([at, bt, ct, dt, et])

min_db = pd.DataFrame(min_db, index = mineral_S, columns=['at', 'bt', 'ct', 'dt', 'et'])
#mineral_S = np.array(mineral_S)


# LINE 260
min_S = "murtf3"
text = read_file(min_S)
nc, na, nm0 = [int(i) for i in text[0]]

#ion0 = text[1:(2+nc+na)]
mineral0_S = [i[0] for i in text[(2+nc+na):]]
nwmin[np.in1d(mineral_S, mineral0_S)] = 1

# Functions
def input_data(filepath):
    
    data = pd.read_csv(filepath, header=None, index_col=0)
    
    label = data.loc['label'][1]
    temp = float(data.loc['t'][1])
    dens = float(data.loc['dens'][1])
    ph = float(data.loc['ph'][1])
    na = float(data.loc['na'][1])
    k = float(data.loc['k'][1])
    li = float(data.loc['li'][1])
    ca = float(data.loc['ca'][1])
    mg = float(data.loc['mg'][1])
    cl = float(data.loc['cl'][1])
    so4 = float(data.loc['so4'][1])
    no3 = float(data.loc['no3'][1])
    b = float(data.loc['b'][1])
    si = float(data.loc['si'][1])
    alk = float(data.loc['alk'][1])
    
    return(Water(label, temp, dens, ph, na, k, li, ca, 
                 mg, cl, so4, no3, b, si, alk))
    

    

class Water:
    'Class for input waters to be used in EQL/EVP simulations'
    
    def __init__(self, label, temp, dens, ph, na, k, li, 
                 ca, mg, cl, so4, no3, b, si, alk, pc=None, 
                 unit_S="molal"):
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
        self.pc = pc
        self.unit_S = unit_S
        
        self.tot = np.zeros(ntot+1)
        self.tot[1] = self.na / 1000
        self.tot[2] = self.k / 1000
        self.tot[3] = self.li / 1000
        self.tot[4] = self.ca / 1000
        self.tot[5] = self.mg / 1000
        self.tot[6] = self.cl / 1000
        self.tot[7] = self.so4 / 1000
        self.tot[8] = self.no3 / 1000
        self.tot[9] = self.b / 1000
        self.tot[10] = self.si / 1000
        self.tot[11] = 10 ** -self.ph
        self.tot[12] = self.alk / 1000
        
        self.tinit = self.temp
        self.ica = np.ones(n+1)
        self.ica[0] = 0
        self.ndepact = 0
        self.nminer = np.count_nonzero(self.tot[1:11]) + 1
        
        if self.dens == 0:
            self.dens = 1
        if self.dens == 1:
            self.unit_S = "molal"
        else:
            self.unit_S = "molar"
        
        self.tot0 = self.tot
        
        self.psc = np.zeros(15)
        self.psc[1:] = 10 ** (complex3['at'] + 
                              complex3['bt'] / 300 * self.temp + 
                              complex3['ct'] / 30000 * self.temp ** 2 + 
                              complex3['dt'] / 3000000 * self.temp ** 3 + 
                              complex3['et'] / 300000000 * self.temp ** 4).values
        self.mu = np.zeros(15)
        self.mu = (ion_db['at'] + 
                   ion_db['bt'] / 300 * self.temp + 
                   ion_db['ct'] / 30000 * self.temp ** 2 + 
                   ion_db['dt'] / 3000000 * self.temp ** 3 + 
                   ion_db['et'] / 300000000 * self.temp ** 4).values
        self.mum = np.zeros(nm+1)
        self.mum = (min_db['at'] + 
                    min_db['bt'] / 300 * self.temp +
                    min_db['ct'] / 30000 * temp ** 2 +
                    min_db['dt'] / 3000000 * temp ** 3 +
                    min_db['et'] / 300000000 * temp ** 4).values
        
        self.psol = np.zeros(nm+1)
        for k in range(1, nm+1):
            u = self.mum[k]
            for i in range(0, nt+1):
                u -= wmin[k, i] * self.mu[i]
                
            self.psol[k] = np.exp(u)

        #self.psol = np.exp(self.mum - np.sum((wmin[:,:] * self.mu), axis=1))
        self.gact = np.zeros(n+1)
        self.molal = np.zeros(n+1)
        self.totinit = np.zeros(n0+1)
        self.act = np.zeros(n+1)
        self.z = np.zeros((n+1, n+1))
        self.zz = np.zeros(n+1)
        self.xx = np.zeros(n+1)
        self.cat = np.zeros(max_cat+1)
        self.nchcat = np.zeros(max_cat+1)
        self.ani = np.zeros(max_an+1)
        self.nchani = np.zeros(max_an+1)
        self.nconv = 2
        self.pk = .1
        self.pk0 = .1
        self.pkf = .0001
        self.pkstd = .1
        self.dph = .2
        self.dil = 0.
        self.diltot = 1
        self.pco2_S = ""
        self.mh2o = 55.51
        self.eps = 1.0e-12
        self.poa = 0.
        self.std = 0.
        self.nch = nch
        
    def temperature(self, at, bt, ct, dt, et):
        return((at + bt * self.temp + ct * self.temp ** 2 + 
                dt * self.temp ** 3 + et * self.temp ** 4))
    
    def j0(self, x):
        ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
        j0 = x / (4. + ya * x ** yb * np.exp(yc * x ** yd))
        return(j0)
            
    def j1(self, x):
        ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
        j1 = ((4. + ya * x ** yb * (1. - yb - yc * yd * x ** yd) * 
          np.exp(yc * x ** yd)) / (4. + ya * x ** yb * 
                                   np.exp(yc * x ** yd)) ** 2)
        return(j1)
    
    def g0(self, x):
        g0 = 2. * (1. - (1. + x) * np.exp(-x)) / x ** 2
        return(g0)
            
    def g1(self, x):
        g1 = -2. * (1. - (1. + x + x ** 2 / 2.) * np.exp(-x)) / x ** 2
        return(g1)
    
    def charge_balance(self, verbose=True):
        sc, cmax = 0, 0
        for i in range(1, ntot+1):
            if self.nch[i] > 0 and i != 11:
                sc += self.tot[i] * self.nch[i]
                if self.tot[i] * self.nch[i] >= cmax:
                    cmax = self.tot[i] * self.nch[i]
                    icat = i
        sa, amax = 0, 0
        for i in range(1, ntot+1):
            if self.nch[i] < 0:
                sa += self.tot[i] * -self.nch[i]
                if self.tot[i] * -self.nch[i] >= amax:
                    amax = self.tot[i] * -self.nch[i]
                    iani = i
        
        if verbose:
            print("Sum of cations = {}".format(sc))
            print("Sum of anions = {}".format(sa))
        
        if sc + sa != 0:
            self.dca = 200 * np.abs(sc - sa) / (sc + sa)
            delta = sc - sa
            self.tot[icat] = self.tot[icat] - delta / 2 / nch[icat]
            self.tot[iani] = self.tot[iani] + delta / 2 / -(nch[iani])
            if verbose:
                print("electrical balance = {}%".format((self.dca * 100 + 
                                                         0.5) / 100))
                print()
              
    def calculate_molalities(self, verbose=True):
        self.molal[1] = self.tot[1]
        self.molal[2] = self.tot[2]
        self.molal[3] = self.tot[3]
        self.molal[6] = self.tot[6]
        self.molal[8] = self.tot[8]
        self.molal[11] = self.tot[11]
        self.molal[13] = self.psc[1] / self.molal[11]
        if self.tot[9] > 0:
            a = 1 + self.molal[13] / self.psc[7]
            b = 3 * self.psc[8] * self.molal[13]
            c = 4 * self.psc[9] * self.molal[13] ** 2
            xu = self.tot[9] / 2
            u = xu
            
            eq = a * xu + b * xu ** 3 + c * xu ** 4
            while (200 * abs(eq - self.tot[9]) / 
                           (eq + self.tot[9]) >= self.pk):
                
                u = u / 2
                if eq > self.tot[9]:
                    xu -= u
                else:
                    xu += u
                
                eq = a * xu + b * xu ** 3 + c * xu ** 4
            
            self.molal[9] = xu
            self.molal[19] = (self.molal[9] * self.molal[13] / self.psc[7])
            self.molal[20] = (self.molal[13] * self.molal[9] ** 3 * 
                              self.psc[8])
            self.molal[21] = (self.molal[13] ** 2 * self.molal[9] ** 4 * 
                              self.psc[9])
        self.molal[14] = (self.tot[12] + self.molal[11] - self.molal[13] - 
                          self.molal[19] - self.molal[20] - 
                          2 * self.molal[21]) / (2 + self.molal[11] / 
                                                 self.psc[2])
        self.molal[12] = (self.tot[12] + self.molal[11] - self.molal[13] - 
                          self.molal[19] - self.molal[20] - 
                          2 * self.molal[21]) / (1 + 2 * self.psc[2] / 
                                                 self.molal[11])
                                                 
        self.molal[15] = self.molal[12] * self.molal[11] / self.psc[3]
        self.molal[4] = self.tot[4] / (1 + self.molal[14] / self.psc[4] +
                                       self.molal[19] / self.psc[10])
        self.molal[16] = self.molal[4] * self.molal[14] / self.psc[4]
        self.molal[22] = self.molal[4] * self.molal[19] / self.psc[10]
        self.molal[5] = self.tot[5] / (1 + self.molal[14] / self.psc[5] + 
                                       self.molal[13] / self.psc[6] + 
                                       self.molal[19] / self.psc[11])
        self.molal[17] = self.molal[5] * self.molal[14] / self.psc[5]
        self.molal[18] = self.molal[5] * self.molal[13] / self.psc[6]
        self.molal[23] = self.molal[5] * self.molal[19] / self.psc[11]
        self.molal[10] = self.tot[10] / (1 + self.psc[12] / self.molal[11])
        self.molal[24] = self.tot[10] / (1 + self.molal[11] / self.psc[12])
        self.molal[7] = self.tot[7] / (1 + self.molal[11] / self.psc[13])
        self.molal[25] = self.molal[7] * self.molal[11] / self.psc[13]
        
        sc, cmax = 0, 0
        for i in range(1, n+1):
            if self.nch[i] > 0:
                sc += self.molal[i] * self.nch[i]
                if self.molal[i] * self.nch[i] > cmax:
                    cmax = self.molal[i] * self.nch[i]
                    self.icat = i
        
        sa, amax = 0, 0
        for i in range(1, n+1):
            if self.nch[i] < 0:
                sa += self.molal[i] * -self.nch[i]
                if self.molal[i] * -self.nch[i] > amax:
                    amax = self.molal[i] * -self.nch[i]
                    self.iani = i
        
        delta = sc - sa
        self.molal[self.icat] = self.molal[self.icat] - delta / 2 / nch[self.icat]
        self.molal[self.iani] = self.molal[self.iani] + delta / 2 / (-nch[self.iani])
        
        sc = (self.molal[1] + self.molal[2] + self.molal[3] + 
              self.molal[4] * 2 + self.molal[5] * 2 + self.molal[11] + 
              self.molal[18] + self.molal[22] + self.molal[23])
        sa = (self.molal[6] + self.molal[7] * 2 + self.molal[8] + 
              self.molal[12] + self.molal[13] + self.molal[14] * 2 +
              self.molal[19] + self.molal[20] + self.molal[21] * 2 +
              self.molal[24] + self.molal[25])
        
        if verbose:
            print("\nSum of cations = {} corrected for {}".format(sc, 
                                                                  aq_S[self.icat]))
            print("Sum of anions = {} corrected for {}\n".format(sa, 
                                                                 aq_S[self.iani]))
            
        s = 0
        for i in range(1, n+1):
            s += self.molal[i] * atom[i]
        
        if self.unit_S == "molar":
            self.ee = 1000 / (1000 * self.dens - s)
            
            for i in range(1, ntot+1):
                if i != 11:
                    self.tot[i] = self.tot[i] * self.ee
            for i in range(1, n+1):
                if i != 11:
                    self.molal[i] = self.molal[i] * self.ee

        elif self.unit_S == "molal":
            self.ee = 1
        
        self.stdi = s * self.ee
    
    def actp(self):
        c = np.zeros(10)
        a = np.zeros(12)
        h = np.zeros(4)
        
        c[1] = self.molal[1]
        c[2] = self.molal[2]
        c[3] = self.molal[3]
        c[4] = self.molal[4]
        c[5] = self.molal[22]
        c[6] = self.molal[5]
        c[7] = self.molal[18]
        c[8] = self.molal[23]
        c[9] = self.molal[11]
        a[1] = self.molal[6]
        a[2] = self.molal[7]
        a[3] = self.molal[25]
        a[4] = self.molal[12]
        a[5] = self.molal[14]
        a[6] = self.molal[13]
        a[7] = self.molal[24]
        a[8] = self.molal[19]
        a[9] = self.molal[20]
        a[10] = self.molal[21]
        a[11] = self.molal[8]
        h[1] = self.molal[10]
        h[2] = self.molal[9]
        h[3] = self.molal[15]

        if self.ndepact == 0:
            text = read_file("coefft4")
        
            (self.nc, self.na, self.nn) = (int(i) for i in text[0])
            self.nzc, self.nza = np.zeros(self.nc+1), np.zeros(self.na+1)
            (self.b0, self.b1, self.b2, self.c0) = (np.zeros((self.nc+1, self.na+1), dtype=np.float64), 
                                np.zeros((self.nc+1, self.na+1), dtype=np.float64), 
                                np.zeros((self.nc+1, self.na+1), dtype=np.float64), 
                                np.zeros((self.nc+1, self.na+1), dtype=np.float64))
            self.sc = np.zeros((self.nc+1, self.nc+1, self.na+1), dtype=np.float64)
            self.sa = np.zeros((self.na+1, self.na+1, self.nc+1), dtype=np.float64)
            self.tc = np.zeros((self.nc+1, self.nc+1), dtype=np.float64)
            self.ta = np.zeros((self.na+1, self.na+1), dtype=np.float64)
            self.lc = np.zeros((self.nn+1, self.nc+1), dtype=np.float64)
            self.la = np.zeros((self.nn+1, self.na+1), dtype=np.float64)
            self.xi = np.zeros((self.nn+1, self.nc+1, self.na+1), dtype=np.float64)
            
            self.nzc[1:] = [int(i[1]) for i in text[1:1+self.nc]]
            self.nza[1:] = [int(i[1]) for i in text[1+self.nc:1+self.nc+self.na]]
            
            (at, bt, ct, dt, et) = (np.float64(i.replace('d', 'e')) 
                                    for i in text[1+self.nc+self.na][1:])
            
            self.ap0 = self.temperature(at, bt, ct, dt, et)
            
            index = 1+self.nc+self.na
            
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.b0[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.b1[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.b2[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.c0[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for i in range(1, self.nc):
                for j in range(i+1, self.nc+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.tc[i, j] = self.temperature(at, bt, ct, dt, et)
                    self.tc[j, i] = self.tc[i, j]
                    
            for i in range(1, self.na):
                for j in range(i+1, self.na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.ta[i, j] = self.temperature(at, bt, ct, dt, et)
                    self.ta[j, i] = self.ta[i, j]
                    
            for k in range(1, self.nc):
                for i in range(k+1, self.nc+1):
                    for j in range(1, self.na+1):
                        index += 1
                        (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                                for k in text[index][1:])
                        self.sc[k, i, j] = self.temperature(at, bt, ct, dt, et)
                        self.sc[i, k, j] = self.sc[k, i, j]
            
            for k in range(1, self.na):
                for i in range(k+1, self.na+1):
                    for j in range(1, self.nc+1):
                        index += 1
                        (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                                for k in text[index][1:])
                        self.sa[k, i, j] = self.temperature(at, bt, ct, dt, et)
                        self.sa[i, k, j] = self.sa[k, i, j]
                        
            for i in range(1, self.nn+1):
                for j in range(1, self.nc+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.lc[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for i in range(1, self.nn+1):
                for j in range(1, self.na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    self.la[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for k in range(1, self.nn+1):
                for i in range(1, self.nc+1):
                    for j in range(1, self.na+1):
                        self.xi[k, i, j] = 0
                        
            self.xi[2, 9, 1] = -0.0102
            self.xi[2, 1, 2] = 0.046
        
        ec, ea = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                  np.zeros((self.na+1, self.na+1), dtype=np.float64))
        fc, fa = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                  np.zeros((self.na+1, self.na+1), dtype=np.float64))
        xc, xa = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                  np.zeros((self.na+1, self.na+1), dtype=np.float64))
        pp, qp = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                  np.zeros((self.na+1, self.na+1), dtype=np.float64))
        p, q = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                np.zeros((self.na+1, self.na+1), dtype=np.float64))
        pf, qf = (np.zeros((self.nc+1, self.nc+1), dtype=np.float64), 
                  np.zeros((self.na+1, self.na+1), dtype=np.float64))
        cc, bf = (np.zeros((self.nc+1, self.na+1), dtype=np.float64), 
                  np.zeros((self.nc+1, self.na+1), dtype=np.float64))
        b, bp = (np.zeros((self.nc+1, self.na+1), dtype=np.float64), 
                 np.zeros((self.nc+1, self.na+1), dtype=np.float64))
        gc, ga, gn = (np.zeros(self.nc+1, dtype=np.float64), 
                      np.zeros(self.na+1, dtype=np.float64), 
                      np.zeros(self.nn+1, dtype=np.float64))
        
        self.bp0 = 1.2e0
        
        u, z = 0, 0
        
        u += np.sum(c * self.nzc ** 2)
        z += np.sum(c * self.nzc)
        u += np.sum(a * self.nza ** 2)
        z += np.sum(a * self.nza)

        self.fi = u / 2
        fj = np.sqrt(self.fi)
        u = 6 * self.ap0 * fj        
        
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                
                if self.nzc[i] == self.nzc[j]:
                    ec[i, j] = 0
                    fc[i, j] = 0
                    
                else:
                    xc[i, j] = 2 * u
                    xc[i, i] = self.nzc[i] ** 2 * u
                    xc[j, j] = self.nzc[j] ** 2 * u
                    ec[i, j] = ((self.j0(xc[i, j]) - self.j0(xc[i, i]) / 2 - 
                                self.j0(xc[j, j]) / 2) / self.fi / 2)
                    fc[i, j] = ((xc[i, j] * self.j1(xc[i, j]) - xc[i, i] * 
                                self.j1(xc[i, i]) / 2 - xc[j, j] * 
                                self.j1(xc[j, j]) / 2) / self.fi ** 2 / 4 - 
                                ec[i, j] / self.fi)
                    ec[j, i] = ec[i, j]
                    fc[j, i] = fc[i, j]
                    
        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                if self.nza[i] == self.nza[j]:
                    ea[i, j] = 0
                    fa[i, j] = 0
                else:
                    xa[i, j] = 2 * u
                    xa[i, i] = self.nza[i] ** 2 * u
                    xa[j, j] = self.nza[j] ** 2 * u
                    ea[i, j] = (self.j0(xa[i, j]) - self.j0(xa[i, i]) / 2 -
                                self.j0(xa[j, j]) / 2) / self.fi / 2
                    fa[i, j] = ((xa[i, j] * self.j1(xa[i, j]) - xa[i, i] * 
                                self.j1(xa[i, i]) / 2 - xa[j,j] * 
                                self.j1(xa[j, j]) / 2) / self.fi ** 2 / 4 - 
                                ea[i, j] / self.fi)
                    ea[j, i] = ea[i, j]
                    fa[j, i] = fa[i, j]
        
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                pp[i, j] = fc[i, j]
                p[i, j] = self.tc[i, j] + ec[i, j]
                pf[i, j] = p[i, j] + pp[i, j] * self.fi
                pp[j, i] = pp[i, j]
                p[j, i] = p[i, j]
                pf[j, i] = pf[i, j]
        
        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                qp[i, j] = fa[i, j]
                q[i, j] = self.ta[i, j] + ea[i, j]
                qf[i, j] = q[i, j] + qp[i, j] * self.fi
                qp[j, i] = qp[i, j]
                q[j, i] = q[i, j]
                qf[j, i] = qf[i, j]
                
        w = fj * 12
        for i in range(1, self.nc+1):
            for j in range(1, self.na+1):
                cc[i, j] = self.c0[i, j] / np.sqrt(self.nzc[i] * self.nza[j]) / 2
                if self.nzc[i] == 2 and self.nza[j] == 2:
                    v = fj * 1.4e0
                if self.nzc[i] == 1 or self.nza[j] == 1:
                    v = fj * 2
                bf[i, j] = (self.b0[i, j] + self.b1[i, j] * np.exp(-v) + 
                            self.b2[i, j] * np.exp(-w))
                b[i, j] = (self.b0[i, j] + self.b1[i, j] * 
                           (self.g0(v)) + self.b2[i, j] * (self.g0(w)))
                bp[i, j] = (self.b1[i, j] * (self.g1(v)) / 
                            self.fi + self.b2[i, j] * (self.g1(w)) / self.fi)
        
        f = -self.ap0 * (fj / (1 + self.bp0 * fj) + 2 / self.bp0 * np.log(1 + self.bp0 * fj))
        
        i = np.repeat(np.arange(1, self.nc+1), self.na)
        j = np.tile(np.arange(1, self.na+1), self.nc)
        f += np.sum(c[i] * a[j] * bp[i, j])
        
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                f += c[i] * c[j] * pp[i, j]

        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                f += a[i] * a[j] * qp[i, j]
                
        for ii in range(1, self.nc+1):
            u = self.nzc[ii] ** 2 * f
            for j in range(1, self.na+1):
                u += a[j] * (b[ii, j] * 2 + z * cc[ii, j])
            for i in range(1, self.nc+1):
                if i != ii:
                    v = 0
                    for j in range(1, self.na+1):
                        v += a[j] * self.sc[ii, i, j]
                    u += c[i] * (p[ii, i] * 2 + v)
            for i in range(1, self.na):
                for j in range(i+1, self.na+1):
                    u += a[i] * a[j] * self.sa[i, j, ii]
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += c[i] * a[j] * cc[i, j] * self.nzc[ii]
            for i in range(1, self.nn+1):
                u += h[i] * self.lc[i, ii] * 2
            for k in range(1, self.nn+1):
                for j in range(1, self.na+1):
                    u += h[k] * a[j] * self.xi[k, ii, j]
            gc[ii] = np.exp(u)
            
        for jj in range(1, self.na+1):
            u = self.nza[jj] ** 2 * f
            for i in range(1, self.nc+1):
                u += c[i] * ((b[i, jj]) * 2 + z * cc[i, jj])
            for i in range(1, self.na+1):
                if i != jj:
                    v = 0
                    for j in range(1, self.nc+1):
                        v += c[j] * self.sa[jj, i, j]
                    u += a[i] * (q[jj, i] * 2 + v)
                        
            for i in range(1, self.nc):
                for j in range(i+1, self.nc+1):
                    u += c[i] * c[j] * self.sc[i, j, jj]
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += c[i] * a[j] * cc[i, j] * self.nza[jj]
            for j in range(1, self.nn+1):
                u += h[j] * self.la[j, jj]
            for k in range(1, self.nn+1):
                for i in range(1, self.nc+1):
                    u += h[k] * c[i] * self.xi[k, i, jj]
            ga[jj] = np.exp(u)
            
        for k in range(1, self.nn+1):
            u = 0
            for i in range(1, self.nc+1):
                u += c[i] * self.lc[k, i] * 2
            for j in range(1, self.na+1):
                u += a[j] * self.la[k, j] * 2
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += c[i] * a[j] * self.xi[k, i, j]
            gn[k] = np.exp(u)
            
        u = -self.ap0 * self.fi ** 1.5e0 / (1 + self.bp0 * fj)
        for i in range(1, self.nc+1):
            for j in range(1, self.na+1):
                u += c[i] * a[j] * (bf[i, j] + z * cc[i, j])
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                v = 0
                for k in range(1, self.na+1):
                    v += a[k] * self.sc[i, j, k]
                u += c[i] * c[j] * (pf[i, j] + v)
        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                v = 0
                for k in range(1, self.nc+1):
                    v += c[k] * self.sa[i, j, k]
                u += a[i] * a[j] * (qf[i, j] + v)
        for k in range(1, self.nn+1):
            for i in range(1, self.nc+1):
                u += h[k] * c[i] * self.lc[k, i]
        for k in range(1, self.nn+1):
            for j in range(1, self.na+1):
                u += h[k] * a[j] * self.la[k, j]
        for k in range(1, self.nn+1):
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += h[k] * c[i] * a[j] * self.xi[k, i, j]
        
        s = 0
        for i in range(1, self.nc+1):
            s += c[i]
        for j in range(1, self.na+1):
            s += a[j]
        co = 1 + 2 * u / s
        self.aw = np.exp(-s * co / self.mh2o)
        self.gact[1] = gc[1]
        self.gact[2] = gc[2]
        self.gact[3] = gc[3]
        self.gact[4] = gc[4]
        self.gact[22] = gc[5]
        self.gact[5] = gc[6]
        self.gact[18] = gc[7]
        self.gact[23] = gc[8]
        self.gact[11] = gc[9]
        self.gact[6] = ga[1]
        self.gact[7] = ga[2]
        self.gact[25] = ga[3]
        self.gact[12] = ga[4]
        self.gact[14] = ga[5]
        self.gact[13] = ga[6]
        self.gact[24] = ga[7]
        self.gact[19] = ga[8]
        self.gact[20] = ga[9]
        self.gact[21] = ga[10]
        self.gact[8] = ga[11]
        self.gact[10] = self.aw * self.aw * gn[1] ** np.log(10)
        self.gact[9] = gn[2]
        self.gact[15] = gn[3]
        self.gact[16] = 1
        self.gact[17] = 1
        self.ndepact = 1
            
    def density(self):
        nc, na = 5, 5
        
        s = np.zeros((nc+1, na+1))
        
        with open("densite", 'r+') as file:
            text = [line.split(',') for line in file.read().split('\n')[:-1]]
            file.close()
        
        
        (ao, bo, au, bu) = (np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)), 
                            np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)))
        index = 0
        for i in range(1, 6):
            for j in range(1, 6):
                (ao[i, j], bo[i, j], au[i, j], bu[i, j]) = (float(i) for i 
                                                            in text[index][1:])
                index += 1
        
        if self.unit_S == "molar":
            self.dens = 1
            u = np.sum(self.ani[1:na+1])
                
            i = np.repeat(np.arange(1, nc+1), na)
            j = np.tile(np.arange(1, na+1), nc)
            s[i, j] = ((self.nchcat[i] + self.nchani[j]) / 2 * self.cat[i] * 
                       self.ani[j] / self.nchcat[i] / self.nchani[j] / u)
            self.dens += np.sum(ao[i, j] * s[i, j] + bo[i, j] * s[i, j] ** 2)
            
        elif self.unit_S == "molal":
            self.dens = 1
            u = np.sum(self.ani[1:na+1] * self.nchani[1:na+1])
            i = np.repeat(np.arange(1, nc+1), na)
            j = np.tile(np.arange(1, na+1), nc)
            s[i, j] = ((self.nchcat[i] + self.nchani[j]) / 2 * self.cat[i] * 
                       self.ani[j] / u)
            self.dens += np.sum(au[i, j] * s[i, j] + bu[i, j] * s[i, j] ** 2)   

    def iterate_molalities(self, verbose=True):
        iterate = True
        while iterate:
            self.nu = 1
            self.ncompt = 0
            while self.nu != 0:
                self.actp()
                
                for i in range(1, n+1):
                    self.act[i] = self.molal[i] * self.gact[i]
                    
                self.act[0] = self.aw
                self.tot[11] = (10 ** (-self.ph)) / self.gact[11]
                self.act[11] = 10 ** (-self.ph)
                
                for i in range(1, 13):
                    for j in range(1, n+1):
                        if self.molal[j] != 0:
                            self.z[i, j] = kmat[i, j]
                    u = 0
                    for j in range(1, n+1):
                        u += kmat[i, j] * self.molal[j]
                    self.zz[i] = self.tot[i] - u
                
                for i in range(13, n+1):
                    for j in range(1, n+1):
                        if self.molal[j] != 0:
                            self.z[i, j] = kmat[i, j] / self.molal[j]
                        elif self.molal[j] == 0:
                            self.z[i, j] = 0
                    u = 0
                    for j in range(0, n+1):
                        if self.act[j] > 0:
                            u += kmat[i, j] * np.log(self.act[j])
                    self.zz[i] = np.log(self.psc[i-12]) - u
                
                for k in range(1, n0+1):
                    if self.tot[k] == 0 and k != 12:
                        self.ica[k] = 0
                        for i in range(k+1, n+1):
                            if kmat[i, k] != 0:
                                self.ica[i] = 0
                        for j in range(k+1, n+1):
                            if kmat[k, j] != 0:
                                self.ica[j] = 0
                ni, nj = n, n
                for k in range(n, 0, -1):
                    if self.ica[k] == 0:
                        for i in range(k, ni):
                            for j in range(1, nj+1):
                                self.z[i, j] = self.z[i+1, j]
                            self.zz[i] = self.zz[i+1]    
                        ni -= 1
                        for j in range(k, nj):
                            for i in range(1, ni+1):
                                self.z[i, j] = self.z[i, j+1]
                        nj -= 1
                for k in range(2, ni+1):
                    for i in range(k, ni+1):
                        if self.z[i, k-1] != 0:
                            u = self.z[k-1, k-1] / self.z[i, k-1]
                            for j in range(k, ni+1):
                                self.z[i, j] = self.z[k-1, j] - self.z[i, j] * u
                            self.zz[i] = self.zz[k-1] - self.zz[i] * u
                self.xx[ni] = self.zz[ni] / self.z[ni, ni]
                
                for i in range(ni-1, 0, -1):
                    s = 0
                    for j in range(i+1, ni+1):
                        s += self.z[i, j] * self.xx[j]
                    self.xx[i] = (self.zz[i] - s) / self.z[i, i]
                    
                for k in range(1, n+1):
                    if self.ica[k] == 0:
                        for i in range(ni, k-1, -1):
                            self.xx[i+1] = self.xx[i]
                        self.xx[k] = 0
                        ni += 1
                
                self.ncompt += 1
                if verbose:
                    print("iteration molalities {}".format(self.ncompt))
                
                if self.ncompt >= 100:
                    for i in range(1, n+1):
                        if self.molal[i] + self.xx[i] / self.nconv < 0:
                            print("the equation set diverges: end of program")
                            return()
                
                for i in range(1, n+1):
                    if self.molal[i] + self.xx[i] / self.nconv < 0:
                        self.molal[i] = self.eps
                    else:
                        self.molal[i] += self.xx[i] / self.nconv
                
                self.nu = 0
                for i in range(1, n+1):
                    if self.ica[i] == 1:
                        if (200 * np.abs(self.xx[i] / self.nconv / 
                                         (2 * self.molal[i] - 
                                          self.xx[i] / self.nconv)) > self.pk):
                            self.nu = 1
            
            self.std = 0
            for i in range(0, n+1):
                self.std += self.molal[i] * atom[i]
            
            if verbose:
                print("tdsi = {}".format(self.stdi))
                print("tds = {}".format(self.std))
            
            if (np.abs(self.std-self.stdi)/
               (self.std+self.stdi)*200 < self.pkstd):
                iterate = False
                break
            
            else:
                if self.unit_S == "molar":
                    self.ef = (1000 + self.std) / self.dens / 1000
                    for i in range(1, ntot+1):
                        if i != 11:
                            self.tot[i] = self.tot[i] / self.ee * self.ef
                    for i in range(0, n+1):
                        if i != 11:
                            self.molal[i] = self.molal[i] / self.ee * self.ef
                    self.ee = self.ef
                
                if verbose:
                    print("iteration TDS")
                self.stdi = self.std
                
        if self.unit_S == "molal" and self.dil == 0:
            self.cat[1:6] = self.tot[1:6]
            self.nchcat[1:6] = nch[1:6]
            self.ani[1:4] = self.tot[6:9]
            self.nchani[1:4] = -nch[6:9]
            self.ani[4] = self.molal[12]
            self.ani[5] = self.molal[14] + self.molal[16] + self.molal[17]
            self.nchani[4] = -nch[12]
            self.nchani[5] = -nch[14]
            
            self.density()  
            
    def iterate_pco2(self, verbose=True):
        
        # calculate pco2 and ph
        self.calculate_pCO2(verbose=verbose)
        
        # if user specifies a pco2, redo the calculations
        if self.pco2_S == "y":
                
            while np.abs(self.po - self.pc) > 0.01:
                if self.po < self.pc and self.poa < self.pc:
                    self.ph -= self.dph
                if self.po < self.pc and self.poa > self.pc:
                    self.dph = self.dph / 2
                    self.ph -= self.dph
                if self.po > self.pc and self.poa > self.pc:
                    self.ph += self.dph
                if self.po > self.pc and self.poa < self.pc:
                    self.dph = self.dph/2
                    self.ph += self.dph
                self.poa = self.po
                
                self.iterate_molalities(verbose=verbose)
                
                self.calculate_pCO2(verbose=verbose)
        
        # if pk is greater than pkf, 
        if self.pk > self.pkf:
            self.pk = self.pkf
            
            if verbose:
                print("last iteration")
            
            self.iterate_molalities(verbose=verbose)
            self.calculate_pCO2(verbose=verbose)
            
        
    def calculate_pCO2(self, verbose=True):
        
        self.po = np.log10(self.act[15] / self.psc[14])
        
        if self.pco2_S == "":
            if self.diltot == 1:
                self.poinit = self.po
                self.phinit = self.ph
            
            if verbose:
                print("LOG PCO2 = {}".format(self.po))
            
            if self.pc:
                self.po0 = self.po
                self.ph0 = self.ph
                self.pco2_S = "y"
                
            elif not self.pc:
                self.pco2_S = "n"
                
        if self.pco2_S == "y":
            if verbose:
                    print("\nLog(PCO2) selected = {}".format(self.pc))
                    print("Log(PCO2) calculated = {}\n".format(self.po))
                
                
    def dilute_solution(self, dil=1, verbose=True, output=True):
        
        if dil > 1:
            self.diltot += dil
            self.pco2_S = ""
            self.pk = self.pk0
            self.dph = 0.2
            
            self.tot[np.arange(0, 13) != 11] = (
                self.tot[np.arange(0, 13) != 11] / dil)
            
            self.cat[1:6] = self.tot[1:6]
            self.nchcat[1:6] = nch[1:6]
            self.ani[1:4] = self.tot[6:9]
            self.nchani[1:4] = -nch[6:9]
            
            self.ani[4] = self.molal[12] / dil
            self.ani[5] = (self.molal[14] + self.molal[16] + 
                           self.molal[17]) / dil
            self.nchani[4] = -nch[12]
            self.nchani[5] = -nch[14]
            
            self.unit_S = "molal"
            
            self.density()
            
            # Recalculate the molalities after dilution
            self.calculate_molalities(verbose=verbose)
        
            # calculate activity coefficients // 500 Loop
            self.iterate_molalities(verbose=verbose)
            
            self.iterate_pco2(verbose=verbose)
            
            self.calculate_alkalinities()
            
            self.print_screen(verbose=verbose, output=output)
            
            self.saturation_state(verbose=verbose)
            
            if verbose:
                print("The initial solution has been "\
                      "diluted {} times".format(self.diltot))

    def invar(self):
        self.kinvariant = 0
        ncm = 14
        nbmin = 10
        ninvar = nbmin + 3
        kinv = np.zeros(ninvar+1)
        minv_S = np.empty(ninvar+1, dtype=np.object)
        psminv = np.zeros(ninvar+1)
        winv = np.zeros((ninvar+1, ncm+1))
        minvar_S = np.empty(ninvar+1, dtype=np.object)
        psminvar = np.zeros(ninvar+1)
        t0 = np.zeros((ninvar+1, ninvar+1))
        t1 = np.zeros((ninvar+1, ninvar+1))
        t2 = np.zeros((ninvar+1, ninvar+1))
        t3 = np.zeros((ninvar+1, ncm+1))
        t4 = np.zeros((ncm+1, ncm+1))
        tt4 = np.zeros(ncm+1)
        for k in range(1, 4):
            psminv[k] = np.log10(self.psc[k])
        winv[1, 11] = 1
        winv[1, 13] = 1
        winv[1, 0] = -1
        winv[2, 11] = 1
        winv[2, 12] = -1
        winv[2, 14] = 1
        winv[3, 11] = 1
        winv[3, 12] = 1
        winv[3, 0] = -1
        n1 = 3
        for k in range(1, nm+1):
            if lmin[k] == 1:
                n1 += 1
                kinv[n1] = k
                minv_S[n1] = mineral_S[k]
                psminv[n1] = np.log10(self.psol[k])
                for j in range(0, ncm+1):
                    winv[n1, j] = wmin[k, j]
        for i in range(1, n1+1):
            winv[i, 0], winv[i, 14] = winv[i, 14], winv[i, 0]
        for i in range(1, n1+1):
            for j in range(i, n1+1):
                t1[i, j] = 0
                for k in range(0, ncm):
                    t1[i, j] = t1[i, j] + winv[i, k] * winv[j, k]
                    t1[j, i] = t1[i, j]
                    t0[i, j] = t1[i, j]
                    t0[j, i] = t0[i, j]
        for k in range(2, n1+1):
            for i in range(k, n1+1):
                if np.abs(t1[i, k-1]) > epsilon:
                    u = t1[k-1, k-1] / t1[i, k-1]
                    for j in range(k, n1+1):
                        t1[i, j] = t1[k-1, j] - t1[i, j] * u
                        if np.abs(t1[i, j]) < epsilon:
                            t1[i, j] = 0
        det = 1
        for i in range(1, n1+1):
            if np.abs(t1[i, i]) < epsilon:
                det = 0
                break
        if det == 0:
            n3 = 0
            n2 = n1 - 1
            for kk in range(1, n1+1):
                ii = 0
                for i in range(1, n1+1):
                    if i != kk:
                        ii += 1
                        jj = 0
                        for j in range(1, n1+1):
                            if j != kk:
                                jj += 1
                                t2[ii, jj] = t0[i, j]
                for k in range(2, n2+1):
                    for i in range(k, n2+1):
                        if np.abs(t2[i, k-1]) > epsilon:
                            u = t2[k-1, k-1] / t2[i, k-1]
                            for j in range(k, n2+1):
                                t2[i, j] = t2[k-1, j] - t2[i, j] * u
                                if np.abs(t2[i, j]) < epsilon:
                                    t2[i, j] = 0
                det1 = 1
                for i in range(1, n2+1):
                    if np.abs(t2[i, i]) < epsilon:
                        det1 = 0
                        break
                if det1 == 1:
                    n3 += 1
                    self.kinvar[n3] = kinv[kk]
                    minvar_S[n3] = minv_S[kk]
                    psminvar[n3] = psminv[kk]
                    for j in range(0, ncm+1):
                        t3[n3, j] = winv[kk, j]
            if n3 == 0:
                self.kinvariant = -1
            elif n3 > 0:
                n4 = ncm
                for j in range(ncm, 0, -1):
                    u = 0
                    for i in range(1, n3+1):
                        u += t3[i, j] ** 2
                    if u < epsilon:
                        for k in range(j+1, n4+1):
                            for i in range(1, n3+1):
                                t3[i, k-1] = t3[i, k]
                        n4 -= 1
                for i in range(1, n4+1):
                    for j in range(1, n4+1):
                        t4[i, j] = 0
                        for k in range(1, n3+1):
                            t4[i, j] = t4[i, j] + t3[k, i] * t3[k, j]
                            t4[j, i] = t4[i, j]
                for i in range(1, n4+1):
                    tt4[i] = 0
                    for k in range(1, n3+1):
                        tt4[i] = tt4[i] + t3[k, i] * psminvar[k]
                for k in range(2, n4+1):
                    for i in range(k, n4+1):
                        if np.abs(t4[i, k-1]) > epsilon:
                            u = t4[k-1, k-1] / t4[i, k-1]
                            for j in range(k, n4+1):
                                t4[i, j] = t4[k-1, j] - t4[i, j] * u
                                if np.abs(t4[i, j]) < epsilon:
                                    t4[i, j] = 0
                            tt4[i] = tt4[k-1] - tt4[i] * u
                if np.abs(t4[n4, n4]) > epsilon:
                    ah2o = 10 ** (tt4[n4] / t4[n4, n4])
                    if ah2o > 1 or ah2o <= 0.01:
                        self.kinvariant = -2
                    else:
                        self.kinvariant = n3
                        for i in range(self.kinvariant, 0, -1):
                            if self.kinvar[i] == 0:
                                for k in range(1, self.kinvariant):
                                    self.kinvar[k] = self.kinvar[k+1]
                                self.kinvariant -= 1
                elif np.abs(t4[n4, n4]) <= epsilon:
                    self.kinvariant = -2

    def modify_database(self, add_min=None, rem_min=None, verbose=True, 
                        output=False):
        
        if add_min:
            if type(add_min) is not list:
                if type(add_min) is str:
                    add_min = [add_min]
                else:
                    add_min = list(add_min)
            add_min = [i.upper() for i in add_min]
            am = [mineral_S.index(add_min[i]) for i in 
                  range(0, len(add_min)) if add_min[i] in mineral_S]
            nwmin[am] = 1
            if verbose:
                print("ADDING MINERALS")
                for i in add_min:
                    if i in mineral_S:
                        print("{} added".format(i))
                    else:
                        print("{} could not be found "\
                              "in the database".format(i))
                print()
        
        if rem_min:
            if type(rem_min) is not list:
                if type(rem_min) is str:
                    rem_min = [rem_min]
                else:
                    rem_min = list(rem_min)
            rem_min = [i.upper() for i in rem_min]
            rm = [mineral_S.index(rem_min[i]) for i in 
                  range(0, len(rem_min)) if rem_min[i] in mineral_S]
            nwmin[rm] = 0
            if verbose:
                print("REMOVING MINERALS")
                for i in rem_min:
                    if i in mineral_S:
                        print("{} removed".format(i))
                    else:
                        print("{} could not be found "\
                              "in the database".format(i))
                print()

        if add_min or rem_min:
            murtf2 = read_file("murtf2")
            (nc, na, a) = (int(i) for i in murtf2[0])
            species = murtf2[1:nc+na+2]
            minerals = murtf2[nc+na+2:]
            minerals_on = [i for i in minerals if i[0] in 
                           np.array(mineral_S)[nwmin == 1].tolist()]
            nmnew = len(minerals_on)

            with open("murtf0", "w") as file:
                file.write(",".join([str(i) for i in (nc, na, nmnew)]))
                file.write('\n')
                for i in species:
                    file.write(",".join(i))
                    file.write('\n')
                for i in minerals_on:
                    file.write(",".join(i))
                    file.write('\n')
                file.close()
                
                self.print_screen(verbose=verbose, output=output)
            
                self.saturation_state(verbose=verbose)
        
    def saturation_state(self, verbose=True):
        self.kinvar = np.zeros(nt+1)
        
        for k in range(1, nm+1):
            lmin[k] = 0
        
        if verbose:
            print("The initial solution is oversaturated in {} mineral(s)" \
                  " of the data base MURTF3:".format(self.nwm))
            print()
            
            oversat = []
            for k in range(1, nm+1):
                if self.pai[k] / self.psol[k] >= 1 and nwmin[k] == 1:
                    oversat.append([mineral_S[k], self.pai[k] / self.psol[k]])
                    lmin[k] = 1
            oversat = pd.DataFrame(oversat)
            with pd.option_context('display.max_rows', None, 
                                               'display.max_columns', None):
                print(oversat.to_string(index=False, header=False))
            print()
            
        if self.nwm > self.nminer:
            print("VIOLATION OF THE PHASE RULE:")
            print("The maximum number of minerals " \
                  "allowed is {}".format(self.nminer))
            print("The evaporation program cannot start with this paragenesis")
            print()
        elif self.nwm <= self.nminer:
            self.invar()
            if self.kinvariant > 0:
                print("The activity of water is constrained by:")
                print()
                for k in range(1, self.kinvariant+1):
                    print(mineral_S[self.kinvar[k]])
                print()
            elif self.kinvariant == -1:
                print("System in thermodynamic desequilibrium")
                print("The activity of water is constrained at "\
                      "different values")
                print("by more than one mineral assemblage")
                print()
            elif self.kinvariant == -2:
                print("System in thermodynamic desequilibrium:")
                print("inconsistant mineral assemblage")
                print()
            elif self.kinvariant == 0:
                print("No invariant paragensis detected")
                print()
            
            if self.kinvariant != 0:
                print("The evaporation program cannot start with this "\
                      "paragenesis")
                print()
                    
            if self.kinvariant == 0 and self.nwmp > 0:
                print("The solution is close to saturation in {} mineral(s) "\
                      "of the data base MURTF3:".format(self.nwmp))
                close_saturation = []
                for k in range(1, nm+1):
                    if (self.pai[k] / self.psol[k] >= 0.9 and 
                        self.pai[k] / self.psol[k] < 1 and nwmin[k] == 1):
                        
                        close_saturation.append([mineral_S[k], self.pai[k] / self.psol[k]])
                        lmin[k] = 1
                        
                close_saturation = pd.DataFrame(close_saturation)
                with pd.option_context('display.max_rows', None, 
                                               'display.max_columns', None):
                    print(close_saturation.to_string(index=False, header=False))
                print()
                
                self.invar()
                
                if self.kinvariant > 0:
                    print("At the start of evaporation, the activity of "\
                          "water may be constrained by: ")
                    
                    for k in range(1, self.kinvariant+1):
                        print(mineral_S[self.kinvar[k]])
                    print()
                elif self.kinvariant == -1:
                    print("System in thermodynamic desequilibrium")
                    print("The activity of water is constrained at "\
                          "different values")
                    print("by more than one mineral assemblage")
                elif self.kinvariant == -2:
                    print("System in thermodynamic desequilibrium:")
                    print("inconsistant mineral assemblage")
                    print()
                elif self.kinvariant == 0:
                    print("No invariant paragensis detected")
                    print()
                if self.kinvariant != 0:
                    print("If the evaporation program does not start")
                    print("slightly dilute the solution again")
                    
                

    def calculate_alkalinities(self):
        self.alcar = self.molal[12] + 2 * (self.molal[14] + self.molal[16] + self.molal[17])
        self.albor = self.molal[19] + self.molal[20] + 2 * self.molal[21] + self.molal[22] + self.molal[23]
        self.alsil = self.molal[24]
        self.aloh = self.molal[13] + self.molal[18]
        self.alh = -1 * self.molal[11] - self.molal[25]
        self.altest = self.alcar + self.albor + self.alsil + self.aloh + self.alh
        
    def print_screen(self, verbose=True, output=False):
        outfile = "{}.log".format(self.label)
        if output:
            with open(outfile, "r+") as file:
                file.truncate(0)
                file.close()
        
        for i in range(1, n0+1):
            self.totinit[i] = 0
            for j in range(1, n+1):
                self.totinit[i] += self.molal[j] * kmat[i, j]
                
        df1 = pd.DataFrame({"molality": self.molal[1:n0+1], 
                            "act coeff": self.gact[1:n0+1],
                            "activity": self.act[1:n0+1],
                            "molal tot": self.totinit[1:n0+1]},
                           index =  aq_S[1:n0+1])
        df1 = df1.loc[df1["molality"] != 0]
        df2 = pd.DataFrame({"molality": self.molal[n0+1:n+1],
                            "act coeff": self.gact[n0+1:n+1],
                            "activity": self.act[n0+1:n+1]},
                           index = aq_S[n0+1:n+1])
        df2 = df2.loc[df2["molality"] != 0]
        
        df = pd.concat([df1, df2])
        
        
        if verbose:
            with pd.option_context('display.max_rows', None, 
                                   'display.max_columns', None):
                print(df.to_string(na_rep=""))
                print()
        
        if output:
            with open(outfile, "a") as file:
                file.write(df.to_string(na_rep=""))
                file.write("\n\n")
                file.close()

        text = ["ELECTRICAL BALANCE     = {} % corrected on {} and {}".format(self.dca, aq_S[self.icat], aq_S[self.iani]),
                "TOTAL DISSOLVED SOLIDS = {} g/kg(H2O)".format(self.std),
                "MOLAL/MOLAR FACTOR     = {}".format(1000 * self.dens / (1000 + self.std)),
                "DENSITY                = {}".format(self.dens),
                "IONIC STRENGTH         = {}".format(self.fi),
                "WATER ACTIVITY         = {}".format(self.aw),
                "CARBONATE ALKALINITY   = {}".format(self.alcar),
                "BORATE ALKALINITY      = {}".format(self.albor),
                "SILICATE ALKALINITY    = {}".format(self.alsil),
                "OH ALKALINITY          = {}".format(self.aloh),
                "H ALKALINITY           = {}".format(self.alh),
                "\nTOTAL ALKALINITY       = {}".format(self.altest),
                "init. alk.             = {}".format(self.tot[12])]
        
        if self.pco2_S == "" or self.pco2_S == "n":
            text.insert(6, "pH                     = {}".format(self.ph))
            text.insert(7, "LOG PCO2               = {}".format(self.po))
        elif self.pco2_S == "y":
            text.insert(6, "INITIAL LOG PCO2       = {}".format(self.poinit))
            text.insert(7, "INITIAL pH             = {}".format(self.phinit))
            text.insert(8, "CALCULATED LOG PCO2    = {}".format(self.po))
            text.insert(9, "CALCULATED pH          = {}".format(self.ph))
        
        if self.diltot > 1:
            text.insert(6, "DILUTION               = {}".format(self.diltot))
            
        if verbose:
            for line in text:
                print(line)
            print()
            
        if output:
            with open(outfile, "a") as file:
                for line in text:
                    file.write(line)
                    file.write("\n")
                file.write("\n")
                file.close()

        condition = np.where(self.ica==1)[0][np.where(np.where(self.ica==1)[0] <= 12)]
        u = np.zeros(condition.shape[0])
        for k in range(0, u.shape[0]):
            i = condition[k] * np.ones(n, dtype=int)
            j = np.arange(1, n+1)
            u[k] = np.sum(self.molal[j] * kmat[i, j])
            
        tests = self.tot[condition]
        d = np.zeros(condition.shape[0])
        d[u + tests != 0] = 200 * np.abs(u - tests) / (u + tests)
        names = np.char.upper(aq_S[condition])
        names[-1] = "ALK"
        
        df = pd.DataFrame({"TESTS": names, "SUM OF SPECIES": u, 
                           "INIT CONC.": tests, "BALANCE %": d})
        
        if verbose:
            with pd.option_context('display.max_rows', None, 
                                               'display.max_columns', None):
                print(df.to_string(index=False))
            print()
        
        if output:
            with open(outfile, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()
            
        data_lst = []
        
        for i in range(13, n+1):
            u = 0
            if self.ica[i] == 1:
                for j in range(0, n+1):
                    if self.act[j] != 0:
                        u += kmat[i, j] * np.log10(self.act[j])
                v = np.log10(self.psc[i-12])
                d = 200 * np.abs(u - v) / (u + v)
                if i == 13:
                    zone_S = "h2o"
                else:
                    zone_S = aq_S[i].lower()
                data_lst.append([zone_S, u, v, d])
        
        df = pd.DataFrame(data_lst, columns=["", "log(IAP)", 
                                             "log(K)", "balance %"])
        
        if verbose:
            with pd.option_context('display.max_rows', None, 
                                               'display.max_columns', None):
                print(df.to_string(index=False))
            print()
            
        if output:
            with open(outfile, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()
            
        self.nwm = 0
        self.nwmp = 0
        
        self.pai = np.zeros(nm+1)
        data_lst = []
        
        for k in range(1, nm+1):
            self.pai[k] = 1
            for i in range(0, nt+1):
                self.pai[k] = self.pai[k] * self.act[i] ** wmin[k, i]
            if self.pai[k] != 0:
                if nwmin[k] == 0:
                    zone_S = " " + mineral_S[k].lower()
                elif nwmin[k] == 1:
                    zone_S = "*" + mineral_S[k].lower()
                x_S = " "
                if self.pai[k] / self.psol[k] >= 1 and nwmin[k] == 1:
                    self.nwm += 1
                    x_S = "*"
                elif self.pai[k] / self.psol[k] >= 0.9 and self.pai[k] / self.psol[k] < 1 and nwmin[k] == 1:
                    self.nwmp += 1
                data_lst.append([zone_S, self.psol[k], self.pai[k], self.pai[k] / self.psol[k], x_S])
        
        df = pd.DataFrame(data_lst, columns=["", "solub prod", "ion act prod", 
                                             "satur ratio", ""])
                
        if verbose:
            with pd.option_context('display.max_rows', None, 
                                               'display.max_columns', None):
                print(df.to_string(index=False))
            print()
        
        if output:
            with open(outfile, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()
                
        if verbose and output:
            print("LOG FILE IS {}".format(outfile))
            print()
            
            
            
    def run_eql(self, syst_S, unit_S, dil=1, add_min=None, rem_min=None, 
                pkmol=None, pkeq=None, incr=0, print_step=1, output_step=1, 
                storage_step=1, stdmax=0, verbose=True, output=False, 
                classic=False):
        
        if output_step == 0:
            p_S = "n"
        else:
            p_S = "y"
            
        if verbose:
            print("\nThis is EQL..............\n")
            print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')
        
        # initial charge balance
        self.charge_balance(verbose=verbose)
    
        # intial molality calculation
        self.calculate_molalities(verbose=verbose)
        
        # calculate activity coefficients // 500 Loop
        self.iterate_molalities(verbose=verbose)
        
        self.iterate_pco2(verbose=verbose)
        
        self.calculate_alkalinities()
        
        self.print_screen(verbose=verbose, output=output)
        
        self.saturation_state(verbose=verbose)
        
        self.dilute_solution(dil=dil)
        
        # Modify mineral database
        self.modify_database(add_min=add_min, rem_min=rem_min, verbose=verbose)
        if add_min or rem_min:
            min_S = "murtf0"
        # Modify convergence limits
        if pkmol:
            self.pkmol = pkmol
        else:
            self.pkmol = 0.001
        
        if pkeq:
            self.pkeq = pkeq
        else:
            self.pkeq = .0000000000001
        
        # Set files for output
        self.chem_file = "{}.j{}&".format(self.label, syst_S)
        self.event_file = "{}.j{}@".format(self.label, syst_S)
        self.min_file = "{}.j{}%".format(self.label, syst_S)
        self.transfer_file = "{}.tra".format(self.label)
        
        # Add heading to chem_file
        
        self.constit_S = ["label", "fc", "eva", "ds", "ph", "alk"]
        
        for i in range(1, 9):
            if self.tot[i] > 0:
                self.constit_S.append(aq_S[i])
        if self.tot[9] != 0:
            self.constit_S.append('b')
        if self.tot[10] != 0:
            self.constit_S.append('si')
        self.constit_S.append('tds')
            
        # Write transfer file name to stockage file
        with open("stockage", "w") as file:
            file.write(self.transfer_file)
            file.close()
        
        # Write transfer information to the transfer 
        lines =[self.temp, self.tinit, self.ph, self.phinit, self.po, 
                self.poinit, self.diltot]
        lines.append(",".join(self.constit_S))
        for i in range(1, 11):
            lines.append(self.tot[i])
        lines.append(self.molal[15])
        for i in range(1, 11):
            lines.append(self.molal[i])
        lines.append(self.mh2o)
        lines.append(self.molal[13])
        lines.append(self.molal[11])
        lines.append(self.molal[14])
        lines.append(self.molal[12])
        for i in range(16, 26):
            lines.append(self.molal[i])
        lines.append(syst_S)
        lines.append(incr)
        lines.append(print_step)
        lines.append(p_S)
        lines.append(output_step)
        lines.append(storage_step)
        lines.append(unit_S)
        lines.append(self.chem_file)
        lines.append(self.event_file)
        lines.append(self.min_file)
        lines.append(min_S)
        lines.append(stdmax)
        lines.append(self.pkmol)
        lines.append(self.pkeq)
        
        lines = [str(i) for i in lines]
        
        with open(self.transfer_file, "w") as file:
            file.write("\n".join(lines))
            file.close()
            
        # Run EVP
        self.run_evp()
        
    def run_evp(self, verbose=True):
        pass
        if verbose:
            print("\nThis is EVP..............")
            print("starting the evaporation program")
            print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')
            
        
        
        
    
    
    
test = Water(label='test', temp=30, dens=1, ph=6.55, na=84.5, k=3.3, li=0,
             ca=2.7, mg=1.3, cl=39.5, so4=0, alk=56.2, no3=0, si=0, b=0, 
             pc=-3)

am = ['CALCITE']
rm = ['BRUCITE', 'MAGNESITE', 'HYDROMAGNESITE', 'NESQUEHONITE',
      'ANTARCTICITE', 'ARAGONITE', 'BURKEITE', 'GLASERITE', 'DOLOMITE']

test.run_eql(syst_S='c', unit_S="molal", dil=1, add_min=am, rem_min=rm, verbose=True)

test.charge_balance()
test.calculate_molalities()
# RUN EVP
# parameters
