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

n = 25
ntot = 12
n0 = 10
max_cat = 5
max_an = 5



repinit_S, fich_S, molemin_S, fichmin_S = ".", ".", ".", "."
nchcat, nchani = 0, 0


kmat = np.zeros((n+1, n+1))


aqu = pd.read_csv("aqu.dat", header=None)
aq_S = aqu.iloc[:,0].values.astype(str)
atom = aqu.iloc[:,1].values
nch = aqu.iloc[:,2].values

kmat[1:, 1:] = pd.read_csv("matrice1", header=None).values
kmat[13, 0], kmat[15, 0], kmat[20, 0], kmat[21, 0] = -1, -1, 3, 5

complex3 = pd.read_csv("complex3", header=None)
complex3.columns = ['species', 'at', 'bt', 'ct', 'dt', 'et']

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
min_db = np.zeros((nm, 5))

for k in range(0, nm):
    ncomp = int(mineral_list[k][1])
    c_ion = np.zeros(ncomp)
    nom_ion_S = np.empty(shape=ncomp, dtype=np.object)
    for i in range(0, ncomp):
        c_ion[i] = float(mineral_list[k][2+2*i])
        nom_ion_S[i] = mineral_list[k][3+2*i]

        j = np.where(aq_S == nom_ion_S[i].lower())
        wmin[k+1, j] = c_ion[i]
    min_db[k, :] = [float(i) for i in mineral_list[k][(1+ncomp)*2:]]
    
min_db = pd.DataFrame(min_db, columns=['at', 'bt', 'ct', 'dt', 'et'])
min_db['species'] = mineral_S[1:]


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
        self.mu = 10 ** (ion_db['at'] + 
                         ion_db['bt'] / 300 * self.temp + 
                         ion_db['ct'] / 30000 * self.temp ** 2 + 
                         ion_db['dt'] / 3000000 * self.temp ** 3 + 
                         ion_db['et'] / 300000000 * self.temp ** 4).values
        self.mum = np.zeros(nm+1)
        self.mum[1:] = (min_db['at'] + 
                        min_db['bt'] / 300 * self.temp +
                        min_db['ct'] / 30000 * temp ** 2 +
                        min_db['dt'] / 3000000 * temp ** 3 +
                        min_db['et'] / 300000000 * temp ** 4).values
        
        # for k in range(1, nm+1):
        #     u = self.mum[k]
        #     for i in range(0, nt+1):
        #         u += wmin[k, i] * self.mu[i]
                
        # self.psol = np.exp(u)

        self.psol = np.exp(self.mum - np.sum((wmin[:,:] * self.mu), axis=1))
        self.gact = np.zeros(n+1)
        self.molal = np.zeros(n+1)
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
            dca = 200 * np.abs(sc - sa) / (sc + sa)
            delta = sc - sa
            self.tot[icat] = self.tot[icat] - delta / 2 / nch[icat]
            self.tot[iani] = self.tot[iani] + delta / 2 / -(nch[iani])
            print("electrical balance = {}%".format((dca * 100 + 0.5) / 100))
              
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
                    icat = i
        
        sa, amax = 0, 0
        for i in range(1, n+1):
            if self.nch[i] < 0:
                sa += self.molal[i] * -self.nch[i]
                if self.molal[i] * -self.nch[i] > amax:
                    amax = self.molal[i] * -self.nch[i]
                    iani = i
        
        delta = sc - sa
        self.molal[icat] = self.molal[icat] - delta / 2 / nch[icat]
        self.molal[iani] = self.molal[iani] + delta / 2 / (-nch[iani])
        
        sc = (self.molal[1] + self.molal[2] + self.molal[3] + 
              self.molal[4] * 2 + self.molal[5] * 2 + self.molal[11] + 
              self.molal[18] + self.molal[22] + self.molal[23])
        sa = (self.molal[6] + self.molal[7] * 2 + self.molal[8] + 
              self.molal[12] + self.molal[13] + self.molal[14] * 2 +
              self.molal[19] + self.molal[20] + self.molal[21] * 2 +
              self.molal[24] + self.molal[25])
        
        if verbose:
            print("\nSum of cations = {} corrected for {}".format(sc, 
                                                                  aq_S[icat]))
            print("Sum of anions = {} corrected for {}\n".format(sa, 
                                                                 aq_S[iani]))
            
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

        fi = u / 2
        fj = np.sqrt(fi)
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
                                self.j0(xc[j, j]) / 2) / fi / 2)
                    fc[i, j] = ((xc[i, j] * self.j1(xc[i, j]) - xc[i, i] * 
                                self.j1(xc[i, i]) / 2 - xc[j, j] * 
                                self.j1(xc[j, j]) / 2) / fi ** 2 / 4 - 
                                ec[i, j] / fi)
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
                                self.j0(xa[j, j]) / 2) / fi / 2
                    fa[i, j] = ((xa[i, j] * self.j1(xa[i, j]) - xa[i, i] * 
                                self.j1(xa[i, i]) / 2 - xa[j,j] * 
                                self.j1(xa[j, j]) / 2) / fi ** 2 / 4 - 
                                ea[i, j] / fi)
                    ea[j, i] = ea[i, j]
                    fa[j, i] = fa[i, j]
        
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                pp[i, j] = fc[i, j]
                p[i, j] = self.tc[i, j] + ec[i, j]
                pf[i, j] = p[i, j] + pp[i, j] * fi
                pp[j, i] = pp[i, j]
                p[j, i] = p[i, j]
                pf[j, i] = pf[i, j]
        
        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                qp[i, j] = fa[i, j]
                q[i, j] = self.ta[i, j] + ea[i, j]
                qf[i, j] = q[i, j] + qp[i, j] * fi
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
                            fi + self.b2[i, j] * (self.g1(w)) / fi)
        
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
            
        u = -self.ap0 * fi ** 1.5e0 / (1 + self.bp0 * fj)
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
                
                
    def dilute_solution(self, dil=1, verbose=True):
        
        if dil > 1:
            self.diltot += dil
            self.pco2_S = ""
            self.pk = self.pk0
            self.dph = 0.2
            
            self.tot[np.arange(1, 13) != 11] = (
                self.tot[np.arange(1, 13) != 11] / dil)
            
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
         
    def run_eql(self, verbose=True):
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

with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    pass
     
test = Water(label='test', temp=30, dens=1, ph=6.55, na=84.5, k=3.3, li=0,
             ca=2.7, mg=1.3, cl=39.5, so4=0, alk=56.2, no3=0, si=0, b=0, pc=-6.73)

#test.charge_balance()
#test.calculate_molalities()
#test.iterate_molalities()


#test.calculate_pCO2()
#test.iterate_pco2()

#test.iterate_molalities()
#test.calculate_pCO2()

test.run_eql()