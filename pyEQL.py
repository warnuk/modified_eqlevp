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

n = 25
ntot = 12
n0 = 10
max_cat = 5
max_an = 5

print("\nThis is EQL..............\n")
print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')

repinit_S, fich_S, molemin_S, fichmin_S = ".", ".", ".", "."
pk, pk0, pkf, pkstd, dph, dil, diltot = .1e0, .1e0, .0001e0, .1e0, .2e0, 0., 1.
mh2o, pco2_S, nconv, eps = 55.51e0, "", 2, 1e-12
nchcat, nchani = 0, 0


kmat = np.zeros((n+1, n+1))
ica = np.zeros(n+1)

aqu = pd.read_csv("aqu.dat", header=None)
aq_S = aqu.iloc[:,0].values.astype(str)
atom = aqu.iloc[:,1].values
nch = aqu.iloc[:,2].values

kmat[1:, 1:] = pd.read_csv("matrice1", header=None).values
kmat[13, 0], kmat[15, 0], kmat[20, 0], kmat[21, 0] = -1, -1, 3, 5

complex3 = pd.read_csv("complex3", header=None)
complex3.columns = ['species', 'at', 'bt', 'ct', 'dt', 'et']

with open("murtf2", 'r+') as file:
    text = [i.split(',') for i in file.read().split('\n')][:-1]
    file.close()

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


# LINE 260
with open("murtf3", 'r+') as file:
    text = [i.split(',') for i in file.read().split('\n')][:-1]
    file.close()
nc, na, nm0 = [int(i) for i in text[0]]

#ion0 = text[1:(2+nc+na)]
mineral0_S = [i[0] for i in text[(2+nc+na):]]


nwmin[np.in1d(mineral_S, mineral0_S)] = 1

min_S = "murtf3"

#500 continue !dilution (LINE 281)






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
    
    return(Water(label, temp, dens, ph, na, k, li, ca, mg, cl, so4, no3, b, si, alk))
    

    

class Water:
    'Class for input waters to be used in EQL/EVP simulations'
    
    def __init__(self, label, temp, dens, ph, na, k, li, 
                 ca, mg, cl, so4, no3, b, si, alk, unit_S="molal"):
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
        self.tot[11] = 10 ** -ph
        self.tot[12] = self.alk / 1000
        
        self.nminer = np.count_nonzero(self.tot[1:11])
        
        self.charge_balance()
        
        self.tot0 = self.tot
        
        
        self.psc = 10 ** (complex3['at'] + 
                          complex3['bt'] / 300 * self.temp + 
                          complex3['ct'] / 30000 * self.temp ** 2 + 
                          complex3['dt'] / 3000000 * self.temp ** 3 + 
                          complex3['et'] / 300000000 * self.temp ** 4).values
        
        self.mu = 10 ** (ion_db['at'] + 
                         ion_db['bt'] / 300 * self.temp + 
                         ion_db['ct'] / 30000 * self.temp ** 2 + 
                         ion_db['dt'] / 3000000 * self.temp ** 3 + 
                         ion_db['et'] / 300000000 * self.temp ** 4).values
        
        self.mum = (min_db['at'] + 
                    min_db['bt'] / 300 * self.temp +
                    min_db['ct'] / 30000 * temp ** 2 +
                    min_db['dt'] / 3000000 * temp ** 3 +
                    min_db['et'] / 300000000 * temp ** 4).values
        
        self.psol = np.exp(self.mum - np.sum((wmin[1:,:] * self.mu), axis=1)) 
        
        self.calculate_molalities()
        
        self.gact = np.zeros(n+1)
        
        self.ndepact = 0
        self.z = np.zeros((n, n))
        self.zz = np.zeros(n)
    
    def charge_balance(self):
        cat = np.array(np.where(nch[1:ntot+1] > 0)).flatten() + 1
        cat = cat[np.where(cat != 11)]
        cat_ch = self.tot[cat] * nch[cat]
        icat = np.argmax(cat_ch)
        sc = np.sum(cat_ch)
        
        ani = np.array(np.where(nch[1:ntot+1] < 0)).flatten() + 1
        ani_ch = self.tot[ani] * -nch[ani]
        iani = np.argmax(-ani_ch)
        sa = np.sum(ani_ch)
        
        print("Sum of cations = {}".format(sc))
        print("Sum of anions = {}".format(sa))
        
        if sc + sa != 0:
            dca = 200 * np.abs(sc - sa) / (sc + sa)
            delta = sc - sa
            self.tot[icat] = self.tot[icat] - delta / 2 / nch[icat+1]
            self.tot[iani] = self.tot[iani] + delta / 2 / -(nch[iani+1])
            print("electrical balance = {}%".format((dca * 100 + 0.5) / 100))
              
        
    def calculate_molalities(self):
        self.molal = np.zeros(n+1)
        self.molal[1] = self.tot[1]
        self.molal[2] = self.tot[2]
        self.molal[3] = self.tot[3]
        self.molal[6] = self.tot[6]
        self.molal[8] = self.tot[8]
        self.molal[11] = self.tot[11]
        self.molal[13] = self.psc[0] / self.molal[11]
        if self.tot[9] > 0:
            a = 1 + self.molal[13] / self.psc[6]
            b = 3 * self.psc[7] * self.molal[13]
            c = 4 * self.psc[8] * self.molal[13] ** 2
            xu = self.tot[9] / 2
            u = xu
            
            eq = a * xu + b * xu ** 3 + c * xu ** 4
            while (200 * abs(eq - self.tot[9]) / (eq + self.tot[9]) >= pk):
                u = u / 2
                if eq > self.tot[9]:
                    xu = xu - u
                else:
                    xu = xu + u
                
                eq = a * xu + b * xu ** 3 + c * xu ** 4
            
            self.molal[9] = xu
            self.molal[19] = (self.molal[9] * self.molal[13] / self.psc[6])
            self.molal[20] = (self.molal[13] * self.molal[9] ** 3 * 
                              self.psc[7])
            self.molal[21] = (self.molal[13] ** 2 * self.molal[9] ** 4 * 
                              self.psc[8])
        self.molal[14] = (self.tot[12] + self.molal[11] - self.molal[13] - 
                          self.molal[19] - self.molal[20] - 
                          2 * self.molal[21]) / (2 + self.molal[11] / 
                                                 self.psc[1])
        self.molal[12] = (self.tot[12] + self.molal[11] - self.molal[13] - 
                          self.molal[19] - self.molal[20] - 
                          2 * self.molal[21]) / (1 + 2 * self.psc[1] / 
                                                 self.molal[11])
                                                 
        self.molal[15] = self.molal[12] * self.molal[11] / self.psc[2]
        self.molal[4] = self.tot[4] / (1 + self.molal[14] / self.psc[3] +
                                       self.molal[19] / self.psc[9])
        self.molal[16] = self.molal[4] * self.molal[14] / self.psc[3]
        self.molal[22] = self.molal[4] * self.molal[19] / self.psc[9]
        self.molal[5] = self.tot[5] / (1 + self.molal[14] / self.psc[4] + 
                                       self.molal[13] / self.psc[5] + 
                                       self.molal[19] / self.psc[10])
        self.molal[17] = self.molal[5] * self.molal[14] / self.psc[4]
        self.molal[18] = self.molal[5] * self.molal[13] / self.psc[5]
        self.molal[23] = self.molal[5] * self.molal[19] / self.psc[10]
        self.molal[10] = self.tot[10] / (1 + self.psc[11] / self.molal[11])
        self.molal[24] = self.tot[10] / (1 + self.molal[11] / self.psc[11])
        self.molal[7] = self.tot[7] / (1 + self.molal[11] / self.psc[12])
        self.molal[25] = self.molal[7] * self.molal[11] / self.psc[12]

        cat = np.array(np.where(nch[1:n+1] > 0)).flatten() + 1
        cat_ch = self.molal[cat] * nch[cat]
        icat = np.where(self.molal == self.molal[cat][np.argmax(cat_ch)])[0]
        sc = np.sum(cat_ch)
        
        ani = np.array(np.where(nch[1:n+1] < 0)).flatten() + 1
        ani_ch = self.molal[ani] * -nch[ani]
        iani = np.where(self.molal == self.molal[ani][np.argmax(ani_ch)])[0]
        sa = np.sum(ani_ch)
        
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
        
        print("\nSum of cations = {} corrected for {}".format(sc, aq_S[icat][0]))
        print("Sum of anions = {} corrected for {}\n".format(sa, aq_S[iani][0]))
        
        s = np.sum(atom * self.molal)
        
        if self.unit_S == "molar":
            ee = 1000 / (1000 * self.dens - s)
            
            self.tot[np.delete(np.arange(0, ntot), 10)] = (ee * 
                                self.tot[np.delete(np.arange(0, ntot), 10)])
            
            self.molal[np.delete(np.arange(0, n), 10)] = (ee * 
                                self.molal[np.delete(np.arange(0, n), 10)])
        elif self.unit_S == "molal":
            ee = 1
        
        self.stdi = s * ee
    
    def actp(self):
        c = np.arange(0, 10)
        a = np.arange(0, 12)
        h = np.arange(0, 4)
        
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
            with open("coefft4", 'r+') as file:
                text = [i.split(",") for i in file.read().split("\n")][:-1]
                file.close()
        
            (nc, na, nn) = (int(i) for i in text[0])
            nzc, nza = np.zeros(nc+1), np.zeros(na+1)
            (b0, b1, b2, c0) = (np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)), 
                                np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)))
            sc = np.zeros((nc+1, nc+1, na+1))
            sa = np.zeros((na+1, na+1, nc+1))
            tc = np.zeros((nc+1, nc+1))
            ta = np.zeros((na+1, na+1))
            lc = np.zeros((nn+1, nc+1))
            la = np.zeros((nn+1, na+1))
            xi = np.zeros((nn+1, nc+1, na+1))
            
            
            
            nzc[1:] = [int(i[1]) for i in text[1:1+nc]]
            nza[1:] = [int(i[1]) for i in text[1+nc:1+nc+na]]
            
            (at, bt, ct, dt, et) = (float(i.replace('d', 'e')) 
                                    for i in text[1+nc+na][1:])
            
            self.ap0 = self.temperature(at, bt, ct, dt, et)
            
            # b0, b1, b2, c0
            
            
            index = 1+nc+na
            
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    b0[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    b1[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    b2[i, j] = self.temperature(at, bt, ct, dt, et)
                    
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    c0[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for i in range(1, nc):
                for j in range(i+1, nc+1):
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    tc[i, j] = self.temperature(at, bt, ct, dt, et)
                    tc[j, i] = tc[i, j]
                    
            for i in range(1, na):
                for j in range(i+1, na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    ta[i, j] = self.temperature(at, bt, ct, dt, et)
                    ta[j, i] = ta[i, j]
                    
            for k in range(1, nc):
                for i in range(k+1, nc+1):
                    for j in range(1, na+1):
                        index += 1
                        (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                                for k in text[index][1:])
                        sc[k, i, j] = self.temperature(at, bt, ct, dt, et)
                        sc[i, k, j] = sc[k, i, j]
            
            for k in range(1, na):
                for i in range(k+1, na+1):
                    for j in range(1, nc+1):
                        index += 1
                        (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                                for k in text[index][1:])
                        sa[k, i, j] = self.temperature(at, bt, ct, dt, et)
                        sa[i, k, j] = sa[k, i, j]
                        
            for i in range(1, nn+1):
                for j in range(1, nc+1):
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    lc[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for i in range(1, nn+1):
                for j in range(1, na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (float(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    la[i, j] = self.temperature(at, bt, ct, dt, et)
            
            for k in range(1, nn+1):
                for i in range(1, nc+1):
                    for j in range(1, na+1):
                        xi[k, i, j] = 0
                        
            xi[2, 9, 1] = -0.0102
            xi[2, 1, 2] = 0.046
        
        ec, ea = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        fc, fa = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        xc, xa = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        pp, qp = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        p, q = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        pf, qf = np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1))
        cc, bf = np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1))
        b, bp = np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1))
        gc, ga, gn = np.zeros(nc+1), np.zeros(na+1), np.zeros(nn+1)
        
        self.bp0 = 1.2e0
        self.mh2o = 55.51e0
        u, z = 0, 0
        
        for i in range(1, nc+1):
            u += c[i] * nzc[i] ** 2
            z += c[i] * nzc[i]
        for j in range(1, na+1):
            u += a[j] * nza[j] ** 2
            z += a[j] * nza[j]
        fi = u / 2
        fj = np.sqrt(fi)
        u = 6 * self.ap0 * fj
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                if nzc[i] == nzc[j]:
                    ec[i, j] = 0
                    fc[i, j] = 0
                else:
                    xc[i, j] = 2 * u
                    xc[i, i] = nzc[i] ** 2 * u
                    xc[j, j] = nzc[j] ** 2 * u
                    ec[i, j] = (self.j0(xc[i, j]) - self.j0(xc[i, i]) / 2 - 
                                self.j0(xc[j, j]) / 2) / fi / 2
                    fc[i, j] = ((xc[i, j] * self.j1(xc[i, j]) - xc[i, i] * 
                                self.j1(xc[i, i]) / 2 - xc[j, j] * 
                                self.j1(xc[j, j]) / 2) / fi ** 2 / 4 - 
                                ec[i, j] / fi)
                    ec[j, i] = ec[i, j]
                    fc[j, i] = fc[i, j]
                    
        for i in range(1, na):
            for j in range(i+1, na+1):
                if nza[i] == nza[j]:
                    ea[i, j] = 0
                    fa[i, j] = 0
                else:
                    xa[i, j] = 2 * u
                    xa[i, i] = nza[i] ** 2 * u
                    xa[j, j] = nza[j] ** 2 * u
                    ea[i, j] = (self.j0(xa[i, j]) - self.j0(xa[i, i]) / 2 -
                                self.j0(xa[j, j]) / 2) / fi / 2
                    fa[i, j] = ((xa[i, j] * self.j1(xa[i, j]) - xa[i, i] * 
                                self.j1(xa[i, i]) / 2 - xa[j,j] * 
                                self.j1(xa[j, j]) / 2) / fi ** 2 / 4 - 
                                ea[i, j] / fi)
                    ea[j, i] = ea[i, j]
                    fa[j, i] = fa[i, j]
        
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                pp[i, j] = fc[i, j]
                p[i, j] = tc[i, j] + ec[i, j]
                pf[i, j] = p[i, j] + pp[i, j] * fi
                pp[j, i] = pp[i, j]
                p[j, i] = p[i, j]
                pf[j, i] = pf[i, j]
        
        for i in range(1, na):
            for j in range(i+1, na+1):
                qp[i, j] = fa[i, j]
                q[i, j] = ta[i, j] + ea[i, j]
                qf[i, j] = q[i, j] + qp[i, j] * fi
                qp[j, i] = qp[i, j]
                q[j, i] = q[i, j]
                qf[j, i] = qf[i, j]
                
        w = fj * 12
        for i in range(1, nc+1):
            for j in range(1, na+1):
                cc[i, j] = c0[i, j] / np.sqrt(nzc[i] * nza[j]) / 2
                if nzc[i] == 2 and nza[j] == 2:
                    v = fj * 1.4e0
                if nzc[i] == 1 or nza[j] == 1:
                    v = fj * 2
                bf[i, j] = (b0[i, j] + b1[i, j] * np.exp(-v) + 
                            b2[i, j] * np.exp(-w))
                b[i, j] = (b0[i, j] + b1[i, j] * 
                           (self.g0(v)) + b2[i, j] * (self.g0(w)))
                bp[i, j] = (b1[i, j] * (self.g1(v)) / 
                            fi + b2[i, j] * (self.g1(w)) / fi)
        
        f = -self.ap0 * (fj / (1 + self.bp0 * fj) + 2 / self.bp0 * np.log(1 + self.bp0 * fj))
        for i in range(1, nc+1):
            for j in range(1, na+1):
                f += c[i] * a[j] * bp[i, j]
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                f += c[i] * c[j] * pp[i, j]
        for i in range(1, na):
            for j in range(i+1, na+1):
                f += a[i] * a[j] * qp[i, j]
        for ii in range(1, nc+1):
            u = nzc[ii] ** 2 * f
            for j in range(1, na+1):
                u += a[j] * (b[ii, j] * 2 + z * cc[ii, j])
            for i in range(1, nc+1):
                if i != ii:
                    v = 0
                    for j in range(1, na+1):
                        v += a[j] * sc[ii, i, j]
                    u += c[i] * (p[ii, i] * 2 + v)
            for i in range(1, na):
                for j in range(i+1, na+1):
                    u += a[i] * a[j] * sa[i, j, ii]
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    u += c[i] * a[j] * cc[i, j] * nzc[ii]
            for i in range(1, nn+1):
                u += h[i] * lc[i, ii] * 2
            for k in range(1, nn+1):
                for j in range(1, na+1):
                    u += h[k] * a[j] * xi[k, ii, j]
            gc[ii] = np.exp(u)
        for jj in range(1, na+1):
            u = nza[jj] ** 2 * f
            for i in range(1, nc+1):
                u += c[i] * (b[i, jj]) * 2 + z * cc[i, jj]
            for i in range(1, na+1):
                if i != jj:
                    v = 0
                    for j in range(1, nc+1):
                        v += c[j] * sa[jj, i, j]
                    u += a[i] * (q[jj, i] * 2 + v)
                        
            for i in range(1, nc):
                for j in range(i+1, nc+1):
                    u += c[i] * c[j] * sc[i, j, jj]
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    u += c[i] * a[j] * cc[i, j] * nza[jj]
            for j in range(1, nn+1):
                u += h[j] * la[j, jj]
            for k in range(1, nn+1):
                for i in range(1, nc+1):
                    u += h[k] * c[i] * xi[k, i, jj]
            ga[jj] = np.exp(u)
            
        for k in range(1, nn+1):
            u = 0
            for i in range(1, nc+1):
                u += c[i] * lc[k, i] * 2
            for j in range(1, na+1):
                u += a[j] * la[k, j] * 2
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    u += c[i] * a[j] * xi[k, i, j]
            gn[k] = np.exp(u)
        u = -self.ap0 * fi ** 1.5e0 / (1 + self.bp0 * fj)
        for i in range(1, nc+1):
            for j in range(1, na+1):
                u += c[i] * a[j] * (bf[i, j] + z * cc[i, j])
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                v = 0
                for k in range(1, na+1):
                    v += a[k] * sc[i, j, k]
                u += c[i] * c[j] * (pf[i, j] + v)
        for i in range(1, na):
            for j in range(i+1, na+1):
                v = 0
                for k in range(1, nc+1):
                    v += c[k] * sa[i, j, k]
                u += a[i] * a[j] * (qf[i, j] + v)
        for k in range(1, nn+1):
            for i in range(1, nc+1):
                u += h[k] * c[i] * lc[k, i]
        for k in range(1, nn+1):
            for j in range(1, na+1):
                u += h[k] * a[j] * la[k, j]
        for k in range(1, nn+1):
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    u += h[k] * c[i] * a[j] * xi[k, i, j]
        
        s = 0
        for i in range(1, nc+1):
            s += c[i]
        for j in range(1, na+1):
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
        self.gact[10] = self.aw * self.aw * gn[1] ** np.log(10e0)
        self.gact[9] = gn[2]
        self.gact[15] = gn[3]
        self.gact[16] = 1
        self.gact[17] = 1
        self.ndepact = 1
            
    def temperature(self, at, bt, ct, dt, et):
        return((at + bt * self.temp + ct * self.temp ** 2 + 
                dt * self.temp ** 3 + et * self.temp ** 4))
    
    def j0(self, x):
        ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
        j0 = x / (4. + ya * x ** yb * np.exp(yc * x ** yd))
        return(j0)
            
    def j1(self, x):
        ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
        j1 = (4. + ya * x ** yb * (1. - yb - yc * yd * x ** yd) * 
          np.exp(yc * x ** yd)) / (4. + ya * x ** yb * 
                                   np.exp(yc * x ** yd)) ** 2
        return(j1)
    
    def g0(self, x):
        g0 = 2. * (1. - (1. + x) * np.exp(-x)) / x ** 2
        return(g0)
            
    def g1(self, x):
        g1 = -2. * (1. - (1. + x + x ** 2 / 2.) * np.exp(-x)) / x ** 2
        return(g1)

    def run(self):
        # 200 continue !reprise
        # aka line 373
        self.nu = 1
        self.ncompt = 0
        while self.nu != 0:
            self.actp()
            self.act[1:] = self.molal * self.gact[1:]
            self.act[0] = self.aw
            self.tot[11] = (10 ** (-self.ph)) / self.gact[11]
            self.act[11] = 10 ** (-self.ph)
            for i in range(1, 13):
                u = 0
                for j in range(1, n+1):
                    if self.molal[j] != 0:
                        self.z[i, j] = kmat[i, j]
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
            for k in range(n, 1, -1):
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
            
            # LINE 453



test = Water(label='test', temp=30, dens=1, ph=6.55, na=84.5, k=3.3, li=0,
             ca=2.7, mg=1.3, cl=39.5, so4=0, alk=56.2, no3=0, si=0, b=0)


test.actp()




# b0, b1, b2, c0


# def temperature(temp, at, bt, ct, dt, et):
#     return((at + bt * temp + ct * temp ** 2 + 
#             dt * temp ** 3 + et * temp ** 4))


# import time

# ta = time.perf_counter_ns() 
 
# index = 2+nc+na
# coordinates = np.zeros((6, nc*na), dtype=int)
# coordinates[0,:] = np.repeat(np.arange(0, nc), na)
# coordinates[1,:] = np.tile(np.arange(0, na), nc)
# coordinates[2,:] = np.arange(0, nc*na*4, 4) + index
# coordinates[3,:] = np.arange(1, nc*na*4, 4) + index
# coordinates[4,:] = np.arange(2, nc*na*4, 4) + index
# coordinates[5,:] = np.arange(3, nc*na*4, 4) + index

# b0 = np.zeros((nc, na))
# b1 = np.zeros((nc, na))
# b2 = np.zeros((nc, na))
# c0 = np.zeros((nc, na))

# for i in range(0, nc*na):
#     x, y = coordinates[0:2, i]
    
#     (at, bt, ct, dt, et) = (float(k) for k in 
#                             text[coordinates[2,i]][1:])
#     b0[x, y] = temperature(test.temp, at, bt, ct, dt, et)
    
#     (at, bt, ct, dt, et) = (float(k) for k in 
#                             text[coordinates[3,i]][1:])
#     b1[x, y] = temperature(test.temp, at, bt, ct, dt, et)
    
#     (at, bt, ct, dt, et) = (float(k) for k in 
#                             text[coordinates[4,i]][1:])
#     b2[x, y] = temperature(test.temp, at, bt, ct, dt, et)
    
#     (at, bt, ct, dt, et) = (float(k) for k in 
#                             text[coordinates[5,i]][1:])
#     c0[x, y] = temperature(test.temp, at, bt, ct, dt, et)
    
# tb = time.perf_counter_ns() 

# b0 = np.zeros((nc, na))
# b1 = np.zeros((nc, na))
# b2 = np.zeros((nc, na))
# c0 = np.zeros((nc, na))

# index = 2+nc+na
# for i in range(0, nc):
#     for j in range(0, na):
#         (at, bt, ct, dt, et) = (float(k) for k in text[index][1:])
#         b0[i, j] = temperature(test.temp, at, bt, ct, dt, et)
        
#         index += 1
#         (at, bt, ct, dt, et) = (float(k) for k in text[index][1:])
#         b1[i, j] = temperature(test.temp, at, bt, ct, dt, et)
        
#         index += 1
#         (at, bt, ct, dt, et) = (float(k) for k in text[index][1:])
#         b2[i, j] = temperature(test.temp, at, bt, ct, dt, et)
        
#         index += 1
#         (at, bt, ct, dt, et) = (float(k) for k in text[index][1:])
#         c0[i, j] = temperature(test.temp, at, bt, ct, dt, et)

# tc = time.perf_counter_ns() 




# d1 = tb - ta
# d2 = tc - tb

