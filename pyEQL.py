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


kmat = np.zeros((n, n))
ica = np.zeros(n+1)

aqu = pd.read_csv("aqu.dat", header=None)
aq_S = aqu.iloc[:,0].values
atom = aqu.iloc[:,1].values
nch = aqu.iloc[:,2].values

kmat = pd.read_csv("matrice1", header=None).values
kmat[13, 0], kmat[15, 0], kmat[20, 0], kmat[21, 0] = -1, -1, 3, 5

complex3 = pd.read_csv("complex3", header=None)
complex3.columns = ['species', 'at', 'bt', 'ct', 'dt', 'et']

with open("murtf2", 'r+') as file:
    text = [i.split(',') for i in file.read().split('\n')][:-1]
    file.close()

(nc, na, nm) = (int(i) for i in text[0])
nt = nc + na

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


mineral_db = text[nt+2:]
mineral_S = [i[0] for i in mineral_db]

for k in range(0, nm):
    ncomp = int(mineral_db[k][1])
    c_ion = np.zeros(ncomp)
    nom_ion_S = np.empty(shape=ncomp, dtype=np.object)
    for i in range(0, ncomp):
        c_ion[i] = float(mineral_db[k][2+2*i])
        nom_ion_S[i] = mineral_db[k][3+2*i]
        for j in range(0, nt+1):
            if nom_ion_S[i].lower() == aq_S[j]:
                wmin[k, j] = c_ion[i]

# LINE 214



#x = ion_db['species'].map(dict(zip(aq_S.tolist(), range(aq_S.shape[0]))))


# Functions

filepath = 'test.dat'

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
                 ca, mg, cl, so4, no3, b, si, alk):
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
        
        self.tot = np.zeros(12)
        self.tot[0] = self.na / 1000
        self.tot[1] = self.k / 1000
        self.tot[2] = self.li / 1000
        self.tot[3] = self.ca / 1000
        self.tot[4] = self.mg / 1000
        self.tot[5] = self.cl / 1000
        self.tot[6] = self.so4 / 1000
        self.tot[7] = self.no3 / 1000
        self.tot[8] = self.b / 1000
        self.tot[9] = self.si / 1000
        self.tot[10] = 10 ** -ph
        self.tot[11] = self.alk / 1000
        
        
        self.psc = 10 ** (complex3['at'] + 
                          complex3['bt'] / 300 * self.temp + 
                          complex3['ct'] / 30000 * self.temp ** 2 + 
                          complex3['dt'] / 3000000 * self.temp ** 3 + 
                          complex3['et'] / 300000000 * self.temp ** 4)
        
        self.mu = 10 ** (ion_db['at'] + 
                         ion_db['bt'] / 300 * self.temp + 
                         ion_db['ct'] / 30000 * self.temp ** 2 + 
                         ion_db['dt'] / 3000000 * self.temp ** 3 + 
                         ion_db['et'] / 300000000 * self.temp ** 4)
        
test = input_data(filepath)

