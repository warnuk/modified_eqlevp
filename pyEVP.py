# -*- coding: utf-8 -*-
"""
Created on Sun May  9 13:10:17 2021

@author: warnu
"""

"""
Simulation of evaporation
System: Na-K-Li-Ca-Mg-Cl-SO4-NO3-CO2-HCO3-CO3-H4SiO4-B(OH)3-H-OH-H2O
Ionic interaction (Pitzer)
Solution of the system of equations by Newton-Raphson method
Temperature dependence included for minerals
  and for some Pitzer's coefficients (0-50C)
Connected files: MATRICE2-MURTF3/MURTF0-COMPLEX3-COEFFT4-DENSITE
Program chained to EQL through transfer file: *.tra
Last revision: march 1999
author: F. RISACHER
translated in FORTRAN 90 by A. CLEMENT
translated in Python 3 by W. ARNUK
"""

import pandas as pd
import numpy as np
from time import perf_counter_ns
from geochemistry import actp, evp_actp
#import pyEQL

def read_file(file):
    with open(file, 'r+') as file:
        text = list(filter(None, file.read().split('\n')))
        file.close()
        values = [line.split(',') for line in text]
    return(values)

# EVP local variables
ncpt = 0
mwev = 0
fc = 1
q0_S = ""
n=25
ntot=12
ncomplex = 14
mh2o=55.51

tot = np.zeros((ntot+1))
totinit = np.zeros((ntot+1))
tot0 = np.zeros((ntot+1))
totest = np.zeros((ntot+1))
psc = np.zeros(ncomplex+1)

gact = np.zeros((n+1))
gact0 = np.zeros((n+1))
gact1 = np.zeros((n+1))
aq_S = np.empty((n+1), np.object)
atom = np.zeros((n+1))
kmat = np.zeros((n+1, n+1))
nch = np.zeros((n+1))

mol = np.zeros((n+1))
mol0 = np.zeros((n+1))
mol1 = np.zeros((n+1))
molal = np.zeros((n+1))
molal0 = np.zeros((n+1))
act = np.zeros((n+1))
act0 = np.zeros((n+1))

kmat = np.zeros((n+1, n+1))
ica = np.zeros(n+1)

# Read aquv data file
text = read_file("aquv.dat")
aq_S[:] = [line[0].lstrip() for line in text]
atom[:] = [float(line[1]) for line in text]
nch[:] = [float(line[2]) for line in text]


# Read k matrix
text = read_file("matrice2")
kmat[1:,1:] = text

# test case for writing code, comes from pyEQL test Water
with open("stockage", "r") as file:
    transfer = read_file(file.read())
    file.close()

temp = float(transfer[0][0])
tinit = float(transfer[1][0])
ph = float(transfer[2][0])
phinit = float(transfer[3][0])
po = float(transfer[4][0])
poinit = float(transfer[5][0])
diltot = float(transfer[6][0])
constit_S = transfer[7]

totinit[1:11] = [float(i[0]) for i in transfer[8:18]]
nbmin = np.count_nonzero(totinit)
ica = np.zeros((n+nbmin+1))
kinvar = np.zeros((nbmin+4))

mol[0:n+1] = [float(i[0]) for i in transfer[18:18+n+1]]
mol0[:] = mol[:]

syst_S = transfer[18+n+1][0]
incr = float(transfer[18+n+2][0])
print_step = float(transfer[18+n+3][0])
p_S = transfer[18+n+4][0]
output_step = float(transfer[18+n+5][0])
storage_step = float(transfer[18+n+6][0])
unit_S = transfer[18+n+7][0]
chem_file = transfer[18+n+8][0]
event_file = transfer[18+n+9][0]
min_file = transfer[18+n+10][0]
miner_S = transfer[18+n+11][0]
stdmax = float(transfer[18+n+12][0])
pkmol = float(transfer[18+n+13][0])
pkeq = float(transfer[18+n+14][0])




# if print_step > 0:
#     logfile = open("{}.log".format(label), "w+")

if incr == 0:
    inc_S = "auto"
else:
    inc_S = "manu"
incr0 = incr

ica[1:n+1] = 1

# Write the constituent headers to the chemistry file
with open(chem_file, "w") as file:
    file.write(",".join(constit_S))
    file.write('\n')
    file.close()
    
# Write the first mineral events 
with open(event_file, 'w') as file:
    file.write("Temperature of solution = {} Deg C".format(tinit))
    file.write("     ")
    file.write("Temperature of simulation = {} Deg C".format(temp))
    file.write("\n")
    if diltot > 1:
        file.write("The initial solution has "\
                          "been diluted {} times".format(diltot))
        file.write("\n")
    if ph != phinit:
        file.write("Initial Log(pco2) = {}     ".format(poinit))
        file.write("Selected Log(pco2) = {}".format(po))
        file.write("\n")
        file.write("Initial pH = {}     ".format(phinit))
        file.write("Calculated pH = {}".format(ph))
        file.write("\n")
    file.close()


with open(min_file, "w") as file:
    file.close()


text = read_file("complex3")
complex3 = np.zeros((ncomplex, 5))
complex3[:, :] = [i[1:] for i in text]

psc[1:] = 10 ** (complex3[:,0] + 
                 complex3[:,1] / 300 * temp + 
                 complex3[:,2] / 30000 * temp ** 2 + 
                 complex3[:,3] / 3000000 * temp ** 3 + 
                 complex3[:,4] / 300000000 * temp ** 4)

psc3 = psc[3]
psc[3] = psc[3] * psc[14] * 10 ** po


text = read_file(miner_S)
(nc, na, nm) = (int(i) for i in text[0])
ncm = nc + na + 1

wmin = np.zeros((nm+1, ncm+1))
mu = np.zeros((ncm+1))
linvar = np.zeros((nm+1))

mineral_S = np.empty((nm+1), np.object)
mum = np.zeros((nm+1))
psol = np.zeros((nm+1))
psol0 = np.zeros((nm+1))
pai = np.zeros((nm+1))
pai0 = np.zeros((nm+1))

lmin = np.zeros((nm+1))
lmin0 = np.zeros((nm+1))
lmin1 = np.zeros((nm+1))

minerals = np.zeros((nm+1))
minerals0 = np.zeros((nm+1))
minp = np.zeros((nm+1))
minp0 = np.zeros((nm+1))

# nom_ion_S = np.empty((ncm+1), np.object)
# c_ion = np.zeros((ncm+1))

#species_df = pd.DataFrame(text[1:ncm+1]).values[:,1:].astype(np.float64)

ion_list = text[1:ncm+1]
ion_list.insert(0, None)
#ion_S[:] = [i[0].lower() for i in ion_list]
#ion_db = np.zeros((ncm, 5))
#ion_db[:,:] = [i[1:] for i in ion_list]

for i in range(1, ncm+1):
    a_S = str(ion_list[i][0]).lower()
    (at, bt, ct, dt, et) = (float(x) for x in ion_list[i][1:])
    for j in range(1, ncm+1):
        if a_S == str(aq_S[j]):
            mu[j] = (at + 
                     bt / 300 * temp + 
                     ct / 30000 * temp ** 2 + 
                     dt / 3000000 * temp ** 3 + 
                     et / 300000000 * temp ** 4)


# mu[:] = (ion_db[:,0] + 
#          ion_db[:,1] / 300 * temp + 
#          ion_db[:,2] / 30000 * temp ** 2 + 
#          ion_db[:,3] / 3000000 * temp ** 3 + 
#          ion_db[:,4] / 300000000 * temp ** 4)

mineral_list = text[ncm+1:]
mineral_list.insert(0, None)
mineral_S = np.empty(nm+1, np.object)

for k in range(1, nm+1):
    line = mineral_list[k]
    mineral_S[k] = line[0]
    ncomp = int(line[1])
    c_ion = np.zeros(ncomp+1)
    nom_ion_S = np.empty(ncomp+1, dtype=np.object)
    for i in range(1, ncomp+1):
        c_ion[i] = float(line[i*2])
        nom_ion_S[i] = str(line[1+i*2])
        for j in range(0, ncm+1):
            x_S = nom_ion_S[i].lower()
            if x_S == aq_S[j]:
                wmin[k, j] = c_ion[i]
    (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
    mum[k] = (at + 
              bt / 300 * temp + 
              ct / 30000 * temp ** 2 + 
              dt / 3000000 * temp ** 3 + 
              et / 300000000 * temp ** 4)


# mineral_S = np.empty(nm+1, np.object)
# min_db = np.zeros((nm+1, 5))

# for k in range(1, nm+1):
#     line = text[ncm+k]
#     mineral_S[k] = line[0]
#     ncomp = int(line[1])
#     c_ion = np.zeros(ncomp+1)
#     nom_ion_S = np.empty(ncomp+1, dtype=np.object)
#     for i in range(1, ncomp+1):
#         c_ion[i] = float(line[i*2])
#         nom_ion_S[i] = str(line[1+i*2])
#         for j in range(0, ncm+1):
#             x_S = nom_ion_S[i].lower()
#             if x_S == aq_S[j]:
#                 wmin[k, j] = c_ion[i]
#     (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
#     min_db[k] = [at, bt, ct, dt, et]

# mum[:] = (min_db[:,0] + 
#           min_db[:,1] / 300 * temp + 
#           min_db[:,2] / 30000 * temp ** 2 + 
#           min_db[:,3] / 3000000 * temp ** 3 + 
#           min_db[:,4] / 300000000 * temp ** 4)


for k in range(1, nm+1):
    u = mum[k]
    for i in range(0, ncm+1):
        u = u - wmin[k, i] * mu[i]
    psol[k] = np.exp(u)
    psol0[k] = psol[k]
#psol[1:] = np.exp(mum - np.sum((wmin[:,:] * mu), axis=1))[1:]
#psol0[:] = psol[:]

with open(min_file, "w") as file:
    file.write("fc,")
    file.write(",".join(mineral_S[1:].tolist()))
    file.write('\n')
    file.close()

# cations = nch > 0
# anions = nch < 0
# mch_prod = mol * nch
# sc = np.sum(mch_prod[cations])
# cmax = np.max(mch_prod)
# ic = np.argmax(mch_prod)
# sa = -np.sum(mch_prod[anions])
# amax = -np.min(mch_prod)
# ia = np.argmin(mch_prod)

sc, cmax = 0, 0
for i in range(1, n+1):
    if nch[i] > 0:
        sc += mol[i] * nch[i]
        if mol[i] * nch[i] > cmax:
            cmax = mol[i] * nch[i]
            ic = i
sa, amax = 0, 0
for i in range(1, n+1):
    if nch[i] < 0:
        sa += mol[i] * -nch[i]
        if mol[i] * -nch[i] > amax:
            amax = mol[i] * -nch[i]
            ia = i

dca = 200 * np.abs(sc-sa) / (sc+sa)
delta = sc-sa
mol[ic] = mol[ic] - delta / 2 / nch[ic]
mol[ia] = mol[ia] + delta / 2 / -nch[ia]

# mch_prod = mol * nch
# sc = np.sum(mch_prod[cations])
# sa = -np.sum(mch_prod[anions])

sc, sa = 0, 0
for i in range(1, n+1):
    if nch[i] > 0:
        sc += mol[i] * nch[i]
    if nch[i] < 0:
        sa += mol[i] * -nch[i]

print("sum of cations = {}".format(sc))
print("sum of anions = {}".format(sa))
print()

for i in range(1, 12):
    totinit[i] = 0
    for j in range(1, n+1):
        totinit[i] += kmat[i, j] * mol[j]
        tot[i] = totinit[i]
        tot0[i] = totinit[i]

tot[12] = 0
ctot0 = mol[0] + mol[14] + mol[15] + mol[16] + mol[17]

(gact, nc, na, nn, nzc, nza, ndepact, 
 ap0, bp0, mh2o, b0, b1, b2, c0, tc, 
 ta, lc, la, sc, sa, xi, fi, aw) = evp_actp(mol, temp)

for i in range(0, n+1):
    molal[i] = mol[i] * mh2o / mol[11]
    act[i] = molal[i] * gact[i]
    
for k in range(1, nm+1):
    pai[k] = 1
    for i in range(1, ncm+1):
        pai[k] = pai[k] * act[i] ** wmin[k, i]
    if pai[k] >= psol[k]:
        lmin[k] = 1

with open(event_file, "a") as file:
    for k in range(1, nm+1):
        if lmin[k] == 1:
            file.write("Initial solution "\
                       "oversaturated in {}".format(mineral_S[k]))
            file.write("\n")
    file.close()
