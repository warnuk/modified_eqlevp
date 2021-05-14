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

# define functions

def actp():
    global mol, gact, temp, fi
    global nc, na, nn, ndepact, nzc, nza, ap0, bp0
    global cat, ani, h, b0, b1, b2, c0, c, ta, tc, lc, la, sc, sa, xi
    
    cat = mol[[1, 2, 3, 4, 22, 5, 18, 23, 13]]
    ani = mol[[6, 7, 25, 15, 14, 12, 24, 19, 9, 20, 21, 8]]
    h = mol[[10, 9, 0]]
    
    cat = cat * mh2o / mol[11]
    ani = ani * mh2o / mol[11]
    h = h * mh2o / mol[11]

    if ndepact == 0:
        with open('coefft4', 'r+') as file:
            text = file.read()
            file.close()
            
        lines = text.split('\n')
        values = [line.replace(" ","").split(',') for line in lines]
        
        index = 0
        (nc, na, nn) = (int(v) for v in values[index])
        
        nzc = np.zeros(nc+1)
        nza = np.zeros(na+1)
        (b0, b1, b2, c0) = (np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)), 
                            np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)))
        (sc, sa) = (np.zeros((nc+1, nc+1, na+1)), np.zeros((na+1, na+1, nc+1)))
        (tc, ta) = (np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1)))
        (lc, la) = (np.zeros((nn+1, nc+1)), np.zeros((nn+1, na+1)))
        xi = np.zeros((nn+1, nc+1, na+1))
        
        for i in range(0, nc):
            index += 1
            row = values[index]
            x_S = row[0]
            nzc[i] = int(row[1])
            
        for i in range(0, na):
            index += 1
            row = values[index]
            x_S = row[0]
            nza[i] = int(row[1])
            
        index += 1
        row = values[index]
        x_S = row[0]
        (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
        
        ap0 = temperature(at, bt, ct, dt, et, temp)
        
        for i in range(0, nc):
            for j in range(0, na):
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                b0[i, j] = temperature(at, bt, ct, dt, et, temp)
                
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                b1[i, j] = temperature(at, bt, ct, dt, et, temp)
                
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                b2[i, j] = temperature(at, bt, ct, dt, et, temp)
                
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                c0[i, j] = temperature(at, bt, ct, dt, et, temp)
                
        for i in range(0, nc-1):
            for j in range(i+1, nc):
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                tc[i, j] = temperature(at, bt, ct, dt, et, temp)
                tc[j, i] = tc[i, j]
                
        for i in range(0, na-1):
            for j in range(i+1, na):
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                ta[i, j] = temperature(at, bt, ct, dt, et, temp)
                ta[j, i] = ta[i, j]
                
        for k in range(0, nc-1):
            for i in range(k+1, nc):
                for j in range(0, na):
                    index += 1
                    row = values[index]
                    x_S = row[0]
                    (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                    sc[k, i, j] = temperature(at, bt, ct, dt, et, temp)
                    sc[i, k, j] = sc[k, i, j] 
                    
        for k in range(0, na-1):
            for i in range(k+1, na):
                for j in range(0, nc):
                    index += 1
                    row = values[index]
                    x_S = row[0]
                    (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                    sa[k, i, j] = temperature(at, bt, ct, dt, et, temp)
                    sa[i, k, j] = sa[k, i, j]
                    
        for i in range(0, nn):
            for j in range(0, nc):
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                lc[i, j] = temperature(at, bt, ct, dt, et, temp)
                
        for i in range(0, nn):
            for j in range(0, na):
                index += 1
                row = values[index]
                x_S = row[0]
                (at, bt, ct, dt, et) = (float(v.replace('d', 'e')) for v in row[1:6])
                la[i, j] = temperature(at, bt, ct, dt, et, temp)
                
        for k in range(0, nn):
            for i in range(0, nc):
                for j in range(0, na):
                    xi[k, i, j] = 0               
        xi[1, 8, 0] = -0.0102
        xi[1, 0, 1] = 0.046
        
    bp0 = 1.2
    
    (ec, fc, xc, ea, fa, xa) = (np.zeros((nc+1, nc+1)), np.zeros((nc+1, nc+1)), 
                                np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1)),
                                np.zeros((na+1, na+1)), np.zeros((na+1, na+1)))
    (pp, p, pf, qp, q, qf) = (np.zeros((nc+1, nc+1)), np.zeros((nc+1, nc+1)), 
                              np.zeros((nc+1, nc+1)), np.zeros((na+1, na+1)),
                              np.zeros((na+1, na+1)), np.zeros((na+1, na+1)))
    (cc, bf, b, bp) = (np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)),
                       np.zeros((nc+1, na+1)), np.zeros((nc+1, na+1)))
    (gc, ga, gn) = (np.zeros(nc+1), np.zeros(na+1), np.zeros(nn+1))
    
    ndepact = 1
    u, z = 0., 0.
    for i in range(0, nc):
        u = u + cat[i] * nzc[i] ** 2
        z = z + cat[i] * nzc[i]
    for j in range(0, na):
        u = u + ani[j] * nza[j] ** 2
        z = z + ani[j] * nza[j]
    fi = u / 2.
    fj = fi ** 0.5
    u = 6. * ap0 * fj
    for i in range(0, nc-1):
        for j in range(i+1, nc):
            if (nzc[i] == nzc[j]):
                ec[i, j] = 0.
                fc[i, j] = 0.
            else:
                xc[i, j] = 2. * u
                xc[i, i] = nzc[i] ** 2 * u
                xc[i, j] = nzc[j] ** 2 * u
                
                ec[i, j] = (j0(xc[i, j]) - j0(xc[i, i]) / 2. - 
                            j0(xc[j, j]) / 2.) / fi / 2. 
                fc[i, j] = (xc[i, j] * j1(xc[i, j]) - xc[i, i] * 
                            j1(xc[i, i]) / 2. - xc[j, j] * 
                            j1(xc[j, j]) / 2.) / fi ** 2 / 4 - ec[i, j] / fi
                ec[j, i] = ec[i, j]
                fc[j, i] = fc[i, j]
    for i in range(0, na-1):
        for j in range(i+1, na):
            if (nza[i] == nza[j]):
                ea[i, j] = 0.
                fa[i, j] = 0.
            else:
                xa[i, j] = 2. * u
                xa[i, i] = nza[i] ** 2 * u
                xa[j, j] = nza[j] ** 2 * u
                
                ea[i, j] = (j0(xa[i, j]) - j0(xa[i, i]) / 2. - 
                            j0(xa[j, j]) / 2.) / fi / 2.
                fa[i, j] = (xa[i, j] * j1(xa[i, j]) - xa[i, i] * 
                            j1(xa[i, i]) / 2. - xa[j, j] * j1(xa[j, j]) 
                            / 2.) / fi ** 2 / 4. - ea[i, j] / fi
                ea[j, i] = ea[i, j]
                fa[j, i] = fa[i, j]
    for i in range(0, nc-1):
        for j in range(i+1, nc):
            pp[i, j] = fc[i, j]
            p[i, j] = tc[i, j] + ec[i, j]
            pf[i, j] = p[i, j] + pp[i, j] * fi
            pp[j, i] = pp[i, j]
            p[j, i] = p[i, j]
            pf[j, i] = pf[i, j]
    for i in range(0, na-1):
        for j in range(i+1, na):
            qp[i, j] = fa[i, j]
            q[i, j] = ta[i, j] + ea[i, j]
            qf[i, j] = q[i, j] + qp[i, j] * fi
            qp[j, i] = qp[i, j]
            q[j, i] = q[i, j]
            qf[j, i] = qf[i, j]
    w = fj * 12.
    
    for i in range(0, nc):
        for j in range(0, na):
            cc[i, j] = c0[i, j] / (float(nzc[i] * nza[j])) ** 0.5 / 2.
            if (nzc[i] == 2 and nza[j] == 2):
                v = fj * 1.4
            if (nzc[i] == 1 or nza[j] == 1):
                v = fj * 2.
            bf[i, j] = b0[i, j] + b1[i, j] * np.exp(-v) + b2[i, j] * np.exp(-w)
            b[i, j] = b0[i, j] + b1[i, j] * (g0(v)) + b2[i, j] * (g0(w))
            bp[i, j] = b1[i, j] * (g1(v)) / fi + b2[i, j] * (g1(w)) / fi
    f = -ap0 * (fj / (1. + bp0 * fj) + 2. / bp0 * np.log(1. + bp0 * fj))
    
    for i in range(0, nc):
        for j in range(0, na):
            f = f + cat[i] * ani[j] * bp[i, j]
    for i in range(0, nc-1):
        for j in range(i+1, nc):
            f = f + cat[i] * cat[j] * pp[i, j]
    for i in range(0, na-1):
        for j in range(i+1, na):
            f = f + ani[i] + ani[j] * qp[i, j]
    for ii in range(0, nc):
        u = nzc[ii] ** 2 * f
        for j in range(0, na):
            u = u + ani[j] * (b[ii, j] * 2. + z * cc[ii, j])
        for i in range(0, nc):
            if i != ii:
                v = 0.
                for j in range(0, na):
                    v = v + ani[j] * sc[ii, i, j]
                u = u +cat[i] * (p[ii, i] * 2. + v)
        for i in range(0, na-1):
            for j in range(i+1, na):
                u = u + ani[i] * ani[j] * sa[i, j, ii]
        for i in range(0, nc):
            for j in range(0, na):
                u = u + cat[i] * ani[j] * cc[i, j] * nzc[ii]
        for i in range(0, nn):
            u = u + h[k] * lc[i, ii] * 2.
        for k in range(0, nn):
            for j in range(0, na):
                u = u + h[k] * ani[j] * xi[k, ii, j]
                
        gc[ii] = np.exp(u)
    
    for jj in range(0, na):
        u = nza[jj] ** 2 * f
        for i in range(0, nc):
            u = u + cat[i] * (b[i, jj] * 2. + z * cc[i, jj])
        for i in range(0, na):
            if i != jj:
                v = 0.
                for j in range(0, nc):
                    v = v + cat[j] * sa[jj, i, j]
                u = u + ani[i] * (q[jj, i] * 2. + v)
        
        for i in range(0, nc-1):
            for j in range(i+1, nc):
                u = u + cat[i] * cat[j] * sc[i, j, jj]
        for i in range(0, nc):
            for j in range(0, na):
                u = u + cat[i] * ani[j] * cc[i, j] * nza[jj]
        for j in range(0, nn):
            u = u + h[j] * la[j, jj]
        
        for k in range(0, nn):
            for i in range(0, nc):
                u = u + h[k] * cat[i] * xi[k, i, jj]
        
        ga[jj] = np.exp(u)
    
    for k in range(0, nn):
        u = 0.
        for i in range(0, nc):
            u = u + cat[i] * lc[k, i] * 2.
        for j in range(0, na):
            u = u + ani[j] * la[k, j] * 2.
        for i in range(0, nc):
            for j in range(0, na):
                u = u + cat[i] * ani[j] * xi[k, i, j]
        gn[k] = np.exp(u)
    
    u = -ap0 * fi ** 1.5 / (1. + bp0 * fj)
    for i in range(0, nc):
        for j in range(0, na):
            u = u + cat[i] * ani[j] * (bf[i, j] + z * cc[i, j])
    for i in range(0, nc-1):
        for j in range(i+1, nc):
            v = 0.
            for k in range(0, na):
                v = v + ani[k] * sc[i, j, k]
            u = u + cat[i] * cat[j] * (pf[i, j] + v)
    for i in range(0, na-1):
        for j in range(i+1, na):
            v = 0.
            for k in range(0, nc):
                v = v + cat[k] * sa[i, j, k]
            u = u + ani[i] * ani[j] * (qf[i, j] + v)
    for k in range(0, nn):
        for i in range(0, nc):
            u = u + h[k] * cat[i] * lc[k, i]
    for k in range(0, nn):
        for j in range(0, na):
            u = u + h[k] * ani[j] * la[k, j]
    for k in range(0, nn):
        for i in range(0, nc):
            for j in range(0, na):
                u = u + h[k] * cat[i] * ani[j] * xi[k, i, j]
    s = 0.
    for i in range(0, nc):
        s = s + cat[i]
    for j in range(0, na):
        s = s + ani[j]
    co = 1. + 2. * u / s
    aw = np.exp(-s * co / mh2o)
    gact[0] = gn[2]
    gact[1] = gc[0]
    gact[2] = gc[1]
    gact[3] = gc[2]
    gact[4] = gc[3]
    gact[22] = gc[4]
    gact[5] = gc[5]
    gact[18] = gc[6]
    gact[23] = gc[7]
    gact[13] = gc[8]
    gact[6] = ga[0]
    gact[7] = ga[1]
    gact[25] = ga[2]
    gact[15] = ga[3]
    gact[14] = ga[4]
    gact[12] = ga[5]
    gact[24] = ga[6]
    gact[19] = ga[7]
    gact[20] = ga[8]
    gact[21] = ga[9]
    gact[8] = ga[10]
    gact[10] = aw * aw * gn[0] ** np.log(10.)
    gact[9] = gn[1]
    gact[16] = 1.
    gact[17] = 1.
    gact[11] = aw / mh2o
    ndepact = 1
            
def temperature(at, bt, ct, dt, et, temp):
    temperature = (at + bt * temp + ct * temp ** 2 + dt * temp ** 3 + 
                   et * temp ** 4)
    return(temperature)

def j0(x):
    ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
    j0 = x / (4. + ya * x ** yb * np.exp(yc * x ** yd))
    return(j0)

def j1(x):
    ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
    j1 = (4. + ya * x ** yb * (1. - yb - yc * yd * x ** yd) * 
          np.exp(yc * x ** yd)) / (4. + ya * x ** yb * 
                                   np.exp(yc * x ** yd)) ** 2
    return(j1)

def g0(x):
    g0 = 2. * (1. - (1. + x) * np.exp(-x)) / x ** 2
    return(g0)

def g1(x):
    g1 = -2. * (1. - (1. + x + x ** 2 / 2.) * np.exp(-x)) / x ** 2
    return(g1)


# Start main program


print()
print("This is EVP..............")
print()

print("Starting the evaporation program")

ndepact = 0
ncpt = 0
mwev = 0
fc = 1
q0_S = ""
n = 25
ntot = 12
ncomplex = 14
mh2o = 55.51

tot = np.zeros(ntot)          
totinit = np.zeros(ntot)
tot0 = np.zeros(ntot)
totest = np.zeros(ntot)
psc = np.zeros(ncomplex)

gact = np.zeros(n+1)
gact0 = np.zeros(n+1)
gact1 = np.zeros(n+1)
aq_S = np.chararray(n+1, 10)
atom = np.zeros(n+1)
kmat = np.zeros((n+1,n+1))
nch = np.zeros(n+1)

mol = np.zeros(n+1)
mol0 = np.zeros(n+1)
moll = np.zeros(n+1)
molal = np.zeros(n+1)
molal0 = np.zeros(n+1)
act = np.zeros(n+1)
act0 = np.zeros(n+1)


aquv = pd.read_csv('aquv.dat', header=None)
aquv.iloc[:,0] = [i.lstrip() for i in aquv.iloc[:,0].values]

aq_S = aquv.iloc[:,0].values
atom = aquv.iloc[:,1].values
nch = aquv.iloc[:,2].values

kmat = pd.read_csv('matrice2', header=None).values

with open('stockage') as stockage:
    transf_S = (stockage.read().lstrip()).strip('\n')
    with open(transf_S) as trace:
        raw_trace = trace.read()
        transf_S = [i.lstrip() for i in raw_trace.split('\n')]
        trace.close()
    stockage.close()

index = 0 
temp = float(transf_S[index])

index += 1
tinit = float(transf_S[index])

index += 1
ph = float(transf_S[index])

index += 1
phinit = float(transf_S[index])

index += 1
po = float(transf_S[index])

index += 1
poinit = float(transf_S[index])

index += 1
diltot = float(transf_S[index])

index += 1
constituant_S = transf_S[index].split(',')

index += 1
totinit[0:10] = [float(transf_S[i]) for i in (np.arange(0, 10)+index).tolist()]
index += 9

nbmin = np.count_nonzero(totinit)
ica = np.zeros(n+nbmin+1)
kinvar = np.zeros(nbmin+1+3)

index += 1
mol[0:n+1] = [float(transf_S[i]) for i in (np.arange(0, n+1)+index).tolist()]
mol0[:] = mol[:]
index += n

index += 1
syst_S = transf_S[index]

index +=1
xi = float(transf_S[index])

index +=1
npasecran = int(transf_S[index])

index +=1
prnt_S = transf_S[index]

index +=1
npasimp = int(transf_S[index])

index +=1
stockfich_S = transf_S[index]

index += 1
npasfich = int(transf_S[index])

index += 1
unit_S = transf_S[index]

index += 1
fich_S = transf_S[index]

index += 1
fichmin_S = transf_S[index]

index += 1
molemin_S = transf_S[index]

index += 1
miner_S = transf_S[index]

index += 1
stdmax = float(transf_S[index])

index += 1
pkmol = float(transf_S[index])

index += 1
pkeq = float(transf_S[index])


if npasimp > 0:
    logfile = open(prnt_S, 'a+')

if xi == 0:
    inc_S = "auto"
else:
    inc_S = "manu"

xi0 = xi

ica[1:n+1] = 1


if stockfich_S == "y":
    with open(fich_S, "w+") as file:
        file.write(",".join(constituant_S))
        file.close()

    with open(fichmin_S, "w+") as file:
        text = "Temperature of solution = {tinit}" \
                "   Deg C.     Temperature of simulation = {temp} " \
                "   Deg C.\n".format(tinit=tinit, temp=temp)
        file.write(text)
        
        if diltot > 1:
            text = "The initial solution has been diluted " \
                "{diltot} times\n".format(diltot=diltot)
            file.write(text)
        
        if ph != phinit:
            text = "Initial Log(pco2) = {poinit}    Selected Log(pco2) = " \
                "{po}\n".format(poinit=poinit, po=po)
            file.write(text)
            
            text = "Initial ph = {phinit}     Calculated ph = " \
                "{ph}\n".format(phinit=phinit, ph=ph)
            file.write(text)
            
        file.close()
        
    with open(molemin_S, "w+") as file:
        file.close()
    
with open("complex3", "r+") as file:
    text = file.read().split('\n')
    file.close()


for i in range(0, ncomplex):
    a_S, at, bt, ct, dt, et = text[i].split(',')
    [at, bt, ct, dt, et] = [float(i) for i in [at, bt, ct, dt, et]]
    
    psc[i] = 10 ** (at + bt / 300. * temp + ct/30000. * temp ** 2 + dt /
                    3000000. * temp ** 3 + et / 300000000. * temp ** 4)
psc2 = psc[2]
psc[2] = psc[2] * psc[13] * 10 ** po

with open(miner_S, 'r+') as file:
    text = file.read()
    file.close()

lines = text.split('\n')
values = [line.replace(" ","").split(',') for line in lines]

index = 0
(nc, na, nm) = (int(i) for i in values[index])

ncm = nc + na + 1

wmin = np.zeros((nm, ncm+1))

mu = {}

for i in range(0, ncm):
    index += 1
    a_S, at, bt, ct, dt, et = values[index]
    (at, bt, ct, dt, et) = (float(x) for x in (at, bt, ct, dt, et))
    a_S = a_S.lower()
    
    j = np.where(aq_S == a_S)[0][0]
    mu[j] = (at + bt/300. * temp + ct / 30000. * temp ** 2 +
            dt / 3000000. * temp ** 3 + et / 300000000. * temp ** 4)

mineral_S = {}
mum = {}
for k in range(0, nm):
    index += 1
    
    mineral_S[k], ncomp = values[index][0:2]
    ncomp = int(ncomp)
    
    c_ion, nom_ion_S = {}, {}
    for component in range(0, ncomp):
        idx = 2 * (component + 1)
        c_ion[component], nom_ion_S[component] = values[index][idx:idx+2]
        c_ion[component] = float(c_ion[component])
    
    (at, bt, ct, dt, et) = (float(x) for x in values[index][-5:])
    
    for i in range(0, ncomp):
        x_S = nom_ion_S[i].lower()
        
        for j in range(0, ncm+1):
            if x_S == aq_S[j]:
                wmin[k, j] = c_ion[i]
        
    mum[k] = (at + bt / 300. * temp + ct / 30000. * temp ** 2 + 
              dt / 3000000. * temp ** 3 + et / 300000000. * temp ** 4)

psol, psol0 = {}, {}

for k in range(0, nm):
    u = mum[k]
    for i in range(0, ncm):
        u = u - wmin[k, i] * mu[i+1]
    psol[k] = np.exp(u)
    psol0[k] = psol[k]
    
if stockfich_S == "y":
    with open(molemin_S, "w+") as file:
        text = ["fc"]
        for mineral in mineral_S.values():
            text.append(mineral)
        file.write(",".join(text))

sc, cmax, sa, amax = 0., 0., 0., 0.
for i in range(0, n):
    if nch[i] > 0:
        sc = sc + mol[i] * nch[i]
        if (mol[i] * nch[i]) > cmax:
            cmax = mol[i] * nch[i]
            ic = i
    elif nch[i] < 0:
        sa = sa + mol[i] * (-nch[i])
        if (mol[i] * (-nch[i])) > amax:
            amax = mol[i] * (-nch[i])
            ia = i

dca = 200. * abs(sc - sa) / (sc + sa)
delta = sc - sa
mol[ic] = mol[ic] - delta / 2. / nch[ic]
mol[ia] = mol[ia] + delta / 2. / (-nch[ia])
sc, sa = 0., 0.
for i in range(0, n):
    if nch[i] > 0:
        sc = sc + mol[i] * nch[i]
    elif nch[i] < 0:
        sa = sa + mol[i] * (-nch[i])
print("Sum of cations = {sc}".format(sc=sc))
print("Sum of anions = {sa}".format(sa=sa))
print()

for i in range(0, 11):
    totinit[i] = 0
    
    for j in range(0, n):
        totinit[i] = totinit[i] + kmat[i, j] * mol[j]
        tot[i] = totinit[i]
        tot0[i] = totinit[i]
        
tot[11] = 0.
ctot0 = mol[0] + mol[14] + mol[15] + mol[16] + mol[17]

""" STOPPING AT LINE 321 OF MAIN_EVP.F90, 4:30PM ON 5/9/2021 """

actp()

for i in range(0, n+1):
    molal[i] = mol[i] * mh2o / mol[11]
    act[i] = molal[i] * gact[i]
for k in range(0, nm):
    pai[k] = 1.
    for i in range(0, ncm):
        pai[k] = pai[k] * act[i] ** wmin[k, i]
    if pai[k] >= psol[k]:
        lmin[k] = 1
if stockfich_S == "y":
    with open(fichmin_S, 'w+') as file:
        for k in range(0, nm):
            if lmin[k] == 1:
                text = "Initial solution oversaturated in " \
                    "{mineral}\n".format(mineral=mineral_S[k])
                file.write(text)
        file.close()

for k in range(0, nm):
    if lmin[k] == 1:
        psol[k] = pai[k]



# 500 CONTINUE !debut:
ncmpt = 0
ix, iy = 1, 2
initdeseq = 0
if mwev == 0.:
    for k in range(0, nm):
        if (psol[k] * 0.95) > psol0[k]:
            psol[k] = psol[k] * 0.95
            print("{a}     {psol}     {psol0}".format(a=mineral_S[k], 
                                                      psol=psol[k], 
                                                      psol0=psol0[k]))
        elif (psol[k] * 0.95) <= psol[k]:
            psol[k] = psol0[k]

nw = 1
while nw != 0:
    
    ncmpt += 1
    m0_S = ""
    for k in range(0, nm):
        if lmin[k] == 1:
            m0_S = "_".join(m0_S, mineral_S[k])
    for i in range(0, n+1):
        gact1[i] = gact[i]
    
    actp()
    
    if kinvariant == 0:
        for i in range(0, n+1):
            gact[i] = (gact[i] + gact1[i] * ix) / iy
    
    for i in range(0, n+1):
        molal[i] = mol[i] * mh2o / mol[11]
        act[i] = molal[i] * gact[i]
    
    for k in range(0, nm):
        pai[k] = 1.
        for i in range(0, ncm):
            pai[k] = pai[k] * act[i] ** wmin[k, i]
        if pai[k] >= psol[k]:
            if min[k] >= 0.:
                lmin[k] = 1
            elif min[k] < 0.:
                lmin[k] = 0
                min[k] = 0.
        elif pai[k] < psol[k]:
            if min[k] <= 0:
                lmin[k] = 0
                min[k] = 0.
            elif min[k] > 0.:
                lmin[k] = 1
    
    for k in range(0, nm):
        if psol[k] == 1:
            if pai[k] < psol0[k]*0.9:
                linvar[k] = 0
            elif pai[k] >= psol0[k]:
                linvar[k] = 1
    mineraux_S = []
    nminer = 0
    
    """ LINE 427 """
    for k in range(0, nm):
        if lmin[k] == 1:
            nminer += 1
            mineraux_S.append(mineral_S[k])
    if (ncpt == 1 or (ncpt%npasecran) == 0):
        if nminer == 0:
            print(ncmpt, "No_minerals", sep=" ")
        else:
            print(ncmpt, "_".join(mineraux_S))
    
    if (mwev > 0. and fc != 1. and nminer - nminer0 >= 2):
        xi = xi / 2
        if xi < epsilon:
            text = "Program unstable\nRestart the initialization" \
                "program (EQL...)\nand lower the limits of convergence\n"
            print(text)
            if stockfich_S == "y":
                with open(fichmin_S, 'w+') as file:
                    file.write(text)
                    file.close()
""" LINE 455 """


# Compact function
def compact(lmin):
    global lmin, fc, q_min
    
    q_min = np.empty(nm+1, dtype=np.object)
    
    print("Compacting mineral file")
    
    for i in range(0, nm):
        lmin[i] = 0
    
    nbmin_comp = 0
    
    with open(molemin_S, 'r+') as file:
        text = file.read()
        file.close()
    text = text.split(',')
    for i in range(1, len(text)):
        q_min[i] = text[i]


# densite function
def densite():
    global dens
    
    ncdens, nadens = 5, 5
    
    s = np.zeros((ncdens, nadens))
    cat, ani = np.zeros(ncdens), np.zeros(nadens)
    ic, ia = np.zeros(ncdens), np.zeros(nadens)
    
    for i in range(0, 8):
        if i+1 <= ncdens:
            cat[i] = mol[i+1] / mol[11] * mh2o
        if i+1 > ncdens:
            ani[i-ncdens] = mol[i+1] / mol[11] * mh2o
    
    ani[3] = mol[15] / mol[11] * mh2o
    ani[4] = mol[14] / mol[11] * mh2o
    
    for i in range(0, ncdens):
        ic[i] = nch[i+1]
    for i in range(0, 3):
        ia[i] = -nch[i+6]
    
    if ncpt == 1:
        au = np.zeros((ncdens, nadens))
        bu = np.zeros((ncdens, nadens))
        
        with open("densite", 'r+') as file:
            text = file.read()
            file.close()
        values = [line.split(',') for line in text.split('\n')]
        
        index = 0
        for i in range(0, 5):
            for j in range(0, 5):
                au[i, j] = float(values[index][3])
                bu[i, j] = float(values[index][4])
                index += 1
        
    dens = 1.
    u = 0.
    
    for j in range(0, nadens):
        u = u + ani[j] * ia[j]
    for i in range(0, ncdens):
        for j in range(0, nadens):
            s[i, j] = int((ic[i] + ia[j]) / 2) * cat[i] * ani[j] / u
            dens = dens + au[i, j] * s[i, j] + bu[i, j] * s[i, j] ** 2


                
# Stop simulation function
