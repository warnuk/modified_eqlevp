# -*- coding: utf-8 -*-
"""
Created on Mon May 24 10:43:01 2021

@author: warnu
"""

import numpy as np
from datetime import datetime

def read_file(file):
    with open(file, 'r+') as file:
        text = list(filter(None, file.read().split('\n')))
        file.close()
        values = [line.split(',') for line in text]
    return(values)

def eql(Water, pCO2, system, units="molal", dilute=1, 
        add_minerals=None, rem_minerals=None, max_salinity=0, 
        verbose=True, output=True, call_evp=True):
    
    
    if verbose:
        print("\nThis is EQL..............\n")
        print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')
        
    
    # Set local variables for simulation
    n = 25
    ntot = 12
    ncomplex = 14
    n0 = 10
    max_cat = 5
    max_an = 5
    epsilon = 1e-8
    nchcat, nchani = 0, 0
    max_constit = 30
    max_nr = 10
    max_conv=100
    
    pk = 0.1
    pk0 = pk
    pkf = 0.0001
    pkstd = 0.1
    dph = 0.2
    diltot = 1
    mh2o = 55.51
    nconv = 2
    eps = 1e-12
    
    # Initialize blank arrays
    constit_S = np.empty(max_constit+1, np.object)
    psc = np.zeros(15)
    tot = np.zeros(ntot+1)
    tot0 = np.zeros(ntot+1)
    totinit = np.zeros(n0+1)
    nch = np.zeros(n+1)
    molal = np.zeros(n+1)
    act = np.zeros(n+1)
    gact = np.zeros(n+1)
    aq_S = np.empty(n+1, np.object)
    atom = np.zeros(n+1)
    kmat = np.zeros((n+1, n+1))
    ica = np.zeros(n+1)
    xx = np.zeros(n+1)
    z = np.zeros((n+1, n+1))
    zz = np.zeros(n+1)
    cat = np.zeros(max_cat+1)
    ani = np.zeros(max_an+1)
    nchcat = np.zeros(max_cat+1)
    nchani = np.zeros(max_an+1)
    
    # Assign chemistry from Water object to local variables
    label = Water.label
    temp = Water.temp
    dens = Water.dens
    pH = Water.ph
    
    
    
    
    tot[1] = Water.na / 1000
    tot[2] = Water.k / 1000
    tot[3] = Water.li / 1000
    tot[4] = Water.ca / 1000
    tot[5] = Water.mg / 1000
    tot[6] = Water.cl / 1000
    tot[7] = Water.so4 / 1000
    tot[8] = Water.no3 / 1000
    tot[9] = Water.b / 1000
    tot[10] = Water.si / 1000
    tot[11] = 10 ** -pH
    tot[12] = Water.alk / 1000
    
    tinit = temp
    
    # Read aqu database
    aqu = read_file("aqu.dat")
    aq_S[:] = [line[0] for line in aqu]
    atom[:] = [line[1] for line in aqu]
    nch[:] = [line[2] for line in aqu]
    
    # Read kmat from matrice1
    kmat[1:,1:] = [line[:] for line in read_file("matrice1")]
    kmat[13, 0] = -1
    kmat[15, 0] = -1
    kmat[20, 0] = 3
    kmat[21, 0] = 5
    
    # Read thermodynamic data for dissociation coefficients from complex3
    complex3 = np.zeros((ncomplex+1, 5))
    complex3[1:,:] = [i[1:] for i in read_file("complex3")]
    psc[:] = 10 ** (complex3[:,0] + 
                    complex3[:,1] / 300 * temp + 
                    complex3[:,2] / 30000 * temp ** 2 + 
                    complex3[:,3] / 3000000 * temp ** 3 + 
                    complex3[:,4] / 300000000 * temp ** 4)
    
    # Read eql mineral database from murtf2
    murtf2 = read_file("murtf2")
    (nc, na, nm) = (int(i) for i in murtf2[0])
    nt = nc + na
    
    wmin = np.zeros((nm+1, nt+1))
    lmin = np.zeros(nm+1)
    nwmin = np.zeros(nm+1)
    mineral_S = np.empty(nm+1, np.object)
    mineral0_S = np.empty(nm+1, np.object)
    mu = np.zeros(nt+1)
    mum = np.zeros(nm+1)
    psol = np.zeros(nm+1)
    pai = np.zeros(nm+1)
    kinvar = np.zeros(nt+1)
    
    ion_list = murtf2[1:nt+2]
    ion_list.insert(0, None)
    
    for i in range(1, nt+2):
        a_S = str(ion_list[i][0]).lower()
        (at, bt, ct, dt, et) = (float(x) for x in ion_list[i][1:])
        for j in range(1, nt+1):
            if a_S == str(aq_S[j]):
                mu[j] = (at + 
                         bt / 300 * temp + 
                         ct / 30000 * temp ** 2 + 
                         dt / 3000000 * temp ** 3 + 
                         et / 300000000 * temp ** 4)
    
    mineral_list = murtf2[nt+2:]
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
            for j in range(0, nt+1):
                x_S = nom_ion_S[i].lower()
                if x_S == aq_S[j]:
                    wmin[k, j] = c_ion[i]
                    
        (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
        mum[k] = (at + 
                  bt / 300 * temp + 
                  ct / 30000 * temp ** 2 + 
                  dt / 3000000 * temp ** 3 + 
                  et / 300000000 * temp ** 4)
        
    for k in range(1, nm+1):
        u = mum[k]
        for i in range(0, nt+1):
            u -= wmin[k, i] * mu[i]
        psol[k] = np.exp(u)
    
    sc, cmax = 0, 0
    for i in range(1, ntot+1):
        if nch[i] > 0:
            sc += tot[i] * nch[i]
            if tot[i] * nch[i] > cmax:
                cmax = tot[i] * nch[i]
                icat = i
    sa, amax = 0, 0
    for i in range(1, ntot+1):
        if nch[i] < 0:
            sa += tot[i] * -nch[i]
            if tot[i] * -nch[i] > amax:
                amax = tot[i] * -nch[i]
                iani = i
    
    if sc + sa != 0:
        dca = 200 * np.abs(sc-sa) / (sc+sa)
    else:
        dca = 0
    delta = sc-sa
    
    if verbose:
        print("sum of cations = {}".format(sc))
        print("sum of anions = {}".format(sa))
        print("Electrical balance = {} %".format((dca * 100 + 0.5) / 100))
        
    tot[icat] = tot[icat] - delta / 2 / nch[icat]
    tot[iani] = tot[iani] + delta / 2 / -nch[iani]
    
    tot0[1:13] = tot[1:13]
    
    
    # Set the default mineral database to murtf3
    min_S = "murtf3"
    murtf3 = read_file(min_S)
    
    nc, na, nm0 = [int(i) for i in murtf3[0]]
    
    mineral0_S[:] = [i[0] for i in murtf3[(2+nc+na):]]
    nwmin[np.in1d(mineral_S, mineral0_S)] = 1
    
    (molal, stdi, ee) = calculate_molalities(molal, units, tot, psc, n, 
                                             ntot, nch, atom, dens, pk, aq_S,
                                             verbose)
    
    
def evp(filename, temp, init, ph, phinit, po, poinit, diltot, constit_S,
        totinit, mol, system, increment, units="molal", mineral_db="murtf3", 
        max_salinity=0, pkmol=0.001, pkeq=.0000000000001, verbose=True, 
        output=True, output_step=1, log=True, log_step=1):
    
    pass

def transfer():
    pass

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

def calculate_molalities(molal, units, tot, psc, n, ntot, nch, atom, 
                         dens, pk, aq_S, verbose=True):
    
    molal[1] = tot[1]
    molal[2] = tot[2]
    molal[3] = tot[3]
    molal[6] = tot[6]
    molal[8] = tot[8]
    molal[11] = tot[11]
    molal[13] = psc[1] / molal[11]
    if tot[9] > 0:
        a = 1 + molal[13] / psc[7]
        b = 3 * psc[8] * molal[13]
        c = 4 * psc[9] * molal[13] ** 2
        xu = tot[9] / 2
        u = xu
        
        eq = a * xu + b * xu ** 3 + c * xu ** 4
        while (200 * abs(eq - tot[9]) / 
                       (eq + tot[9]) >= pk):
            
            u = u / 2
            if eq > tot[9]:
                xu -= u
            else:
                xu += u
            
            eq = a * xu + b * xu ** 3 + c * xu ** 4
        
        molal[9] = xu
        molal[19] = (molal[9] * molal[13] / psc[7])
        molal[20] = (molal[13] * molal[9] ** 3 * 
                          psc[8])
        molal[21] = (molal[13] ** 2 * molal[9] ** 4 * 
                          psc[9])
    molal[14] = (tot[12] + molal[11] - molal[13] - 
                      molal[19] - molal[20] - 
                      2 * molal[21]) / (2 + molal[11] / 
                                             psc[2])
    molal[12] = (tot[12] + molal[11] - molal[13] - 
                      molal[19] - molal[20] - 
                      2 * molal[21]) / (1 + 2 * psc[2] / 
                                             molal[11])
                                             
    molal[15] = molal[12] * molal[11] / psc[3]
    molal[4] = tot[4] / (1 + molal[14] / psc[4] +
                                   molal[19] / psc[10])
    molal[16] = molal[4] * molal[14] / psc[4]
    molal[22] = molal[4] * molal[19] / psc[10]
    molal[5] = tot[5] / (1 + molal[14] / psc[5] + 
                                   molal[13] / psc[6] + 
                                   molal[19] / psc[11])
    molal[17] = molal[5] * molal[14] / psc[5]
    molal[18] = molal[5] * molal[13] / psc[6]
    molal[23] = molal[5] * molal[19] / psc[11]
    molal[10] = tot[10] / (1 + psc[12] / molal[11])
    molal[24] = tot[10] / (1 + molal[11] / psc[12])
    molal[7] = tot[7] / (1 + molal[11] / psc[13])
    molal[25] = molal[7] * molal[11] / psc[13]
    
    sc, cmax = 0, 0
    for i in range(1, n+1):
        if nch[i] > 0:
            sc += molal[i] * nch[i]
            if molal[i] * nch[i] > cmax:
                cmax = molal[i] * nch[i]
                icat = i
    
    sa, amax = 0, 0
    for i in range(1, n+1):
        if nch[i] < 0:
            sa += molal[i] * -nch[i]
            if molal[i] * -nch[i] > amax:
                amax = molal[i] * -nch[i]
                iani = i
    
    delta = sc - sa
    molal[icat] = molal[icat] - delta / 2 / nch[icat]
    molal[iani] = molal[iani] + delta / 2 / (-nch[iani])
    
    sc = (molal[1] + molal[2] + molal[3] + 
          molal[4] * 2 + molal[5] * 2 + molal[11] + 
          molal[18] + molal[22] + molal[23])
    sa = (molal[6] + molal[7] * 2 + molal[8] + 
          molal[12] + molal[13] + molal[14] * 2 +
          molal[19] + molal[20] + molal[21] * 2 +
          molal[24] + molal[25])
    
    if verbose:
        print("\nSum of cations = {} corrected for {}".format(sc, 
                                                              aq_S[icat]))
        print("Sum of anions = {} corrected for {}\n".format(sa, 
                                                             aq_S[iani]))
        
    s = 0
    for i in range(1, n+1):
        s += molal[i] * atom[i]
    
    if units == "molar":
        ee = 1000 / (1000 * dens - s)
        
        for i in range(1, ntot+1):
            if i != 11:
                tot[i] = tot[i] * ee
        for i in range(1, n+1):
            if i != 11:
                molal[i] = molal[i] * ee

    elif units == "molal":
        ee = 1
    
    stdi = s * ee
    
    return(molal, stdi, ee)

def eql_actp(molal, temp, gact=np.zeros(26), ndepact=0, nc=None, na=None,
             nn=None, nzc=None, nza=None, ap0=None, bp0=None, mh2o=55.51, 
             b0=None, b1=None, b2=None, c0=None, tc=None, ta=None, lc=None, 
             la=None, sc=None, sa=None, xi=None):
    
    """ Function for calculating the activity coefficients using Pitzer 
    Equations. 
    
    For the first iteration of this function (ndepact=0), the only parameters 
    needed are molalities and the temperature of solution. The function will
    read data from the "coefft4" file and calculate the coefficients for 
    all chemical species necessary to calculate the activity coefficients 
    using Pitzer equations.
    
    For the second iteration onward, this function should be called with all
    parameters, so as to utilize the existing data arrays rather than 
    repeating the process of reading and calculating all the necessary inputs
    for the Pitzer equations."""
    
    # Make local arrays to store the molalities of the cations, anions, and
    # water species
    c = np.zeros(10)
    a = np.zeros(12)
    h = np.zeros(4)
    
    # Move the molalities from the input array "molal" to their respective
    # local arrays
    c[1] = molal[1]
    c[2] = molal[2]
    c[3] = molal[3]
    c[4] = molal[4]
    c[5] = molal[22]
    c[6] = molal[5]
    c[7] = molal[18]
    c[8] = molal[23]
    c[9] = molal[11]
    a[1] = molal[6]
    a[2] = molal[7]
    a[3] = molal[25]
    a[4] = molal[12]
    a[5] = molal[14]
    a[6] = molal[13]
    a[7] = molal[24]
    a[8] = molal[19]
    a[9] = molal[20]
    a[10] = molal[21]
    a[11] = molal[8]
    h[1] = molal[10]
    h[2] = molal[9]
    h[3] = molal[15]

    # First iteration: read all data from the "coefft4" file and calculate
    # inputs for Pitzer equations
    if ndepact == 0:
        text = read_file("coefft4")
    
        (nc, na, nn) = (int(i) for i in text[0])
        nzc, nza = np.zeros(nc+1), np.zeros(na+1)
        (b0, b1, b2, c0) = (np.zeros((nc+1, na+1), dtype=np.float64), 
                            np.zeros((nc+1, na+1), dtype=np.float64), 
                            np.zeros((nc+1, na+1), dtype=np.float64), 
                            np.zeros((nc+1, na+1), dtype=np.float64))
        sc = np.zeros((nc+1, nc+1, na+1), dtype=np.float64)
        sa = np.zeros((na+1, na+1, nc+1), dtype=np.float64)
        tc = np.zeros((nc+1, nc+1), dtype=np.float64)
        ta = np.zeros((na+1, na+1), dtype=np.float64)
        lc = np.zeros((nn+1, nc+1), dtype=np.float64)
        la = np.zeros((nn+1, na+1), dtype=np.float64)
        xi = np.zeros((nn+1, nc+1, na+1), dtype=np.float64)
        
        nzc[1:] = [int(i[1]) for i in text[1:1+nc]]
        nza[1:] = [int(i[1]) for i in text[1+nc:1+nc+na]]
        
        (at, bt, ct, dt, et) = (np.float64(i.replace('d', 'e')) 
                                for i in text[1+nc+na][1:])
        
        ap0 = temperature(temp, at, bt, ct, dt, et)
        
        index = 1+nc+na
        
        for i in range(1, nc+1):
            for j in range(1, na+1):
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                b0[i, j] = temperature(temp, at, bt, ct, dt, et)
                
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                b1[i, j] = temperature(temp, at, bt, ct, dt, et)
                
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                b2[i, j] = temperature(temp, at, bt, ct, dt, et)
                
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                c0[i, j] = temperature(temp, at, bt, ct, dt, et)
        
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                tc[i, j] = temperature(temp, at, bt, ct, dt, et)
                tc[j, i] = tc[i, j]
                
        for i in range(1, na):
            for j in range(i+1, na+1):
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                ta[i, j] = temperature(temp, at, bt, ct, dt, et)
                ta[j, i] = ta[i, j]
                
        for k in range(1, nc):
            for i in range(k+1, nc+1):
                for j in range(1, na+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    sc[k, i, j] = temperature(temp, at, bt, ct, dt, et)
                    sc[i, k, j] = sc[k, i, j]
        
        for k in range(1, na):
            for i in range(k+1, na+1):
                for j in range(1, nc+1):
                    index += 1
                    (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                            for k in text[index][1:])
                    sa[k, i, j] = temperature(temp, at, bt, ct, dt, et)
                    sa[i, k, j] = sa[k, i, j]
                    
        for i in range(1, nn+1):
            for j in range(1, nc+1):
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                lc[i, j] = temperature(temp, at, bt, ct, dt, et)
        
        for i in range(1, nn+1):
            for j in range(1, na+1):
                index += 1
                (at, bt, ct, dt, et) = (np.float64(k.replace('d', 'e')) 
                                        for k in text[index][1:])
                la[i, j] = temperature(temp, at, bt, ct, dt, et)
        
        for k in range(1, nn+1):
            for i in range(1, nc+1):
                for j in range(1, na+1):
                    xi[k, i, j] = 0
                    
        xi[2, 9, 1] = -0.0102
        xi[2, 1, 2] = 0.046
    
    # Initialize arrays for the Pitzer summations
    ec, ea = (np.zeros((nc+1, nc+1), dtype=np.float64), 
              np.zeros((na+1, na+1), dtype=np.float64))
    fc, fa = (np.zeros((nc+1, nc+1), dtype=np.float64), 
              np.zeros((na+1, na+1), dtype=np.float64))
    xc, xa = (np.zeros((nc+1, nc+1), dtype=np.float64), 
              np.zeros((na+1, na+1), dtype=np.float64))
    pp, qp = (np.zeros((nc+1, nc+1), dtype=np.float64), 
              np.zeros((na+1, na+1), dtype=np.float64))
    p, q = (np.zeros((nc+1, nc+1), dtype=np.float64), 
            np.zeros((na+1, na+1), dtype=np.float64))
    pf, qf = (np.zeros((nc+1, nc+1), dtype=np.float64), 
              np.zeros((na+1, na+1), dtype=np.float64))
    cc, bf = (np.zeros((nc+1, na+1), dtype=np.float64), 
              np.zeros((nc+1, na+1), dtype=np.float64))
    b, bp = (np.zeros((nc+1, na+1), dtype=np.float64), 
             np.zeros((nc+1, na+1), dtype=np.float64))
    gc, ga, gn = (np.zeros(nc+1, dtype=np.float64), 
                  np.zeros(na+1, dtype=np.float64), 
                  np.zeros(nn+1, dtype=np.float64))
    
    bp0 = 1.2e0
    
    u, z = 0, 0
    
    u += np.sum(c * nzc ** 2)
    z += np.sum(c * nzc)
    u += np.sum(a * nza ** 2)
    z += np.sum(a * nza)

    fi = u / 2
    fj = np.sqrt(fi)
    u = 6 * ap0 * fj        
    
    for i in range(1, nc):
        for j in range(i+1, nc+1):
            
            if nzc[i] == nzc[j]:
                ec[i, j] = 0
                fc[i, j] = 0
                
            else:
                xc[i, j] = 2 * u
                xc[i, i] = nzc[i] ** 2 * u
                xc[j, j] = nzc[j] ** 2 * u
                ec[i, j] = ((j0(xc[i, j]) - j0(xc[i, i]) / 2 - 
                            j0(xc[j, j]) / 2) / fi / 2)
                fc[i, j] = ((xc[i, j] * j1(xc[i, j]) - xc[i, i] * 
                            j1(xc[i, i]) / 2 - xc[j, j] * 
                            j1(xc[j, j]) / 2) / fi ** 2 / 4 - 
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
                ea[i, j] = (j0(xa[i, j]) - j0(xa[i, i]) / 2 -
                            j0(xa[j, j]) / 2) / fi / 2
                fa[i, j] = ((xa[i, j] * j1(xa[i, j]) - xa[i, i] * 
                            j1(xa[i, i]) / 2 - xa[j,j] * 
                            j1(xa[j, j]) / 2) / fi ** 2 / 4 - 
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
                       (g0(v)) + b2[i, j] * (g0(w)))
            bp[i, j] = (b1[i, j] * (g1(v)) / 
                        fi + b2[i, j] * (g1(w)) / fi)
    
    f = -ap0 * (fj / (1 + bp0 * fj) + 2 / bp0 * np.log(1 + bp0 * fj))
    
    i = np.repeat(np.arange(1, nc+1), na)
    j = np.tile(np.arange(1, na+1), nc)
    f += np.sum(c[i] * a[j] * bp[i, j])
    
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
            u += c[i] * ((b[i, jj]) * 2 + z * cc[i, jj])
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
        
    u = -ap0 * fi ** 1.5e0 / (1 + bp0 * fj)
    
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
    aw = np.exp(-s * co / mh2o)
    gact[1] = gc[1]
    gact[2] = gc[2]
    gact[3] = gc[3]
    gact[4] = gc[4]
    gact[22] = gc[5]
    gact[5] = gc[6]
    gact[18] = gc[7]
    gact[23] = gc[8]
    gact[11] = gc[9]
    gact[6] = ga[1]
    gact[7] = ga[2]
    gact[25] = ga[3]
    gact[12] = ga[4]
    gact[14] = ga[5]
    gact[13] = ga[6]
    gact[24] = ga[7]
    gact[19] = ga[8]
    gact[20] = ga[9]
    gact[21] = ga[10]
    gact[8] = ga[11]
    gact[10] = aw * aw * gn[1] ** np.log(10)
    gact[9] = gn[2]
    gact[15] = gn[3]
    gact[16] = 1
    gact[17] = 1
    ndepact = 1
    
    return(gact, nc, na, nn, nzc, nza, ndepact, ap0, bp0, mh2o, b0, b1, b2, 
           c0, tc, ta, lc, la, sc, sa, xi, fi, aw)

def iterate_activities(temp, dens, molal, tot, ph, units, n, act, kmat, z, 
                       zz, psc, n0, ica, xx, nconv, eps, pk, pkstd, atom, 
                       stdi, ee, ntot, dil, cat, nchcat, nch, ani, nchani, 
                       verbose=True):
    
    iterate = True
    while iterate:
        nu = 1
        ncompt = 0
        while nu != 0:
            (gact, nc, na, nn, nzc, nza, ndepact, 
             ap0, bp0, mh2o, b0, b1, b2, c0, tc, 
             ta, lc, la, sc, sa, xi, fi, aw) = eql_actp(molal, temp)
            
            for i in range(1, n+1):
                act[i] = molal[i] * gact[i]
                
            act[0] = aw
            tot[11] = (10 ** (-ph)) / gact[11]
            act[11] = 10 ** (-ph)
            
            for i in range(1, 13):
                for j in range(1, n+1):
                    if molal[j] != 0:
                        z[i, j] = kmat[i, j]
                u = 0
                for j in range(1, n+1):
                    u += kmat[i, j] * molal[j]
                zz[i] = tot[i] - u
            
            for i in range(13, n+1):
                for j in range(1, n+1):
                    if molal[j] != 0:
                        z[i, j] = kmat[i, j] / molal[j]
                    elif molal[j] == 0:
                        z[i, j] = 0
                u = 0
                for j in range(0, n+1):
                    if act[j] > 0:
                        u += kmat[i, j] * np.log(act[j])
                zz[i] = np.log(psc[i-12]) - u
            
            for k in range(1, n0+1):
                if tot[k] == 0 and k != 12:
                    ica[k] = 0
                    for i in range(k+1, n+1):
                        if kmat[i, k] != 0:
                            ica[i] = 0
                    for j in range(k+1, n+1):
                        if kmat[k, j] != 0:
                            ica[j] = 0
            ni, nj = n, n
            for k in range(n, 0, -1):
                if ica[k] == 0:
                    for i in range(k, ni):
                        for j in range(1, nj+1):
                            z[i, j] = z[i+1, j]
                        zz[i] = zz[i+1]    
                    ni -= 1
                    for j in range(k, nj):
                        for i in range(1, ni+1):
                            z[i, j] = z[i, j+1]
                    nj -= 1
            for k in range(2, ni+1):
                for i in range(k, ni+1):
                    if z[i, k-1] != 0:
                        u = z[k-1, k-1] / z[i, k-1]
                        for j in range(k, ni+1):
                            z[i, j] = z[k-1, j] - z[i, j] * u
                        zz[i] = zz[k-1] - zz[i] * u
            xx[ni] = zz[ni] / z[ni, ni]
            
            for i in range(ni-1, 0, -1):
                s = 0
                for j in range(i+1, ni+1):
                    s += z[i, j] * xx[j]
                xx[i] = (zz[i] - s) / z[i, i]
                
            for k in range(1, n+1):
                if ica[k] == 0:
                    for i in range(ni, k-1, -1):
                        xx[i+1] = xx[i]
                    xx[k] = 0
                    ni += 1
            
            ncompt += 1
            if verbose:
                print("iteration molalities {}".format(ncompt))
            
            if ncompt >= 100:
                for i in range(1, n+1):
                    if molal[i] + xx[i] / nconv < 0:
                        print("the equation set diverges: end of program")
                        return()
            
            for i in range(1, n+1):
                if molal[i] + xx[i] / nconv < 0:
                    molal[i] = eps
                else:
                    molal[i] += xx[i] / nconv
            
            nu = 0
            for i in range(1, n+1):
                if ica[i] == 1:
                    if (200 * np.abs(xx[i] / nconv / 
                                     (2 * molal[i] - 
                                      xx[i] / nconv)) > pk):
                        nu = 1
        
        std = 0
        for i in range(0, n+1):
            std += molal[i] * atom[i]
        
        if verbose:
            print("tdsi = {}".format(stdi))
            print("tds = {}".format(std))
        
        if (np.abs(std-stdi)/
           (std+stdi)*200 < pkstd):
            iterate = False
            break
        
        else:
            if units == "molar":
                ef = (1000 + std) / dens / 1000
                for i in range(1, ntot+1):
                    if i != 11:
                        tot[i] = tot[i] / ee * ef
                for i in range(0, n+1):
                    if i != 11:
                        molal[i] = molal[i] / ee * ef
                ee = ef
            
            if verbose:
                print("iteration TDS")
            stdi = std
            
    if units == "molal" and dil == 0:
        cat[1:6] = tot[1:6]
        nchcat[1:6] = nch[1:6]
        ani[1:4] = tot[6:9]
        nchani[1:4] = -nch[6:9]
        ani[4] = molal[12]
        ani[5] = molal[14] + molal[16] + molal[17]
        nchani[4] = -nch[12]
        nchani[5] = -nch[14]
        
        dens = density()
        
