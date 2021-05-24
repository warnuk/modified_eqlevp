# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:07:27 2021

@author: Will Arnuk
"""
import numpy as np

def read_file(file):
    with open(file, 'r+') as file:
        text = list(filter(None, file.read().split('\n')))
        file.close()
        values = [line.split(',') for line in text]
    return(values)

def temperature(temp, at, bt, ct, dt, et):
    return((at + bt * temp + ct * temp ** 2 + 
            dt * temp ** 3 + et * temp ** 4))

def j0(x):
    ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
    j0 = x / (4. + ya * x ** yb * np.exp(yc * x ** yd))
    return(j0)
        
def j1(x):
    ya, yb, yc, yd = 4.581, -0.7237, -0.012, 0.528
    j1 = ((4. + ya * x ** yb * (1. - yb - yc * yd * x ** yd) * 
      np.exp(yc * x ** yd)) / (4. + ya * x ** yb * 
                               np.exp(yc * x ** yd)) ** 2)
    return(j1)

def g0(x):
    g0 = 2. * (1. - (1. + x) * np.exp(-x)) / x ** 2
    return(g0)
        
def g1(x):
    g1 = -2. * (1. - (1. + x + x ** 2 / 2.) * np.exp(-x)) / x ** 2
    return(g1)

def calculate_molalities(tot, psc, pk, nch, aq_S, atom, dens, ntot, unit_S, n=25, verbose=True):
    molal = np.zeros(n+1)
    
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
    
    if unit_S == "molar":
        ee = 1000 / (1000 * dens - s)
        
        for i in range(1, ntot+1):
            if i != 11:
                tot[i] = tot[i] * ee
        for i in range(1, n+1):
            if i != 11:
                molal[i] = molal[i] * ee

    elif unit_S == "molal":
        ee = 1
    
    stdi = s * ee
    return(stdi)

def actp(molal, temp, gact=np.zeros(26), ndepact=0, nc=None, na=None,
         nn=None, nzc=None, nza=None, ap0=None, bp0=None, mh2o=55.51, b0=None,
         b1=None, b2=None, c0=None, tc=None, ta=None, lc=None, la=None,
         sc=None, sa=None, xi=None):
    
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

def density(unit_S, cat, ani, nchcat, nchani):
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
    
    if unit_S == "molar":
        dens = 1
        u = np.sum(ani[1:na+1])
            
        i = np.repeat(np.arange(1, nc+1), na)
        j = np.tile(np.arange(1, na+1), nc)
        s[i, j] = ((nchcat[i] + nchani[j]) / 2 * cat[i] * 
                   ani[j] / nchcat[i] / nchani[j] / u)
        dens += np.sum(ao[i, j] * s[i, j] + bo[i, j] * s[i, j] ** 2)
        
    elif unit_S == "molal":
        dens = 1
        u = np.sum(ani[1:na+1] * nchani[1:na+1])
        i = np.repeat(np.arange(1, nc+1), na)
        j = np.tile(np.arange(1, na+1), nc)
        s[i, j] = ((nchcat[i] + nchani[j]) / 2 * cat[i] * 
                   ani[j] / u)
        dens += np.sum(au[i, j] * s[i, j] + bu[i, j] * s[i, j] ** 2)
        
    return(dens)

def evp_actp(mol, temp, gact=np.zeros(26), ndepact=0, nc=None, na=None,
             nn=None, nzc=None, nza=None, ap0=None, bp0=None, mh2o=55.51, 
             b0=None, b1=None, b2=None, c0=None, tc=None, ta=None, lc=None, 
             la=None, sc=None, sa=None, xi=None):
    
    """ Function for calculating the activity coefficients using Pitzer 
    Equations. 
    
    For the first iteration of this function (ndepact=0), the only parameters 
    needed are molities and the temperature of solution. The function will
    read data from the "coefft4" file and calculate the coefficients for 
    all chemical species necessary to calculate the activity coefficients 
    using Pitzer equations.
    
    For the second iteration onward, this function should be called with all
    parameters, so as to utilize the existing data arrays rather than 
    repeating the process of reading and calculating all the necessary inputs
    for the Pitzer equations."""
    
    # Make local arrays to store the molities of the cations, anions, and
    # water species
    cat = np.zeros(10)
    ani = np.zeros(12)
    h = np.zeros(4)
    
    # Move the molities from the input array "mol" to their respective
    # local arrays
    cat[1] = mol[1]
    cat[2] = mol[2]
    cat[3] = mol[3]
    cat[4] = mol[4]
    cat[5] = mol[22]
    cat[6] = mol[5]
    cat[7] = mol[18]
    cat[8] = mol[23]
    cat[9] = mol[13]
    ani[1] = mol[6]
    ani[2] = mol[7]
    ani[3] = mol[25]
    ani[4] = mol[15]
    ani[5] = mol[14]
    ani[6] = mol[12]
    ani[7] = mol[24]
    ani[8] = mol[19]
    ani[9] = mol[20]
    ani[10] = mol[21]
    ani[11] = mol[8]
    h[1] = mol[10]
    h[2] = mol[9]
    h[3] = mol[0]

    cat[1:10] = cat[1:10] * mh2o / mol[11]
    ani[1:12] = ani[1:12] * mh2o / mol[11]
    h[1:4] = h[1:4] * mh2o / mol[11]

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
    
    bp0 = 1.2e0
    
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
    
    ndepact = 1
    
    u, z = 0, 0
    
    u += np.sum(cat * nzc ** 2)
    z += np.sum(cat * nzc)
    u += np.sum(ani * nza ** 2)
    z += np.sum(ani * nza)

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
    
    for i in range(1, nc+1):
        for j in range(1, na+1):
            f += cat[i] * ani[j] * bp[i, j]
            
    # i = np.repeat(np.arange(1, nc+1), na)
    # j = np.tile(np.arange(1, na+1), nc)
    # f += np.sum(cat[i] * ani[j] * bp[i, j])
    
    for i in range(1, nc):
        for j in range(i+1, nc+1):
            f += cat[i] * cat[j] * pp[i, j]

    for i in range(1, na):
        for j in range(i+1, na+1):
            f += ani[i] * ani[j] * qp[i, j]
            
    for ii in range(1, nc+1):
        u = nzc[ii] ** 2 * f
        for j in range(1, na+1):
            u += ani[j] * (b[ii, j] * 2 + z * cc[ii, j])
        for i in range(1, nc+1):
            if i != ii:
                v = 0
                for j in range(1, na+1):
                    v += ani[j] * sc[ii, i, j]
                u += cat[i] * (p[ii, i] * 2 + v)
        for i in range(1, na):
            for j in range(i+1, na+1):
                u += ani[i] * ani[j] * sa[i, j, ii]
        for i in range(1, nc+1):
            for j in range(1, na+1):
                u += cat[i] * ani[j] * cc[i, j] * nzc[ii]
        for i in range(1, nn+1):
            u += h[i] * lc[i, ii] * 2
        for k in range(1, nn+1):
            for j in range(1, na+1):
                u += h[k] * ani[j] * xi[k, ii, j]
        gc[ii] = np.exp(u)
        
    for jj in range(1, na+1):
        u = nza[jj] ** 2 * f
        for i in range(1, nc+1):
            u += cat[i] * ((b[i, jj]) * 2 + z * cc[i, jj])
        for i in range(1, na+1):
            if i != jj:
                v = 0
                for j in range(1, nc+1):
                    v += cat[j] * sa[jj, i, j]
                u += ani[i] * (q[jj, i] * 2 + v)
                    
        for i in range(1, nc):
            for j in range(i+1, nc+1):
                u += cat[i] * cat[j] * sc[i, j, jj]
        for i in range(1, nc+1):
            for j in range(1, na+1):
                u += cat[i] * ani[j] * cc[i, j] * nza[jj]
        for j in range(1, nn+1):
            u += h[j] * la[j, jj]
        for k in range(1, nn+1):
            for i in range(1, nc+1):
                u += h[k] * cat[i] * xi[k, i, jj]
        ga[jj] = np.exp(u)
        
    for k in range(1, nn+1):
        u = 0
        for i in range(1, nc+1):
            u += cat[i] * lc[k, i] * 2
        for j in range(1, na+1):
            u += ani[j] * la[k, j] * 2
        for i in range(1, nc+1):
            for j in range(1, na+1):
                u += cat[i] * ani[j] * xi[k, i, j]
        gn[k] = np.exp(u)
        
    u = -ap0 * fi ** 1.5e0 / (1 + bp0 * fj)
    
    for i in range(1, nc+1):
        for j in range(1, na+1):
            u += cat[i] * ani[j] * (bf[i, j] + z * cc[i, j])
    for i in range(1, nc):
        for j in range(i+1, nc+1):
            v = 0
            for k in range(1, na+1):
                v += ani[k] * sc[i, j, k]
            u += cat[i] * cat[j] * (pf[i, j] + v)
    for i in range(1, na):
        for j in range(i+1, na+1):
            v = 0
            for k in range(1, nc+1):
                v += cat[k] * sa[i, j, k]
            u += ani[i] * ani[j] * (qf[i, j] + v)
    for k in range(1, nn+1):
        for i in range(1, nc+1):
            u += h[k] * cat[i] * lc[k, i]
    for k in range(1, nn+1):
        for j in range(1, na+1):
            u += h[k] * ani[j] * la[k, j]
    for k in range(1, nn+1):
        for i in range(1, nc+1):
            for j in range(1, na+1):
                u += h[k] * cat[i] * ani[j] * xi[k, i, j]
    
    s = 0
    for i in range(1, nc+1):
        s += cat[i]
    for j in range(1, na+1):
        s += ani[j]
    co = 1 + 2 * u / s
    aw = np.exp(-s * co / mh2o)
    gact[0] = gn[3]
    gact[1] = gc[1]
    gact[2] = gc[2]
    gact[3] = gc[3]
    gact[4] = gc[4]
    gact[22] = gc[5]
    gact[5] = gc[6]
    gact[18] = gc[7]
    gact[23] = gc[8]
    gact[13] = gc[9]
    gact[6] = ga[1]
    gact[7] = ga[2]
    gact[25] = ga[3]
    gact[15] = ga[4]
    gact[14] = ga[5]
    gact[12] = ga[6]
    gact[24] = ga[7]
    gact[19] = ga[8]
    gact[20] = ga[9]
    gact[21] = ga[10]
    gact[8] = ga[11]
    gact[10] = aw * aw * gn[1] ** np.log(10)
    gact[9] = gn[2]
    gact[16] = 1
    gact[17] = 1
    gact[11] = aw / mh2o
    ndepact = 1
    
    return(gact, nc, na, nn, nzc, nza, ndepact, ap0, bp0, mh2o, b0, b1, b2, 
           c0, tc, ta, lc, la, sc, sa, xi, fi, aw)