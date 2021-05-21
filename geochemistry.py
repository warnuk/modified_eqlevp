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

def iterate_molalities(molal, temp, ndepact=0, n=25, verbose=True):
    act = np.zeros(n+1)
    
    iterate = True
    while iterate:
        nu = 1
        ncompt = 0
        while nu != 0:
            if ndepact == 0:
                (gact, nc, na, nn, nzc, nza, ndepact, ap0, bp0, mh2o, b0, b1, 
                 b2, c0, tc, ta, lc, la, sc, sa, xi, fi, aw) = actp(molal, 
                                                                    temp)
            else:
                (gact, nc, na, nn, nzc, nza, 
                 ndepact, ap0, bp0, mh2o, b0, 
                 b1, b2, c0, tc, ta, lc, la, 
                 sc, sa, xi, fi, aw) = actp(molal, temp, gact, ndepact, nc, 
                                            na, nn, nzc, nza, ap0, bp0, mh2o, 
                                            b0, b1, b2, c0, tc, ta, lc, la, 
                                            sc, sa, xi)
            
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
            if unit_S == "molar":
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
            
    if unit_S == "molal" and dil == 0:
        cat[1:6] = tot[1:6]
        nchcat[1:6] = nch[1:6]
        ani[1:4] = tot[6:9]
        nchani[1:4] = -nch[6:9]
        ani[4] = molal[12]
        ani[5] = molal[14] + molal[16] + molal[17]
        nchani[4] = -nch[12]
        nchani[5] = -nch[14]
        
        density()  
        
def iterate_pco2(self, verbose=True):
    
    # calculate pco2 and ph
    calculate_pCO2(verbose=verbose)
    
    # if user specifies a pco2, redo the calculations
    if pco2_S == "y":
            
        while np.abs(po - pc) > 0.01:
            if po < pc and poa < pc:
                ph -= dph
            if po < pc and poa > pc:
                dph = dph / 2
                ph -= dph
            if po > pc and poa > pc:
                ph += dph
            if po > pc and poa < pc:
                dph = dph/2
                ph += dph
            poa = po
            
            iterate_molalities(verbose=verbose)
            
            calculate_pCO2(verbose=verbose)
    
    # if pk is greater than pkf, 
    if pk > pkf:
        pk = pkf
        
        if verbose:
            print("last iteration")
        
        iterate_molalities(verbose=verbose)
        calculate_pCO2(verbose=verbose)
        
    
def calculate_pCO2(act, psc, pco2_S, diltot, ph, pc, verbose=True):
    
    po = np.log10(act[15] / psc[14])
    
    if pco2_S == "":
        if diltot == 1:
            poinit = po
            phinit = ph
        
        if verbose:
            print("LOG PCO2 = {}".format(po))
        
        if pc:
            po0 = po
            ph0 = ph
            pco2_S = "y"
            
        elif not pc:
            pco2_S = "n"
            
        
    if pco2_S == "y":
        if verbose:
                print("\nLog(PCO2) selected = {}".format(pc))
                print("Log(PCO2) calculated = {}\n".format(po))
                
    return(pco2_S, po, po0, ph0, poinit, phinit)
            
            
def dilute_solution(tot, nch, molal, cat, ani, nchcat, nchani, dil, diltot, 
                    pk0, verbose=True, output=True):
    
    if dil > 1:
        diltot += dil
        pco2_S = ""
        pk = pk0
        dph = 0.2
        
        tot[np.arange(0, 13) != 11] = (
            tot[np.arange(0, 13) != 11] / dil)
        
        cat[1:6] = tot[1:6]
        nchcat[1:6] = nch[1:6]
        ani[1:4] = tot[6:9]
        nchani[1:4] = -nch[6:9]
        
        ani[4] = molal[12] / dil
        ani[5] = (molal[14] + molal[16] + 
                       molal[17]) / dil
        nchani[4] = -nch[12]
        nchani[5] = -nch[14]
        
        unit_S = "molal"
        
        density()
        
        # Recalculate the molalities after dilution
        calculate_molalities(verbose=verbose)
    
        # calculate activity coefficients // 500 Loop
        iterate_molalities(verbose=verbose)
        
        iterate_pco2(verbose=verbose)
        
        calculate_alkalinities()
        
        print_screen(verbose=verbose, output=output)
        
        saturation_state(verbose=verbose)
        
    return()

def invar(self):
    kinvariant = 0
    ncm = 14
    nbmin = 10
    ninvar = nbmin + 3
    kinv = np.zeros(ninvar+1)
    minv_S = np.empty(ninvar+1, dtype=np.object)
    psminv = np.zeros(ninvar+1)
    winv = np.zeros((ninvar+1, ncm+1))
    minvar_S = np.zeros(ninvar+1)
    psminvar = np.zeros(ninvar+1)
    t0 = np.zeros((ninvar+1, ninvar+1))
    t1 = np.zeros((ninvar+1, ninvar+1))
    t2 = np.zeros((ninvar+1, ninvar+1))
    t3 = np.zeros((ninvar+1, ncm+1))
    t4 = np.zeros((ncm+1, ncm+1))
    tt4 = np.zeros(ncm+1)
    for k in range(1, 4):
        psminv[k] = np.log10(psc[k])
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
            psminv[n1] = np.log10(psol[k])
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
                kinvar[n3] = kinv[kk]
                minvar_S[n3] = minv_S[kk]
                psminvar[n3] = psminv[kk]
                for j in range(0, ncm+1):
                    t3[n3, j] = winv[kk, j]
        if n3 == 0:
            kinvariant = -1
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
                    kinvariant = -2
                else:
                    kinvariant = n3
                    for i in range(kinvariant, 0, -1):
                        if kinvar[i] == 0:
                            for k in range(1, kinvariant):
                                kinvar[k] = kinvar[k+1]
                            kinvariant -= 1
            elif np.abs(t4[n4, n4]) <= epsilon:
                kinvariant = -2


def run_eql(self, verbose=True, output=False, classic=False, dil=1):
    if verbose:
        print("\nThis is EQL..............\n")
        print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')
    
    # initial charge balance
    charge_balance(verbose=verbose)

    # intial molality calculation
    calculate_molalities(verbose=verbose)
    
    # calculate activity coefficients // 500 Loop
    iterate_molalities(verbose=verbose)
    
    iterate_pco2(verbose=verbose)
    
    calculate_alkalinities()
    
    print_screen(verbose=verbose, output=output)
    
    saturation_state(verbose=verbose)
    
    dilute_solution(dil=dil)
    
    # Modify mineral database
    
    # Modify convergence limits
    
    # Add heading to chemistry file
    
    # Write transfer information to stockage file
    
    # Run EVP
    
def saturation_state(self, verbose=True):
    kinvar = np.zeros(nt+1)
    
    for k in range(1, nm+1):
        lmin[k] = 0
    
    if verbose:
        print("The initial solution is oversaturated in {} mineral(s)" \
              " of the data base MURTF3:".format(nwm))
        print()
        
        oversat = []
        for k in range(1, nm+1):
            if pai[k] / psol[k] >= 1 and nwmin[k] == 1:
                oversat.append([mineral_S[k], pai[k] / psol[k]])
                lmin[k] = 1
        oversat = pd.DataFrame(oversat)
        with pd.option_context('display.max_rows', None, 
                                           'display.max_columns', None):
            print(oversat.to_string(index=False, header=False))
        print()
        
    if nwm > nminer:
        print("VIOLATION OF THE PHASE RULE:")
        print("The maximum number of minerals " \
              "allowed is {}".format(nminer))
        print("The evaporation program cannot start with this paragenesis")
        print()
    elif nwm <= nminer:
        invar()
        if kinvariant > 0:
            print("The activity of water is constrained by:")
            print()
            for k in range(1, kinvariant+1):
                print(mineral_S[kinvar[k]])
            print()
        elif kinvariant == -1:
            print("System in thermodynamic desequilibrium")
            print("The activity of water is constrained at "\
                  "different values")
            print("by more than one mineral assemblage")
            print()
        elif kinvariant == -2:
            print("System in thermodynamic desequilibrium:")
            print("inconsistant mineral assemblage")
            print()
        elif kinvariant == 0:
            print("No invariant paragensis detected")
            print()
        
        if kinvariant != 0:
            print("The evaporation program cannot start with this "\
                  "paragenesis")
            print()
                
        if kinvariant == 0 and nwmp > 0:
            print("The solution is close to saturation in {} mineral(s) "\
                  "of the data base MURTF3:".format(nwmp))
            close_saturation = []
            for k in range(1, nm+1):
                if (pai[k] / psol[k] >= 0.9 and 
                    pai[k] / psol[k] < 1 and nwmin[k] == 1):
                    
                    close_saturation.append([mineral_S[k], pai[k] / psol[k]])
                    lmin[k] = 1
                    
            close_saturation = pd.DataFrame(close_saturation)
            with pd.option_context('display.max_rows', None, 
                                           'display.max_columns', None):
                print(close_saturation.to_string(index=False, header=False))
            print()
            
            invar()
            
            if kinvariant > 0:
                print("At the start of evaporation, the activity of "\
                      "water may be constrained by: ")
                
                for k in range(1, kinvariant+1):
                    print(mineral_S[kinvar[k]])
                print()
            elif kinvariant == -1:
                print("System in thermodynamic desequilibrium")
                print("The activity of water is constrained at "\
                      "different values")
                print("by more than one mineral assemblage")
            elif kinvariant == -2:
                print("System in thermodynamic desequilibrium:")
                print("inconsistant mineral assemblage")
                print()
            elif kinvariant == 0:
                print("No invariant paragensis detected")
                print()
            if kinvariant != 0:
                print("If the evaporation program does not start")
                print("slightly dilute the solution again")
                
            

def calculate_alkalinities(self):
    alcar = molal[12] + 2 * (molal[14] + molal[16] + molal[17])
    albor = molal[19] + molal[20] + 2 * molal[21] + molal[22] + molal[23]
    alsil = molal[24]
    aloh = molal[13] + molal[18]
    alh = -1 * molal[11] - molal[25]
    altest = alcar + albor + alsil + aloh + alh
    
def print_screen(self, verbose=True, output=False):
    outfile = "{}.log".format(label)
    if output:
        with open(outfile, "r+") as file:
            file.truncate(0)
            file.close()
    
    for i in range(1, n0+1):
        totinit[i] = 0
        for j in range(1, n+1):
            totinit[i] += molal[j] * kmat[i, j]
            
    df1 = pd.DataFrame({"molality": molal[1:n0+1], 
                        "act coeff": gact[1:n0+1],
                        "activity": act[1:n0+1],
                        "molal tot": totinit[1:n0+1]},
                       index =  aq_S[1:n0+1])
    df1 = df1.loc[df1["molality"] != 0]
    df2 = pd.DataFrame({"molality": molal[n0+1:n+1],
                        "act coeff": gact[n0+1:n+1],
                        "activity": act[n0+1:n+1]},
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

    text = ["ELECTRICAL BALANCE     = {} % corrected on {} and {}".format(dca, aq_S[icat], aq_S[iani]),
            "TOTAL DISSOLVED SOLIDS = {} g/kg(H2O)".format(std),
            "MOLAL/MOLAR FACTOR     = {}".format(1000 * dens / (1000 + std)),
            "DENSITY                = {}".format(dens),
            "IONIC STRENGTH         = {}".format(fi),
            "WATER ACTIVITY         = {}".format(aw),
            "CARBONATE ALKALINITY   = {}".format(alcar),
            "BORATE ALKALINITY      = {}".format(albor),
            "SILICATE ALKALINITY    = {}".format(alsil),
            "OH ALKALINITY          = {}".format(aloh),
            "H ALKALINITY           = {}".format(alh),
            "\nTOTAL ALKALINITY       = {}".format(altest),
            "init. alk.             = {}".format(tot[12])]
    
    if pco2_S == "" or pco2_S == "n":
        text.insert(6, "pH                     = {}".format(ph))
        text.insert(7, "LOG PCO2               = {}".format(po))
    elif pco2_S == "y":
        text.insert(6, "INITIAL LOG PCO2       = {}".format(poinit))
        text.insert(7, "INITIAL pH             = {}".format(phinit))
        text.insert(8, "CALCULATED LOG PCO2    = {}".format(po))
        text.insert(9, "CALCULATED pH          = {}".format(ph))
    
    if diltot > 1:
        text.insert(6, "DILUTION               = {}".format(diltot))
        
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

    condition = np.where(ica==1)[0][np.where(np.where(ica==1)[0] <= 12)]
    u = np.zeros(condition.shape[0])
    for k in range(0, u.shape[0]):
        i = condition[k] * np.ones(n, dtype=int)
        j = np.arange(1, n+1)
        u[k] = np.sum(molal[j] * kmat[i, j])
        
    tests = tot[condition]
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
        if ica[i] == 1:
            for j in range(0, n+1):
                if act[j] != 0:
                    u += kmat[i, j] * np.log10(act[j])
            v = np.log10(psc[i-12])
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
        
    nwm = 0
    nwmp = 0
    
    pai = np.zeros(nm+1)
    data_lst = []
    
    for k in range(1, nm+1):
        pai[k] = 1
        for i in range(0, nt+1):
            pai[k] = pai[k] * act[i] ** wmin[k, i]
        if pai[k] != 0:
            if nwmin[k] == 0:
                zone_S = " " + mineral_S[k].lower()
            elif nwmin[k] == 1:
                zone_S = "*" + mineral_S[k].lower()
            x_S = " "
            if pai[k] / psol[k] >= 1 and nwmin[k] == 1:
                nwm += 1
                x_S = "*"
            elif pai[k] / psol[k] >= 0.9 and pai[k] / psol[k] < 1 and nwmin[k] == 1:
                nwmp += 1
            data_lst.append([zone_S, psol[k], pai[k], pai[k] / psol[k], x_S])
    
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
