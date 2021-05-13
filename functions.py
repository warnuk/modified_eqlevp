# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:25:26 2021

@author: Will Arnuk
"""
import numpy as np
import pandas as pd

def actp(mol, gact, fi, temp):
    global nc, na, nn, ndepact, nzc, nza, ap0, bp0
    global mh2o, cat, ani, h, b0, b1, b2, c0, c, ta, tc, lc, la, sc, sa, xi
    
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
    
    return(mol, gact, fi, temp)
            
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