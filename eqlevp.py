# -*- coding: utf-8 -*-
"""
Created on Mon May 24 10:43:01 2021

@author: warnu
"""

import numpy as np
import pandas as pd
from datetime import datetime

def read_file(file):
    with open(file, 'r+') as file:
        text = list(filter(None, file.read().split('\n')))
        file.close()
        values = [line.split(',') for line in text]
    return(values)

class simulation:

    def __init__(self, label, temp, dens, ph, na, k,
                 li, ca, mg, cl, so4, no3, b, si, alk):

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

    def run_eql(self, pco2, system, units="molal", dilute=1,
                add_minerals=None, rem_minerals=None, max_salinity=0,
                pkmol=None, pkeq=None, verbose=True, output=True,
                call_evp=True):

        # Set additional object attributes to the
        # parameters used in this simulation
        self.system = system
        self.units = units
        self.dil = dilute
        self.pco2 = pco2
        self.stdmax = max_salinity

        if verbose:
            print("\nThis is EQL..............\n")
            print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')

        self.initialize_eql(system=system, verbose=verbose, output=output)

        self.calculate_molalities(verbose=verbose)

        self.iterate_activities(verbose=verbose)

        self.iterate_pco2(verbose=verbose)

        self.calculate_alkalinities()

        self.print_screen(verbose=verbose, output=output)

        if call_evp:
            self.saturation_state(verbose=verbose)
            self.dilute_solution(dilute=dilute)
            self.modify_database(add_min=add_minerals,
                                 rem_min=rem_minerals,
                                 verbose=verbose, output=output)
            if add_minerals or rem_minerals:
                self.min_S = "murtf0"
            # Modify convergence limits
            if pkmol:
                self.pkmol = pkmol
            else:
                self.pkmol = 0.001
            if pkeq:
                self.pkeq = pkeq
            else:
                self.pkeq = .0000000000001
            # Add heading to chem_file
            self.constit_S = ["label", "fc", "eva", "ds", "ph", "alk"]

            for i in range(1, 9):
                if self.tot[i] > 0:
                    self.constit_S.append(self.aq_S[i])
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
            lines.append(system)
            lines.append(1)
            lines.append(1)
            lines.append('y')
            lines.append(1)
            lines.append(1)
            lines.append(self.units)
            lines.append(self.chem_file)
            lines.append(self.event_file)
            lines.append(self.min_file)
            lines.append(self.min_S)
            lines.append(max_salinity)
            lines.append(self.pkmol)
            lines.append(self.pkeq)

            lines = [str(i) for i in lines]

            with open(self.transfer_file, "w") as file:
                file.write("\n".join(lines))
                file.close()
                    
            # Clear out unnecessary simulation attributs
            old_attributes = ['albor', 'alcar', 'alh', 'aloh', 'alsil', 
                              'altest', 'ani', 'ap0', 'aq_S', 'atom', 
                              'b0', 'b1', 'b2', 'bp0', 'c0', 'cat', 'iani',
                              'icat', 'kmat', 'la', 'lc']
            
            for i in old_attributes:
                delattr(self, i)
             
            # Run EVP
            self.run_evp(verbose=verbose, output=output)

    def initialize_eql(self, system, verbose, output):
        
        # Initialize files used by EQL/EVP
        if output:
                self.log_file = "{}.log".format(self.label)
                self.chem_file = "{}.j{}&".format(self.label, system)
                self.event_file = "{}.j{}@".format(self.label, system)
                self.min_file = "{}.j{}%".format(self.label, system)
                self.transfer_file = "{}.tra".format(self.label)
                
                for filename in [self.log_file, self.chem_file, 
                                 self.event_file, self.min_file, 
                                 self.transfer_file]:
                    with open(filename, "w") as file:
                        file.close()
        
        # Set parameters for simulation
        self.n = 25
        self.ntot = 12
        self.ncomplex = 14
        self.n0 = 10
        self.max_cat = 5
        self.max_an = 5
        self.epsilon = 1e-8

        self.max_constit = 30
        self.max_nr = 10
        self.max_conv=100

        self.ndepact = 0
        self.pco2_S = ""
        self.pk = 0.1
        self.pk0 = self.pk
        self.pkf = 0.0001
        self.pkstd = 0.1
        self.dph = 0.2
        self.diltot = 1
        self.mh2o = 55.51
        self.nconv = 2
        self.eps = 1e-12
        self.poa = 0.

        # Initialize blank arrays
        self.constit_S = np.empty(self.max_constit+1, object)
        self.psc = np.zeros(15)
        self.tot = np.zeros(self.ntot+1)
        self.tot0 = np.zeros(self.ntot+1)
        self.totinit = np.zeros(self.n0+1)
        self.nch = np.zeros(self.n+1)
        self.molal = np.zeros(self.n+1)
        self.act = np.zeros(self.n+1)
        self.gact = np.zeros(self.n+1)
        self.aq_S = np.empty(self.n+1, object)
        self.atom = np.zeros(self.n+1)
        self.kmat = np.zeros((self.n+1, self.n+1))
        self.ica = np.ones(self.n+1)
        self.xx = np.zeros(self.n+1)
        self.z = np.zeros((self.n+1, self.n+1))
        self.zz = np.zeros(self.n+1)
        self.cat = np.zeros(self.max_cat+1)
        self.ani = np.zeros(self.max_an+1)
        self.nchcat = np.zeros(self.max_cat+1)
        self.nchani = np.zeros(self.max_an+1)


        # get concentrations from parent inputs
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

        # Set initial temperature to input temperature
        self.tinit = self.temp

        # Get number of components
        self.nminer = np.count_nonzero(self.tot[1:11])+1

        # Read aqu database
        aqu = read_file("aqu.dat")
        self.aq_S[:] = [line[0] for line in aqu]
        self.atom[:] = [line[1] for line in aqu]
        self.nch[:] = [line[2] for line in aqu]

        # Read kmat from matrice1
        self.kmat[1:,1:] = [line[:] for line in read_file("matrice1")]
        self.kmat[13, 0] = -1
        self.kmat[15, 0] = -1
        self.kmat[20, 0] = 3
        self.kmat[21, 0] = 5

        # Read thermodynamic data for dissociation coefficients from complex3
        complex3 = np.zeros((self.ncomplex+1, 5))
        complex3[1:,:] = [i[1:] for i in read_file("complex3")]
        self.psc[:] = 10 ** (complex3[:,0] +
                             complex3[:,1] / 300 * self.temp +
                             complex3[:,2] / 30000 * self.temp ** 2 +
                             complex3[:,3] / 3000000 * self.temp ** 3 +
                             complex3[:,4] / 300000000 * self.temp ** 4)
        self.psc[0] = 0

        # Read eql mineral database from murtf2
        murtf2 = read_file("murtf2")
        (self.nc, self.na, self.nm) = (int(i) for i in murtf2[0])
        self.nt = self.nc + self.na

        self.wmin = np.zeros((self.nm+1, self.nt+1))
        self.lmin = np.zeros(self.nm+1)
        self.nwmin = np.zeros(self.nm+1)
        self.mineral_S = np.empty(self.nm+1, object)
        self.mu = np.zeros(self.nt+1)
        self.psol = np.zeros(self.nm+1)
        self.pai = np.zeros(self.nm+1)
        self.kinvar = np.zeros(self.nt+1)

        ion_list = murtf2[1:self.nt+2]
        ion_list.insert(0, None)

        for i in range(1, self.nt+2):
            a_S = str(ion_list[i][0]).lower()
            (at, bt, ct, dt, et) = (float(x) for x in ion_list[i][1:])
            for j in range(0, self.nt+1):
                if a_S == str(self.aq_S[j]):
                    self.mu[j] = (at +
                                  bt / 300 * self.temp +
                                  ct / 30000 * self.temp ** 2 +
                                  dt / 3000000 * self.temp ** 3 +
                                  et / 300000000 * self.temp ** 4)

        mineral_list = murtf2[self.nt+2:]
        mineral_list.insert(0, None)
        min_db = np.zeros((self.nm+1, 5))

        self.mineral_S = np.empty(self.nm+1, object)

        for k in range(1, self.nm+1):
            line = mineral_list[k]
            self.mineral_S[k] = line[0]
            ncomp = int(line[1])
            c_ion = np.zeros(ncomp+1)
            nom_ion_S = np.empty(ncomp+1, dtype=object)
            for i in range(1, ncomp+1):
                c_ion[i] = float(line[i*2])
                nom_ion_S[i] = str(line[1+i*2])
                for j in range(0, self.nt+1):
                    x_S = nom_ion_S[i].lower()
                    if x_S == self.aq_S[j]:
                        self.wmin[k, j] = c_ion[i]

            (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
            min_db[k] = np.array([at, bt, ct, dt, et])

        self.mum = (min_db[:,0] +
                    min_db[:,1] / 300 * self.temp +
                    min_db[:,2] / 30000 * self.temp ** 2 +
                    min_db[:,3] / 3000000 * self.temp ** 3 +
                    min_db[:,4] / 300000000 * self.temp ** 4)

        for k in range(1, self.nm+1):
            u = self.mum[k]
            for i in range(0, self.nt+1):
                u -= self.wmin[k, i] * self.mu[i]
            self.psol[k] = np.exp(u)

        # Charge Balance the initial water chemistry
        sc, cmax = 0, 0
        for i in range(1, self.ntot+1):
            if self.nch[i] > 0 and i != 11:
                sc += self.tot[i] * self.nch[i]
                if self.tot[i] * self.nch[i] > cmax:
                    cmax = self.tot[i] * self.nch[i]
                    icat = i
        sa, amax = 0, 0
        for i in range(1, self.ntot+1):
            if self.nch[i] < 0:
                sa += self.tot[i] * -self.nch[i]
                if self.tot[i] * -self.nch[i] > amax:
                    amax = self.tot[i] * -self.nch[i]
                    iani = i

        if sc + sa != 0:
            self.dca = 200 * np.abs(sc-sa) / (sc+sa)
        else:
            self.dca = 0
        delta = sc-sa

        if verbose:
            print("sum of cations = {}".format(sc))
            print("sum of anions = {}".format(sa))
            print("Electrical balance = {} %".format((self.dca * 100 +
                                                      0.5) / 100))

        self.tot[icat] = self.tot[icat] - delta / 2 / self.nch[icat]
        self.tot[iani] = self.tot[iani] + delta / 2 / -self.nch[iani]

        self.tot0[1:13] = self.tot[1:13]

        # Set the default mineral database to murtf3
        self.min_S = "murtf3"
        murtf3 = read_file(self.min_S)

        self.nc, self.na, self.nm0 = [int(i) for i in murtf3[0]]

        self.mineral0_S = np.empty(self.nm0+1, object)

        self.mineral0_S[1:] = [i[0] for i in murtf3[(2+self.nc+self.na):]]
        self.nwmin[np.in1d(self.mineral_S, self.mineral0_S)] = 1

        return()

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
        for i in range(1, self.n+1):
            if self.nch[i] > 0:
                sc += self.molal[i] * self.nch[i]
                if self.molal[i] * self.nch[i] > cmax:
                    cmax = self.molal[i] * self.nch[i]
                    self.icat = i

        sa, amax = 0, 0
        for i in range(1, self.n+1):
            if self.nch[i] < 0:
                sa += self.molal[i] * -self.nch[i]
                if self.molal[i] * -self.nch[i] > amax:
                    amax = self.molal[i] * -self.nch[i]
                    self.iani = i

        delta = sc - sa
        self.molal[self.icat] = (self.molal[self.icat] - delta
                                 / 2 / self.nch[self.icat])
        self.molal[self.iani] = (self.molal[self.iani] + delta
                                 / 2 / (-self.nch[self.iani]))

        sc = (self.molal[1] + self.molal[2] + self.molal[3] +
              self.molal[4] * 2 + self.molal[5] * 2 + self.molal[11] +
              self.molal[18] + self.molal[22] + self.molal[23])
        sa = (self.molal[6] + self.molal[7] * 2 + self.molal[8] +
              self.molal[12] + self.molal[13] + self.molal[14] * 2 +
              self.molal[19] + self.molal[20] + self.molal[21] * 2 +
              self.molal[24] + self.molal[25])

        if verbose:
            print("\nSum of cations = {} "\
                  "corrected for {}".format(sc, self.aq_S[self.icat]))
            print("Sum of anions = {} "\
                  "corrected for {}\n".format(sa, self.aq_S[self.iani]))

        s = 0
        for i in range(1, self.n+1):
            s += self.molal[i] * self.atom[i]

        if self.units == "molar":
            self.ee = 1000 / (1000 * self.dens - s)

            for i in range(1, self.ntot+1):
                if i != 11:
                    self.tot[i] = self.tot[i] * self.ee
            for i in range(1, self.n+1):
                if i != 11:
                    self.molal[i] = self.molal[i] * self.ee

        elif self.units == "molal":
            self.ee = 1

        self.stdi = s * self.ee





    def eql_actp(self):
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
        
    def evp_actp(self):
        cat = np.zeros(10)
        ani = np.zeros(12)
        h = np.zeros(4)

        cat[1] = self.mol[1]
        cat[2] = self.mol[2]
        cat[3] = self.mol[3]
        cat[4] = self.mol[4]
        cat[5] = self.mol[22]
        cat[6] = self.mol[5]
        cat[7] = self.mol[18]
        cat[8] = self.mol[23]
        cat[9] = self.mol[13]
        ani[1] = self.mol[6]
        ani[2] = self.mol[7]
        ani[3] = self.mol[25]
        ani[4] = self.mol[15]
        ani[5] = self.mol[14]
        ani[6] = self.mol[12]
        ani[7] = self.mol[24]
        ani[8] = self.mol[19]
        ani[9] = self.mol[20]
        ani[10] = self.mol[21]
        ani[11] = self.mol[8]
        h[1] = self.mol[10]
        h[2] = self.mol[9]
        h[3] = self.mol[0]
        
        cat[1:10] *= self.mh2o / self.mol[11]
        ani[1:12] *= self.mh2o / self.mol[11]
        h[1:4] *= self.mh2o / self.mol[11]

        if self.ndepact == 0:
            text = read_file("coefft4")

            (self.nc, self.na, self.nn) = (int(i) for i in text[0])
            self.nzc, self.nza = np.zeros(self.nc+1, dtype=int), np.zeros(self.na+1, dtype=int)
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

        u += np.sum(cat * self.nzc ** 2)
        z += np.sum(cat * self.nzc)
        u += np.sum(ani * self.nza ** 2)
        z += np.sum(ani * self.nza)

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
        f += np.sum(cat[i] * ani[j] * bp[i, j])

        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                f += cat[i] * cat[j] * pp[i, j]

        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                f += ani[i] * ani[j] * qp[i, j]

        for ii in range(1, self.nc+1):
            u = self.nzc[ii] ** 2 * f
            for j in range(1, self.na+1):
                u += ani[j] * (b[ii, j] * 2 + z * cc[ii, j])
            for i in range(1, self.nc+1):
                if i != ii:
                    v = 0
                    for j in range(1, self.na+1):
                        v += ani[j] * self.sc[ii, i, j]
                    u += cat[i] * (p[ii, i] * 2 + v)
            for i in range(1, self.na):
                for j in range(i+1, self.na+1):
                    u += ani[i] * ani[j] * self.sa[i, j, ii]
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += cat[i] * ani[j] * cc[i, j] * self.nzc[ii]
            for i in range(1, self.nn+1):
                u += h[i] * self.lc[i, ii] * 2
            for k in range(1, self.nn+1):
                for j in range(1, self.na+1):
                    u += h[k] * ani[j] * self.xi[k, ii, j]
            gc[ii] = np.exp(u)

        for jj in range(1, self.na+1):
            u = self.nza[jj] ** 2 * f
            for i in range(1, self.nc+1):
                u += cat[i] * ((b[i, jj]) * 2 + z * cc[i, jj])
            for i in range(1, self.na+1):
                if i != jj:
                    v = 0
                    for j in range(1, self.nc+1):
                        v += cat[j] * self.sa[jj, i, j]
                    u += ani[i] * (q[jj, i] * 2 + v)

            for i in range(1, self.nc):
                for j in range(i+1, self.nc+1):
                    u += cat[i] * cat[j] * self.sc[i, j, jj]
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += cat[i] * ani[j] * cc[i, j] * self.nza[jj]
            for j in range(1, self.nn+1):
                u += h[j] * self.la[j, jj]
            for k in range(1, self.nn+1):
                for i in range(1, self.nc+1):
                    u += h[k] * cat[i] * self.xi[k, i, jj]
            ga[jj] = np.exp(u)

        for k in range(1, self.nn+1):
            u = 0
            for i in range(1, self.nc+1):
                u += cat[i] * self.lc[k, i] * 2
            for j in range(1, self.na+1):
                u += ani[j] * self.la[k, j] * 2
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += cat[i] * ani[j] * self.xi[k, i, j]
            gn[k] = np.exp(u)

        u = -self.ap0 * self.fi ** 1.5e0 / (1 + self.bp0 * fj)
        for i in range(1, self.nc+1):
            for j in range(1, self.na+1):
                u += cat[i] * ani[j] * (bf[i, j] + z * cc[i, j])
        for i in range(1, self.nc):
            for j in range(i+1, self.nc+1):
                v = 0
                for k in range(1, self.na+1):
                    v += ani[k] * self.sc[i, j, k]
                u += cat[i] * cat[j] * (pf[i, j] + v)
        for i in range(1, self.na):
            for j in range(i+1, self.na+1):
                v = 0
                for k in range(1, self.nc+1):
                    v += cat[k] * self.sa[i, j, k]
                u += ani[i] * ani[j] * (qf[i, j] + v)
        for k in range(1, self.nn+1):
            for i in range(1, self.nc+1):
                u += h[k] * cat[i] * self.lc[k, i]
        for k in range(1, self.nn+1):
            for j in range(1, self.na+1):
                u += h[k] * ani[j] * self.la[k, j]
        for k in range(1, self.nn+1):
            for i in range(1, self.nc+1):
                for j in range(1, self.na+1):
                    u += h[k] * cat[i] * ani[j] * self.xi[k, i, j]

        s = 0
        for i in range(1, self.nc+1):
            s += cat[i]
        for j in range(1, self.na+1):
            s += ani[j]
        co = 1 + 2 * u / s
        self.aw = np.exp(-s * co / self.mh2o)
        self.gact[0] = gn[3]
        self.gact[1] = gc[1]
        self.gact[2] = gc[2]
        self.gact[3] = gc[3]
        self.gact[4] = gc[4]
        self.gact[22] = gc[5]
        self.gact[5] = gc[6]
        self.gact[18] = gc[7]
        self.gact[23] = gc[8]
        self.gact[13] = gc[9]
        self.gact[6] = ga[1]
        self.gact[7] = ga[2]
        self.gact[25] = ga[3]
        self.gact[15] = ga[4]
        self.gact[14] = ga[5]
        self.gact[12] = ga[6]
        self.gact[24] = ga[7]
        self.gact[19] = ga[8]
        self.gact[20] = ga[9]
        self.gact[21] = ga[10]
        self.gact[8] = ga[11]
        self.gact[10] = self.aw * self.aw * gn[1] ** np.log(10)
        self.gact[9] = gn[2]
        self.gact[16] = 1
        self.gact[17] = 1
        self.gact[11] = self.aw / self.mh2o
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

        if self.units == "molar":
            self.dens = 1
            u = np.sum(self.ani[1:na+1])

            i = np.repeat(np.arange(1, nc+1), na)
            j = np.tile(np.arange(1, na+1), nc)
            s[i, j] = ((self.nchcat[i] + self.nchani[j]) / 2 * self.cat[i] *
                       self.ani[j] / self.nchcat[i] / self.nchani[j] / u)
            self.dens += np.sum(ao[i, j] * s[i, j] + bo[i, j] * s[i, j] ** 2)

        elif self.units == "molal":
            self.dens = 1
            u = np.sum(self.ani[1:na+1] * self.nchani[1:na+1])
            i = np.repeat(np.arange(1, nc+1), na)
            j = np.tile(np.arange(1, na+1), nc)
            s[i, j] = ((self.nchcat[i] + self.nchani[j]) / 2 * self.cat[i] *
                       self.ani[j] / u)
            self.dens += np.sum(au[i, j] * s[i, j] + bu[i, j] * s[i, j] ** 2)

    def iterate_activities(self, verbose=True):
        iterate = True
        while iterate:
            self.nu = 1
            self.ncompt = 0
            while self.nu != 0:
                self.eql_actp()

                for i in range(1, self.n+1):
                    self.act[i] = self.molal[i] * self.gact[i]

                self.act[0] = self.aw
                self.tot[11] = (10 ** (-self.ph)) / self.gact[11]
                self.act[11] = 10 ** (-self.ph)

                for i in range(1, 13):
                    for j in range(1, self.n+1):
                        if self.molal[j] != 0:
                            self.z[i, j] = self.kmat[i, j]
                    u = 0
                    for j in range(1, self.n+1):
                        u += self.kmat[i, j] * self.molal[j]
                    self.zz[i] = self.tot[i] - u

                for i in range(13, self.n+1):
                    for j in range(1, self.n+1):
                        if self.molal[j] != 0:
                            self.z[i, j] = self.kmat[i, j] / self.molal[j]
                        elif self.molal[j] == 0:
                            self.z[i, j] = 0
                    u = 0
                    for j in range(0, self.n+1):
                        if self.act[j] > 0:
                            u += self.kmat[i, j] * np.log(self.act[j])
                    self.zz[i] = np.log(self.psc[i-12]) - u

                for k in range(1, self.n0+1):
                    if self.tot[k] == 0 and k != 12:
                        self.ica[k] = 0
                        for i in range(k+1, self.n+1):
                            if self.kmat[i, k] != 0:
                                self.ica[i] = 0
                        for j in range(k+1, self.n+1):
                            if self.kmat[k, j] != 0:
                                self.ica[j] = 0
                ni, nj = self.n, self.n
                for k in range(self.n, 0, -1):
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

                for k in range(1, self.n+1):
                    if self.ica[k] == 0:
                        for i in range(ni, k-1, -1):
                            self.xx[i+1] = self.xx[i]
                        self.xx[k] = 0
                        ni += 1

                self.ncompt += 1
                if verbose:
                    print("iteration molalities {}".format(self.ncompt))

                if self.ncompt >= 100:
                    for i in range(1, self.n+1):
                        if self.molal[i] + self.xx[i] / self.nconv < 0:
                            print("the equation set diverges: end of program")
                            return()

                for i in range(1, self.n+1):
                    if self.molal[i] + self.xx[i] / self.nconv < 0:
                        self.molal[i] = self.eps
                    else:
                        self.molal[i] += self.xx[i] / self.nconv

                self.nu = 0
                for i in range(1, self.n+1):
                    if self.ica[i] == 1:
                        if (200 * np.abs(self.xx[i] / self.nconv /
                                         (2 * self.molal[i] -
                                          self.xx[i] / self.nconv)) > self.pk):
                            self.nu = 1

            self.std = 0
            for i in range(0, self.n+1):
                self.std += self.molal[i] * self.atom[i]

            if verbose:
                print("tdsi = {}".format(self.stdi))
                print("tds = {}".format(self.std))

            if (np.abs(self.std-self.stdi)/
               (self.std+self.stdi)*200 < self.pkstd):
                iterate = False
                break

            else:
                if self.units == "molar":
                    self.ef = (1000 + self.std) / self.dens / 1000
                    for i in range(1, self.ntot+1):
                        if i != 11:
                            self.tot[i] = self.tot[i] / self.ee * self.ef
                    for i in range(0, self.n+1):
                        if i != 11:
                            self.molal[i] = self.molal[i] / self.ee * self.ef
                    self.ee = self.ef

                if verbose:
                    print("iteration TDS")
                self.stdi = self.std

        if self.units == "molal" and self.dil == 0:
            self.cat[1:6] = self.tot[1:6]
            self.nchcat[1:6] = self.nch[1:6]
            self.ani[1:4] = self.tot[6:9]
            self.nchani[1:4] = -self.nch[6:9]
            self.ani[4] = self.molal[12]
            self.ani[5] = self.molal[14] + self.molal[16] + self.molal[17]
            self.nchani[4] = -self.nch[12]
            self.nchani[5] = -self.nch[14]

            self.density()

    def iterate_pco2(self, verbose):
        # calculate pco2 and ph
        self.calculate_pCO2(verbose=verbose)

        # if user specifies a pco2, redo the calculations
        if self.pco2_S == "y":

            while np.abs(self.po - self.pco2) > 0.01:
                if self.po < self.pco2 and self.poa < self.pco2:
                    self.ph -= self.dph
                if self.po < self.pco2 and self.poa > self.pco2:
                    self.dph = self.dph / 2
                    self.ph -= self.dph
                if self.po > self.pco2 and self.poa > self.pco2:
                    self.ph += self.dph
                if self.po > self.pco2 and self.poa < self.pco2:
                    self.dph = self.dph/2
                    self.ph += self.dph
                self.poa = self.po

                self.iterate_activities(verbose=verbose)

                self.calculate_pCO2(verbose=verbose)

        # if pk is greater than pkf,
        if self.pk > self.pkf:
            self.pk = self.pkf

            if verbose:
                print("last iteration")

            self.iterate_activities(verbose=verbose)
            self.calculate_pCO2(verbose=verbose)

    def calculate_pCO2(self, verbose=True):

        self.po = np.log10(self.act[15] / self.psc[14])

        if self.pco2_S == "":
            if self.diltot == 1:
                self.poinit = self.po
                self.phinit = self.ph

            if verbose:
                print("LOG PCO2 = {}".format(self.po))

            if self.pco2:
                self.po0 = self.po
                self.ph0 = self.ph
                self.pco2_S = "y"

            elif not self.pc:
                self.pco2_S = "n"

        if self.pco2_S == "y":
            if verbose:
                    print("\nLog(PCO2) selected = {}".format(self.pco2))
                    print("Log(PCO2) calculated = {}\n".format(self.po))

    def calculate_alkalinities(self):
        self.alcar = self.molal[12] + 2 * (self.molal[14] + self.molal[16] + self.molal[17])
        self.albor = self.molal[19] + self.molal[20] + 2 * self.molal[21] + self.molal[22] + self.molal[23]
        self.alsil = self.molal[24]
        self.aloh = self.molal[13] + self.molal[18]
        self.alh = -1 * self.molal[11] - self.molal[25]
        self.altest = self.alcar + self.albor + self.alsil + self.aloh + self.alh

    def print_screen(self, verbose=True, output=False):
        for i in range(1, self.n0+1):
            self.totinit[i] = 0
            for j in range(1, self.n+1):
                self.totinit[i] += self.molal[j] * self.kmat[i, j]

        df1 = pd.DataFrame({"molality": self.molal[1:self.n0+1],
                            "act coeff": self.gact[1:self.n0+1],
                            "activity": self.act[1:self.n0+1],
                            "molal tot": self.totinit[1:self.n0+1]},
                           index =  self.aq_S[1:self.n0+1])
        df1 = df1.loc[df1["molality"] != 0]
        df2 = pd.DataFrame({"molality": self.molal[self.n0+1:self.n+1],
                            "act coeff": self.gact[self.n0+1:self.n+1],
                            "activity": self.act[self.n0+1:self.n+1]},
                           index = self.aq_S[self.n0+1:self.n+1])
        df2 = df2.loc[df2["molality"] != 0]

        df = pd.concat([df1, df2])


        if verbose:
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None):
                print(df.to_string(na_rep=""))
                print()

        if output:
            with open(self.log_file, "a") as file:
                file.write(df.to_string(na_rep=""))
                file.write("\n\n")
                file.close()

        text = ["ELECTRICAL BALANCE     = {} % corrected on {} and {}".format(self.dca, self.aq_S[self.icat], self.aq_S[self.iani]),
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
            with open(self.log_file, "a") as file:
                for line in text:
                    file.write(line)
                    file.write("\n")
                file.write("\n")
                file.close()

        condition = np.where(self.ica==1)[0][np.where(np.where(self.ica==1)[0] <= 12)]

        u = np.zeros(condition.shape[0])

        for k in range(0, u.shape[0]):
            i = condition[k] * np.ones(self.n, dtype=int)
            j = np.arange(1, self.n+1)
            u[k] = np.sum(self.molal[j] * self.kmat[i, j])

        tests = self.tot[condition]
        d = np.zeros(condition.shape[0])

        for i in range(0, d.shape[0]):
            if u[i] + tests[i] != 0:
                d[i] = 200 * np.abs(u[i] + tests[i]) / (u[i] + tests[i])
        #d[u + tests != 0] = 200 * np.abs(u - tests) / (u + tests)



        names = np.char.upper(self.aq_S[condition].astype(str))
        names[-1] = "ALK"

        df = pd.DataFrame({"TESTS": names, "SUM OF SPECIES": u,
                           "INIT CONC.": tests, "BALANCE %": d})

        if verbose:
            with pd.option_context('display.max_rows', None,
                                               'display.max_columns', None):
                print(df.to_string(index=False))
            print()

        if output:
            with open(self.log_file, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()

        data_lst = []

        for i in range(13, self.n+1):
            u = 0
            if self.ica[i] == 1:
                for j in range(0, self.n+1):
                    if self.act[j] != 0:
                        u += self.kmat[i, j] * np.log10(self.act[j])
                v = np.log10(self.psc[i-12])
                d = 200 * np.abs(u - v) / (u + v)
                if i == 13:
                    zone_S = "h2o"
                else:
                    zone_S = self.aq_S[i].lower()
                data_lst.append([zone_S, u, v, d])

        df = pd.DataFrame(data_lst, columns=["", "log(IAP)",
                                             "log(K)", "balance %"])

        if verbose:
            with pd.option_context('display.max_rows', None,
                                               'display.max_columns', None):
                print(df.to_string(index=False))
            print()

        if output:
            with open(self.log_file, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()

        self.nwm = 0
        self.nwmp = 0

        self.pai = np.zeros(self.nm+1)
        data_lst = []

        for k in range(1, self.nm+1):
            self.pai[k] = 1
            for i in range(0, self.nt+1):
                self.pai[k] = self.pai[k] * self.act[i] ** self.wmin[k, i]
            if self.pai[k] != 0:
                if self.nwmin[k] == 0:
                    zone_S = " " + self.mineral_S[k].lower()
                elif self.nwmin[k] == 1:
                    zone_S = "*" + self.mineral_S[k].lower()
                x_S = " "
                if self.pai[k] / self.psol[k] >= 1 and self.nwmin[k] == 1:
                    self.nwm += 1
                    x_S = "*"
                elif self.pai[k] / self.psol[k] >= 0.9 and self.pai[k] / self.psol[k] < 1 and self.nwmin[k] == 1:
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
            with open(self.log_file, "a") as file:
                file.write(df.to_string(index=False))
                file.write("\n\n")
                file.close()

        if verbose and output:
            print("LOG FILE IS {}".format(self.log_file))
            print()

    def saturation_state(self, verbose):
        self.kinvar = np.zeros(self.nt+1)

        for k in range(1, self.nm+1):
            self.lmin[k] = 0

        if verbose:
            print("The initial solution is oversaturated in {} mineral(s)" \
                  " of the data base MURTF3:".format(self.nwm))
            print()

            oversat = []
            for k in range(1, self.nm+1):
                if self.pai[k] / self.psol[k] >= 1 and self.nwmin[k] == 1:
                    oversat.append([self.mineral_S[k], self.pai[k] / self.psol[k]])
                    self.lmin[k] = 1
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
                    print(self.mineral_S[self.kinvar[k]])
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
                for k in range(1, self.nm+1):
                    if (self.pai[k] / self.psol[k] >= 0.9 and
                        self.pai[k] / self.psol[k] < 1 and self.nwmin[k] == 1):

                        close_saturation.append([self.mineral_S[k], self.pai[k] / self.psol[k]])
                        self.lmin[k] = 1

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
                        print(self.mineral_S[self.kinvar[k]])
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

    def dilute_solution(self, dilute):
        pass

    def modify_database(self, add_min, rem_min, verbose, output):
        minerals = self.mineral_S.astype(str).tolist()
        if add_min:
            if type(add_min) is not list:
                if type(add_min) is str:
                    add_min = [add_min]
                else:
                    add_min = list(add_min)
            add_min = [i.upper() for i in add_min]
            am = [minerals.index(add_min[i]) for i in
                  range(0, len(add_min)) if add_min[i] in minerals]
            self.nwmin[am] = 1
            if verbose:
                print("ADDING MINERALS")
                for i in add_min:
                    if i in minerals:
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
            rm = [minerals.index(rem_min[i]) for i in
                  range(0, len(rem_min)) if rem_min[i] in minerals]
            self.nwmin[rm] = 0
            if verbose:
                print("REMOVING MINERALS")
                for i in rem_min:
                    if i in minerals:
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
                           self.mineral_S[self.nwmin == 1].tolist()]
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

    def invar(self):
        self.kinvariant = 0
        ncm = 14
        nbmin = 10
        ninvar = nbmin + 3
        kinv = np.zeros(ninvar+1)
        minv_S = np.empty(ninvar+1, dtype=object)
        psminv = np.zeros(ninvar+1)
        winv = np.zeros((ninvar+1, ncm+1))
        minvar_S = np.empty(ninvar+1, dtype=object)
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
        for k in range(1, self.nm+1):
            if self.lmin[k] == 1:
                n1 += 1
                kinv[n1] = k
                minv_S[n1] = self.mineral_S[k]
                psminv[n1] = np.log10(self.psol[k])
                for j in range(0, ncm+1):
                    winv[n1, j] = self.wmin[k, j]
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
                if np.abs(t1[i, k-1]) > self.epsilon:
                    u = t1[k-1, k-1] / t1[i, k-1]
                    for j in range(k, n1+1):
                        t1[i, j] = t1[k-1, j] - t1[i, j] * u
                        if np.abs(t1[i, j]) < self.epsilon:
                            t1[i, j] = 0
        det = 1
        for i in range(1, n1+1):
            if np.abs(t1[i, i]) < self.epsilon:
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
                        if np.abs(t2[i, k-1]) > self.epsilon:
                            u = t2[k-1, k-1] / t2[i, k-1]
                            for j in range(k, n2+1):
                                t2[i, j] = t2[k-1, j] - t2[i, j] * u
                                if np.abs(t2[i, j]) < self.epsilon:
                                    t2[i, j] = 0
                det1 = 1
                for i in range(1, n2+1):
                    if np.abs(t2[i, i]) < self.epsilon:
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
                    if u < self.epsilon:
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
                        if np.abs(t4[i, k-1]) > self.epsilon:
                            u = t4[k-1, k-1] / t4[i, k-1]
                            for j in range(k, n4+1):
                                t4[i, j] = t4[k-1, j] - t4[i, j] * u
                                if np.abs(t4[i, j]) < self.epsilon:
                                    t4[i, j] = 0
                            tt4[i] = tt4[k-1] - tt4[i] * u
                if np.abs(t4[n4, n4]) > self.epsilon:
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
                elif np.abs(t4[n4, n4]) <= self.epsilon:
                    self.kinvariant = -2

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


    def run_evp(self, verbose, output):

        if verbose:
            print("\nThis is EVP..............\n")
            print("STARTING THE EVAPORATION PROGRAM\n")
            print(datetime.now().strftime("%a %b %d %H:%M:%S %Y"), '\n')
            print()
        
        self.initialize_evp(verbose=verbose)
        
        self.loop_500(verbose, output)


    def initialize_evp(self, verbose):
        # temporary parameters for testing
        self.output_step = 1
        self.print_step = 1
        self.increment = 0


        # Set parameters for simulation
        
        self.ncpt = 0
        self.mwev = 0
        self.fc = 1
        self.q0_S = ""
        self.n = 25
        self.ntot = 12
        self.ncomplex = 14
        self.mh2o = 55.51
        
        self.ksupprim = 0
        self.ndepact = 0
        self.npasi = 0
        self.npasf = 0
        self.initdeseq = 0
        self.kinvariant = 0
        self.ninv = 0
        self.xinv  = 0
        
        self.my_S = ""
        self.my0_S = ""
        
        # Initialize blank arrays
        self.totinit = np.zeros(self.ntot+1)
        self.tot0 = np.zeros(self.ntot+1)
        self.totest = np.zeros(self.ntot+1)
        self.psc = np.zeros(self.ncomplex+1)

        # Initialize blank arrays
        self.nch = np.zeros(self.n+1, dtype=int)
        self.mol = np.zeros(self.n+1)
        self.mol0 = np.zeros(self.n+1)
        self.mol1 = np.zeros(self.n+1)
        self.molal0 = np.zeros(self.n+1)
        self.act = np.zeros(self.n+1)
        self.act0 = np.zeros(self.n+1)
        self.gact0 = np.zeros(self.n+1)
        self.gact1 = np.zeros(self.n+1)
        self.aq_S = np.empty(self.n+1, object)
        self.atom = np.zeros(self.n+1)
        self.kmat = np.zeros((self.n+1, self.n+1))

        # Read aquv database
        aquv = read_file("aquv.dat")
        self.aq_S[:] = [line[0].lstrip() for line in aquv]
        self.atom[:] = [line[1] for line in aquv]
        self.nch[:] = [line[2] for line in aquv]

        # Read kmat from matrice2
        self.kmat[1:,1:] = [line[:] for line in read_file("matrice2")]

        # transfer output from EQL into input for EVP
        self.totinit[1:11] = self.tot[1:11]
        self.tot[:] = 0
        
        self.nbmin = np.count_nonzero(self.totinit != 0)
        self.ica = np.zeros(self.n+self.nbmin+1, dtype=int)
        self.kinvar = np.zeros(self.nbmin+4, dtype=int)
        
        self.mol[0] = self.molal[15]
        self.mol[1:11] = self.molal[1:11]
        self.mol[11] = self.mh2o
        self.mol[12] = self.molal[13]
        self.mol[13] = self.molal[11]
        self.mol[14] = self.molal[14]
        self.mol[15] = self.molal[12]
        self.mol[16:26] = self.molal[16:26]
        
        self.mol0[:] = self.mol[:]
        self.molal[:] = 0

        # Set the increment mode (automatic or manual)
        if self.increment == 0:
            self.inc_S = "auto"
        else:
            self.inc_S = "manu"

        # Set the initial increment to the input increment
        self.inc0 = self.increment

        # Turn on all components by setting them to 1
        self.ica[1:self.n+1] = 1
        
        # Write the heading into the chemistry file
        with open(self.chem_file, "a") as file:
            file.write(",".join(self.constit_S))
            file.write("\n")
            file.close()

        # Write the initial events into the events file
        with open(self.event_file, 'a') as file:
            file.write("Temperature of solution = {} Deg C".format(self.tinit))
            file.write("     ")
            file.write("Temperature of simulation = {} Deg C".format(self.temp))
            file.write("\n")
            if self.diltot > 1:
                file.write("The initial solution has "\
                                  "been diluted {} times".format(self.diltot))
                file.write("\n")
            if self.ph != self.phinit:
                file.write("Initial Log(pco2) = {}     ".format(self.poinit))
                file.write("Selected Log(pco2) = {}".format(self.po))
                file.write("\n")
                file.write("Initial pH = {}     ".format(self.phinit))
                file.write("Calculated pH = {}".format(self.ph))
                file.write("\n")
            file.close()

        # BELOW THIS LINE IS OLD CODE FROM EQL METHOD. MAKE SURE TO UPDATE
        # Read thermodynamic data for dissociation coefficients from complex3
        complex3 = np.zeros((self.ncomplex+1, 5))
        complex3[1:,:] = [i[1:] for i in read_file("complex3")]
        self.psc[1:] = 10 ** (complex3[1:,0] +
                              complex3[1:,1] / 300 * self.temp +
                              complex3[1:,2] / 30000 * self.temp ** 2 +
                              complex3[1:,3] / 3000000 * self.temp ** 3 +
                              complex3[1:,4] / 300000000 * self.temp ** 4)
        self.psc3 = self.psc[3]
        self.psc[3] = self.psc[3] * self.psc[14] * 10 ** self.po




        # Read eql mineral database from murtf2/murtf0
        murtf = read_file(self.min_S)

        (self.nc, self.na, self.nm) = (int(i) for i in murtf[0])
        self.ncm = self.nc + self.na + 1

        self.wmin = np.zeros((self.nm+1, self.ncm+1))
        self.mu = np.zeros(self.ncm+1)
        self.linvar = np.zeros(self.nm+1, dtype=int)

        self.mineral_S = np.empty(self.nm+1, object)
        self.mum = np.zeros(self.nm+1)
        self.psol = np.zeros(self.nm+1)
        self.psol0 = np.zeros(self.nm+1)
        self.pai = np.zeros(self.nm+1)
        self.pai0 = np.zeros(self.nm+1)

        self.lmin = np.zeros(self.nm+1, dtype=int)
        self.lmin0 = np.zeros(self.nm+1, dtype=int)
        self.lmin1 = np.zeros(self.nm+1, dtype=int)

        self.min = np.zeros(self.nm+1)
        self.min0 = np.zeros(self.nm+1)
        self.minp = np.zeros(self.nm+1)
        self.minp0 = np.zeros(self.nm+1)

        ion_list = murtf[1:self.ncm+1]
        ion_list.insert(0, None)
        #ion_S = [i[0].lower() for i in ion_list]

        for i in range(1, self.ncm+1):
            a_S = str(ion_list[i][0]).lower()
            (at, bt, ct, dt, et) = (float(k) for k in ion_list[i][1:])
            for j in range(1, self.ncm+1):
                if a_S == str(self.aq_S[j]):
                    self.mu[j] = (at +
                                  bt / 300 * self.temp +
                                  ct / 30000 * self.temp ** 2 +
                                  dt / 3000000 * self.temp ** 3 +
                                  et / 300000000 * self.temp ** 4)

        mineral_list = murtf[self.ncm+1:]
        mineral_list.insert(0, None)
        self.mineral_S = np.empty(self.nm+1, object)

        for k in range(1, self.nm+1):
            line = mineral_list[k]
            self.mineral_S[k] = line[0]
            ncomp = int(line[1])
            c_ion = np.zeros(ncomp+1)
            nom_ion_S = np.empty(ncomp+1, object)
            for i in range(1, ncomp+1):
                c_ion[i] = float(line[i*2])
                nom_ion_S[i] = str(line[1+i*2])
                x_S = nom_ion_S[i].lower()
                for j in range(0, self.ncm+1):
                    if x_S == self.aq_S[j]:
                        self.wmin[k, j] = c_ion[i]
            (at, bt, ct, dt, et) = (float(i) for i in line[2+ncomp*2:])
            self.mum[k] = (at +
                           bt / 300 * self.temp +
                           ct / 30000 * self.temp ** 2 +
                           dt / 3000000 * self.temp ** 3 +
                           et / 300000000 * self.temp ** 4)

        for k in range(1, self.nm+1):
            u = self.mum[k]
            for i in range(0, self.ncm+1):
                u -= self.wmin[k, i] * self.mu[i]
            self.psol[k] = np.exp(u)
            self.psol0[k] = self.psol[k]

        with open(self.min_file, "w") as file:
            file.write("fc,")
            file.write(",".join(self.mineral_S[1:].tolist()))
            file.write("\n")
            file.close()

        # Charge Balance
        sc, cmax = 0, 0
        for i in range(1, self.n+1):
            if self.nch[i] > 0:
                sc += self.mol[i] * self.nch[i]
                if self.mol[i] * self.nch[i] > cmax:
                    cmax = self.mol[i] * self.nch[i]
                    ic = i
        sa, amax = 0, 0
        for i in range(1, self.n+1):
            if self.nch[i] < 0:
                sa += self.mol[i] * -self.nch[i]
                if self.mol[i] * -self.nch[i] > amax:
                    amax = self.mol[i] * -self.nch[i]
                    ia = i

        self.dca = 200 * np.abs(sc-sa) / (sc+sa)
        delta = sc-sa

        self.mol[ic] = self.mol[ic] - delta / 2 / self.nch[ic]
        self.mol[ia] = self.mol[ia] + delta / 2 / -self.nch[ia]

        sc, sa = 0, 0
        for i in range(1, self.n+1):
            if self.nch[i] > 0:
                sc += self.mol[i] * self.nch[i]
            if self.nch[i] < 0:
                sa += self.mol[i] * -self.nch[i]

        if verbose:
            print("sum of cations = {}".format(sc))
            print("sum of anions = {}".format(sa))
            print()

        for i in range(1, 12):
            self.totinit[i] = 0
            for j in range(1, self.n+1):
                self.totinit[i] += self.kmat[i, j] * self.mol[j]
                self.tot[i] = self.totinit[i]
                self.tot0[i] = self.totinit[i]

        self.tot[12] = 0
        self.ctot0 = (self.mol[0] + self.mol[14] + self.mol[15] +
                      self.mol[16] + self.mol[17])
        
        self.evp_actp()
        
        self.molal[:] = self.mol * self.mh2o / self.mol[11]
        self.act[:] = self.molal * self.gact[i]
        
        for k in range(1, self.nm+1):
            self.pai[k] = 1
            for i in range(1, self.ncm+1):
                self.pai[k] *= self.act[i] ** self.wmin[k, i]
            if self.pai[k] >= self.psol[k]:
                self.lmin[k] = 1
        
        with open(self.event_file, "a") as file:
            for k in range(1, self.nm+1):
                if self.lmin[k] == 1:
                    file.write("Initial solution "\
                               "oversaturated in {}".format(self.mineral_S[k]))
                    file.write("\n")
            file.close()
    
        for k in range(1, self.nm+1):
            if self.lmin[k] == 1:
                self.psol[k] = self.pai[k]
                
    def loop_500(self, verbose, output):
        self.ncmpt = 0
        ix, iy = 1, 2
        self.initdeseq = 0
        if self.mwev == 0:
            for k in range(1, self.nm+1):
                if self.psol[k] > self.psol0[k]:
                    self.initdeseq = 1
            if self.initdeseq == 1:
                for k in range(1, self.nm+1):
                    if self.psol[k] * 0.95 > self.psol0[k]:
                        self.psol[k] *= 0.95
                        if verbose:
                            print(self.mineral_S[k], self.psol[k], self.psol0[k])
                    elif self.psol[k] * 0.95 <= self.psol0[k]:
                        self.psol[k] = self.psol0[k]
        
        self.nw = 1
        while self.nw != 0:
            
            self.ncmpt += 1
            self.m0_S = "_".join(self.mineral_S[self.lmin == 1])
            
            self.gact1 = self.gact
            self.evp_actp()
            
            if self.kinvariant == 0:
                self.gact = ((self.gact + self.gact1 * ix) / iy)
            
            self.molal = self.mol * self.mh2o / self.mol[11]
            self.act = self.molal * self.gact
            
            for k in range(1, self.nm+1):
                self.pai[k] = 1
                for i in range(1, self.ncm+1):
                    self.pai[k] *= self.act[i] ** self.wmin[k, i]
                
                if self.pai[k] >= self.psol[k]:
                    if self.min[k] >= 0:
                        self.lmin[k] = 1
                    elif self.min[k] < 0:
                        self.lmin[k] = 0
                        self.min[k] = 0
                elif self.pai[k] < self.psol[k]:
                    if self.min[k] <= 0:
                        self.lmin[k] = 0
                        self.min[k] = 0
                    elif self.min[k] > 0:
                        self.lmin[k] = 1
            
            for k in range(1, self.nm+1):
                if self.psol[k] == 1e+50:
                    if self.pai[k] < self.psol0[k] * 0.9:
                        self.linvar[k] = 0
                    elif self.pai[k] >= self.psol0[k]:
                        self.linvar = 1
            
            self.nminer = np.count_nonzero(self.lmin == 1)
            self.mineraux_S = "_".join(self.mineral_S[self.lmin == 1])
            
            if (self.ncpt == 1) or (self.ncpt % self.output_step == 0):
                if self.nminer == 0:
                    print(self.ncmpt, "No_minerals")
                else:
                    print(self.ncmpt, self.mineraux_S)
            
            if ((self.mwev > 0) and (self.fc != 1) and 
                (self.nminer - self.nminer0 >= 2)):
                
                self.increment /= 2
                
                if self.increment < self.epsilon:
                    if verbose:
                        print()
                        print("Program unstable")
                        print("Restart the initialization "\
                              "program (EQL...)")
                        print("and lower the limits of convergence")
                    
                    if output:
                        with open(self.event_file, 'a') as file:
                            file.write("Program unstable\n")
                            file.write("Restart the initialization "\
                                       "program (EQL...)\n")
                            file.write("and lower the limits "\
                                       "of convergence\n")
                            file.close()
                
                        self.compact()
                
                    self.stop_simulation()
                    
                    return
                    
                if verbose:
                    print("reduction at increment {}".format(self.increment))
                
                self.mol = self.mol0
                self.tot = self.tot0
                self.lmin[1:] = self.lmin0[1:]
                self.min[1:] = self.min0[1:]
                self.minp[1:] = self.minp0[1:]
                
                self.mwev = self.mwev0
                self.nminer = self.nminer0
                self.mwev += self.mol[11] * self.increment
                self.tot[11] -= 2 * self.mol[11] * self.increment
                
                # start the loop over, exit the parent call when done
                self.loop_500(verbose, output)
                return
            
            if (self.nminer > 1) and (self.mineraux_S != self.m0_S):
                ix, iy = 2, 3
                
                self.evp_invar(verbose, output)
                
                if self.kinvariant > 0:
                    if self.system == "o":
                        for i in range(1, self.kinvariant+1):
                            if ((self.minp[self.kinvar[i]] == 0) and
                                (self.min[self.kinvar[i]] == 0)):
                                    
                                self.kneuf = self.kinvar[i]
                        
                        self.loop_2000()
                        return
                
                    elif self.system == "c":
                        self.mol = self.mol0
                        self.molal = self.molal0
                        self.gact = self.gact0
                        self.act = self.act0
                    
                        self.tot[1:] = self.tot0[1:]
                        self.pai[1:] = self.pai0[1:]
                        self.min[1:] = self.min0[1:]
                        
                        self.mwev = self.mwev0
                        self.nminer = self.nminer0
                        self.fi = self.fi0
                        
            self.mol1 = self.mol
            
            if self.kinvariant == 0:
                reseq = self.reseq(verbose, output) 
                if reseq == 999:
                    return
                
            elif self.kinvariant > 0:
                self.reseqinv()
                self.mwev += self.xinv/2
                self.tot[11] -= self.xinv
            
            self.mol[0] = (self.mol[15] * self.gact[15] * self.mol[13] * 
                           self.gact[13] / self.mol[11] / self.gact[11] / 
                           self.psc3 / self.gact[0])
        
            self.nw = 0
            
            for i in range(1, self.n+1):
                if self.mol[i] > 0:
                    if (200 * np.abs(self.mol[i] - self.mol1[i]) / 
                        (self.mol[i] + self.mol1[i]) > self.pkmol):
                        
                        self.nw = 1
                        
            self.ki = self.kinvariant
            
            if self.kinvariant > 0:
                for k in range(1, self.kinvariant+1):
                    if self.min[self.kinvar[k]] <= 0:
                        self.lmin[self.kinvar[k]] = 0
                        self.min[self.kinvar[k]] = 0
                        self.psol[self.kinvar[k]] = 1e+50
                        self.mwev += self.mol[11] * self.increment
                        self.tot[11] -= 2 * self.mol[11] * self.increment
                        self.ki, self.nw = 0, 1
                        
            self.kinvariant = self.ki
            
            if self.nw == 1:
                self.mol = (self.mol + self.mol1) / 2
            
            if self.ncmpt == 500:
                if verbose:
                    print()
                    print("Program unstable")
                    print("Restart the initialization program (EQL...)")
                    print("and lower the limits of convergence.")
                    print("Set the increment in manual mode at a lower "\
                          "value than .5")
                if output:
                    with open(self.event_file, 'a') as file:
                        file.write("Program unstable\n")
                        file.write("Restart the initialization "\
                                   "program (EQL...)\n")
                        file.write("and lower the limits of "\
                                   "convergence.\n")
                        file.write("Set the increment in manual mode "\
                                   "at a lower value than .5\n")
                        file.close()
                
                    self.compact()
                
                self.stop_simulation()
                return
            
        for k in range(1, self.nm+1):
            if self.psol[k] == 1e+50 and self.linvar[k] == 0:
                self.psol[k] = self.psol0[k]
                if verbose:
                    print("resetting: {}".format(self.mineral_S[k]))
            self.linvar[k] = 0
        
        self.npasi += 1
        self.npasf += 1
        self.ncpt += 1
        
        if self.system == "o":
            self.minp[1:] += self.min[1:]
            
            for i in range(1, 11):
                for k in range(1, self.nm+1):
                    self.tot[i] -= self.wmin[k, i] * self.min[k]
            
            for j in range(1, self.ncm+1):
                for k in range(1, self.nm+1):
                    self.tot[11] -= (self.wmin[k, j] * 
                                     self.kmat[11, j] * 
                                     self.min[k])
            
        for i in range(1, self.ntot):
            self.totest[i] = 0
            for j in range(1, self.n+1):
                self.totest[i] += self.kmat[i, j] * self.mol[j]
        
        self.totinit[12] = 0
        self.totest[12] = 0
        for j in range(1, self.n+1):
            if self.kmat[12, j] > 0:
                self.totest[12] += self.kmat[12, j] * self.mol[j]
            elif self.kmat[12, j] < 0:
                self.totinit[12] += self.kmat[12, j] * self.mol[j]
        self.totinit[12] = -self.totinit[12]
        
        for i in range(1, 11):
            for k in range(1, self.nm+1):
                if self.system == "c":
                    self.totest[i] += self.min[k] * self.wmin[k, i]
                if self.system == "o":
                    self.totest[i] += self.minp[k] * self.wmin[k, i]
        
        for j in range(1, self.ncm+1):
            for k in range(1, self.nm+1):
                if self.system == "c":
                    self.totest[11] += (self.wmin[k, j] * 
                                        self.kmat[11, j] * 
                                        self.min[k])
                if self.system == "o":
                    self.totest[11] += (self.wmin[k, j] * 
                                        self.kmat[11, j] * 
                                        self.minp[k])
                    
        self.totest[11] += self.mwev * 2
        
        self.evp_density()
        
        self.fc = self.mh2o / self.mol[11]
        self.alc = (self.molal[12] - self.molal[13] + 
                    self.molal[14] * 2 + self.molal[15] + 
                    (self.molal[16] + self.molal[17]) * 2)
        self.alc += (self.molal[18] + self.molal[19] + 
                     self.molal[20] + self.molal[21] * 2)
        self.alc += (self.molal[22] + self.molal[23] + 
                     self.molal[24] - self.molal[25])
        
        self.std = (np.sum(self.molal[1:] * self.atom[1:]) - 
                    self.molal[11] * self.atom[11])
        
        self.ee = 1000000e0 * self.dens / (1000e0 + self.std)
        
        self.hco3, self.co3 = 0, 0
        for k in range(1, self.nm+1):
            if self.system == "c":
                self.hco3 += self.wmin[k, 15] * self.min[k]
                self.co3 += self.wmin[k, 14] * self.min[k]
            elif self.system == "o":
                self.hco3 += self.wmin[k, 15] * self.minp[k]
                self.co3 += self.wmin[k, 14] * self.minp[k]
        
        self.ctot = (self.mol[0] + self.mol[14] + self.mol[15] + 
                     self.mol[16] + self.mol[17] + 
                     self.hco3 + self.co3)
        
        self.my_S = "_".join(self.mineral_S[(self.lmin == 1) | 
                                            (self.min != 0)])
        
        self.loop_600(verbose, output)
        return
                
                    
    
    def loop_600(self, verbose, output):
        
        if self.ncpt == 1 or self.ncpt % self.output_step == 0:
            
            min_name = []
            moles_prec = []
            moles_1 = []
            moles_tot = []
            tests = []
            
            if self.nminer > 0:
                for i in range(1, self.nm+1):
                    if self.lmin[i] == 1 or self.min[i] != 0:
                        
                        u = (200 * np.abs(self.pai[i] - self.psol[i]) / 
                             (self.pai[i] + self.psol[i]))
                        
                        if self.system == "c":
                            if self.min[i] > self.min0[i]:
                                x_S = "{} {}".format(self.mineral_S[i], "(p)")
                            elif self.min[i] < self.min0[i]:
                                x_S = "{} {}".format(self.mineral_S[i], "(d)")
                            elif self.min[i] == self.min0[i]:
                                x_S = "{} {}".format(self.mineral_S[i], "(=)")
                                
                            min_name.append(x_S)
                            moles_prec.append(self.min[i])
                            tests.append(u)
                            
                        else:
                            x_S = self.mineral_S[i]
                            min_name.append(x_S)
                            moles_1.append(self.min[i])
                            moles_tot.append(self.minp[i])
                            tests.append(u)
                            
                if self.system == "c":
                    db1 = pd.DataFrame(data={" ": min_name,
                                             "MOLES PREC": moles_prec,
                                             "TESTS": u})
                else:
                    db1 = pd.DataFrame(data={" ": min_name,
                                             "MOLES 1 STEP": moles_1,
                                             "MOLES TOT": moles_tot,
                                             "TESTS": tests})
                    
                if verbose:
                    with pd.option_context('display.max_rows', None,
                                           'display.max_columns', None):
                        print(db1.to_string(index=False))
                    print()
                    
                    
                if output:
                    with open(self.log_file, 'a') as file:
                        file.write("\n")
                        with pd.option_context('display.max_rows', None,
                                           'display.max_columns', None):
                            file.write(db1.to_string(index=False))
                        file.write("\n")
                        file.close()
            if self.my_S == "": 
                if verbose:
                    print("No_minerals")
                    
                if output:
                    with open(self.log_file, 'a') as file:
                        file.write("No_minerals")
                        file.write("\n")
                        file.close()
            
            species = []
            moles =[]
            molalities = []
            act_coef = []
            molal_tot = []
            tests = []
            
            for i in range(1, self.ntot+1):
                if self.tot[i] > 0 or i == 12:
                    u = (200 * np.abs(self.totest[i] - self.totinit[i]) / 
                         (self.totest[i] + self.totinit[i]))
                    
                    species.append(self.aq_S[i])
                    moles.append(self.mol[i])
                    molalities.append(self.molal[i])
                    act_coef.append(self.gact[i])
                    tests.append(u)
                    
                    if i <= 10:
                        s = 0
                        for j in range(1, self.n+1):
                            s += self.molal[j] * self.kmat[i, j]
                        molal_tot.append(s)
                    else:
                        molal_tot.append(np.nan)
            
            for i in range(self.ntot+1, self.n+1):
                if self.mol[i] > 0:
                    p = 1
                    for j in range(1, self.n+1):
                        p *= self.act[j] ** self.kmat[i, j]
                    u = (200 * np.abs(p - self.psc[i-12]) / 
                         (p + self.psc[i-12]))
                    species.append(self.aq_S[i])
                    moles.append(self.mol[i])
                    molalities.append(self.molal[i])
                    act_coef.append(self.gact[i])
                    tests.append(u)
                    molal_tot.append(np.nan)
            
            self.pco = np.log10(self.act[15] * self.act[13] / self.act[11] / 
                                self.psc3 / self.psc[14])
            u = (200 * np.abs(self.pco - self.po) / (self.pco + self.po))
            
            species.append(self.aq_S[0])
            moles.append(self.mol[0])
            molalities.append(self.molal[0])
            act_coef.append(self.gact[0])
            tests.append(u)
            molal_tot.append(np.nan)
            
            db2 = pd.DataFrame(data={" ": species,
                                     "MOLES": moles,
                                     "MOLALITIES": molalities,
                                     "ACT COEFF": act_coef,
                                     "MOLAL TOT": molal_tot,
                                     "TESTS": tests})
            
            if verbose:
                    with pd.option_context('display.max_rows', None,
                                           'display.max_columns', None):
                        print(db2.to_string(index=False, na_rep=""))
                    print()
                    
                    
            if output:
                with open(self.log_file, 'a') as file:
                    file.write("\n")
                    with pd.option_context('display.max_rows', None,
                                       'display.max_columns', None):
                        file.write(db2.to_string(index=False, na_rep=""))
                    file.write("\n")
                    file.close()
                    
            lines = []
            lines.append("concentration factor = {}".format(self.fc))
            lines.append("ionic strength       = {}".format(self.fi))
            lines.append("salinity (g/kg)      = {}".format(self.std))
            lines.append("activity of water    = {}".format(self.act[11]))
            lines.append("water evapor. (mol)  = {}".format(self.mwev))
            lines.append("pH                   = {}".format(-np.log10(self.act[13])))
            lines.append("CO2 exchanged (mol)  = {}".format(self.ctot - self.ctot0))
            lines.append("alkalinity (eq/kg)   = {}".format(self.alc))
            lines.append("Log PCO2             = {}".format(self.pco))
            lines.append("number of steps      = {}".format(self.ncpt))
            lines.append("molal/molar factor   = {}".format(self.ee / 1000))
            if self.kinvariant == 0:
                lines.append("increment (%)        = {}".format(self.increment * 100))
            else:
                lines.append("increment (moles)    = {}".format(self.xinv))
            lines.append("density              = {}".format(self.dens))
            
            if verbose:
                for line in lines:
                    print(line)
                print()
            
            if output:
                with open(self.log_file, 'a') as file:
                    for line in lines:
                        file.write(line)
                        file.write("\n")
                    file.write("\n")
                    file.close()
            
        if output:
            lines = []
            for k in range(1, self.nm+1):
                if self.system == "c":
                    if self.lmin[k] == 1 and self.lmin0[k] == 0:
                        lines.append("start of precipitation of {} at fc = {}".format(self.mineral_S[k], self.fc))
                    elif self.lmin[k] == 1 and self.lmin1[k] == 0:
                        if self.min[k] < self.min0[k]:
                            self.lmin1[k] = 1
                            lines.append("end of precipitation and start of precipitation of {} at fc = {}".format(self.mineral_S[k], self.fc))
                    
                    elif self.lmin[k] == 1 and self.lmin1[k] == 1:
                        if self.min[k] > self.min0[k]:
                            self.lmin1[k] = 0
                            lines.append("end of dissolution and of precipitation of {} at fc = {}".format(self.mineral_S[k], self.fc))
                    
                    elif (self.lmin[k] == 0 and 
                          self.lmin1[k] == 1 and 
                          self.lmin0[k] == 1):
                        
                        self.lmin[k] = 0
                        lines.append("end of dissolution and of saturation of {} at fc = {}".format(self.mineral_S[k], self.fc))
                    
                    elif (self.lmin[k] == 0 and
                          self.lmin1[k] == 0 and 
                          self.lmin0[k] == 1):
                        
                        self.lmin1[k] = 0
                        lines.append("end of saturation of {} at fc = {}".format(self.mineral_S[k], self.fc))
                elif self.system == "o":
                    if self.lmin[k] == 1 and self.lmin0[k] == 0:
                        lines.append("start of precipitation of {} at fc = {}".format(self.mineral_S[k], 
                                                      self.fc))
                    elif self.lmin[k] == 0 and self.lmin0[k] == 1:
                        lines.append("end of precipitation of {} at fc = {}".format(self.mineral_S[k], self.fc))
            with open(self.event_file, 'a') as file:
                for line in lines:
                    file.write(line)
                    file.write("\n")
                file.write("\n")
                file.close()
            
            if (self.ncpt == 1 or 
                self.my_S != self.my0_S or 
                self.npasf == self.output_step):
                
                if self.units == "molal":
                    self.ee = 1000
                if self.units == "molar":
                    self.ee = 1000000 * self.dens / (1000 + self.std)
                    
                line = []
                if self.my_S == "":
                    line.append("No_minerals")
                else:
                    line.append(self.my_S)
                line.append(self.fc)
                line.append(self.mwev)
                line.append(self.dens)
                line.append(-np.log10(self.act[13]))
                line.append(self.alc * self.ee)
                for i in range(1, 11):
                    if self.tot[i] > 0:
                        s = 0
                        for j in range(1, self.n+1):
                            s += self.molal[j] * self.kmat[i, j]
                        line.append(s * self.ee)
                line.append(self.std * self.ee)
                
                
                with open(self.chem_file, 'a') as file:
                    file.write(",".join([str(i) for i in line]))
                    file.write("\n")
                    file.close()
                
                line = []
                line.append(self.fc)
                if self.system == "c":
                    for i in range(1, self.nm+1):
                        if self.min[i] >= 0:
                            line.append(self.min[i])
                        elif self.min[i] < 0:
                            line.append(0)
                elif self.system == "o":
                    for i in range(1, self.nm+1):
                        if self.minp[i] >= 0:
                            line.append(self.minp[i])
                        elif self.minp[i] < 0:
                            line.append(0)
                
                with open(self.min_file, 'a') as file:
                    file.write(",".join([str(i) for i in line]))
                    file.write("\n")
                    file.close()
                
                
        if self.stdmax > 0 and self.std * self.ee >= self.stdmax:
            if output:
                self.compact()
            self.stop_simulation()
            return
        
        if (self.mwev > 0 and 
            self.kinvariant == 0 and 
            (np.abs(self.act[11] - self.act0[11]) / 
             (self.act[11] + self.act0[11])) < .0000000001):
            
            if self.system == "c":
                nu = 0
                for k in range(1, self.nm+1):
                    if self.min[k] > 0 and self.min[k] < self.min0[k]:
                        nu = 1
                if nu == 0:
                    if self.nminer == self.nbmin:
                        self.q_S = "invariant system / eutectic point / end"
                    if self.nminer < self.nbmin:
                        self.q_S = "invariant system / pseudo-eutectic point / end"
                elif nu == 1:
                    if self.nperitec == 0:
                        self.peritec()
                    if self.nperitec == 0:
                        if self.nminer == self.nbmin:
                            self.q_S = "invariant system / peritectic point / end"
                        if self.nminer < self.nbmin:
                            self.q_S = "invariant system / pseudo-peritectic point / end"
                    elif self.nperitec == 1:
                        if self.nminer == self.nbmin:
                            self.q_S = "invariant system / peritectic / passing over"
                        if self.nminer < self.nbmin:
                            self.q_S = "invariant system / pseudo-peritectic / passing over"
            if self.system == "o":
                if self.nminer == self.nbmin:
                    self.q_S = "invariant system / eutectic / end"
                if self.nminer < self.nbmin:
                    self.q_S = "invariant system / pseudo-eutectic / end"
            
            if "pseudo" in self.q_S:
                self.q1_S = ("maximum number of minerals allowed by the "\
                    "phase rule = {}".format(self.nbmin))
                self.q2_S = ("number of minerals in equilibrium with the "\
                    "invariant system = {}".format(self.nminer))
                self.q_S = "\n".join([self.q_S, self.q1_S, self.q2_S])
            
            self.q0_S = self.q_S
            
            if verbose:
                print(self.q_S)
            
            if output:
                with open(self.log_file, 'a') as file:
                    file.write(self.q_S)
                    file.write("\n")
                    file.close()
                    
            if "end" in self.q_S:
                if output:
                    self.compact()
                self.stop_simulation()
                return
            
            # Line 1004
        elif (self.kinvariant == 0 and
              (np.abs(self.act[11] - self.act0[11]) / 
               (self.act[11] + self.act0[11])) > .0000000001):
            
            self.nperitec = 0
            self.q_S = ""
            self.q0_S = ""
            self.q1_S = ""
            self.q2_S = ""
        
        if self.npasi == self.print_step:
            self.npasi = 0
        if self.npasf == self.output_step:
            self.npasf = 0
            
        if self.my_S != self.my0_S:
            self.my0_S = self.my_S
            
        self.mol0 = self.mol
        self.molal0 = self.molal
        self.gact0 = self.gact
        self.act0 = self.act
        
        self.tot0[1:] = self.tot[1:]
        
        self.lmin0[1:] = self.lmin[1:]
        self.pai0[1:] = self.pai[1:]
        self.min0[1:] = self.min[1:]
        self.minp0[1:] = self.minp[1:]
        
        self.fi0 = self.fi
        self.mwev0 = self.mwev
        self.nminer0 = self.nminer
        
        if self.kinvariant == 0 and self.initdeseq == 0:
            if self.inc_S == "auto":
                self.increment = (51 - 8 * np.log10(self.std * 1000)) / 700
                self.inc0 = self.increment
            elif self.inc_S == "manu":
                self.increment = self.inc0
            self.mwev += self.mol[11] * self.increment
            self.tot[11] -= 2 * self.mol[11] * self.increment
        
        self.loop_500(verbose, output)
        return


            
    
    def loop_2000(self, verbose, output):
        
        if verbose:
            print("Free energy minimization")
        
        if self.kinvariant == 2:
            if (self.wmin[self.kinvar[1], self.ncm] > 
                self.wmin[self.kinvar[2], self.ncm]):
                
                self.ksupprim = self.kinvar[1]
            else:
                self.ksupprim = self.kinvar[2]
        elif self.kinvariant > 2:
            self.gmin = 0
            for ii in range(1, self.kinvariant+1):
                self.mol = self.mol0
                self.tot[1:] = self.tot0[1:]
                
                self.mwev += self.mol[11] * self.increment
                self.tot[11] -= 2 * self.mol[11] * self.increment
                
                self.lmin[1:] = self.lmin0[1:]
                self.min[1:] = self.min0[1:]
                self.minp[1:] = self.minp0[1:]
                
                if (self.kinvar[ii] > 0 and 
                    self.kinvar[ii] != self.kneuf and 
                    self.min[self.kinvar[ii]] >= 0):
                    
                    l = self.lmin[self.kinvar[ii]]
                    m = self.min[self.kinvar[ii]]
                    p = self.psol[self.kinvar[ii]]
                    
                    self.lmin[self.kinvar[ii]] = 0
                    self.min[self.kinvar[ii]] = 0
                    self.psol[self.kinvar[ii]] = 1e+50
                    
                    if verbose:
                        print("mineral removed: {}"\
                              "".format(self.mineral_S[self.kinvar[ii]]))
                    
                    self.ncmptinv = 0
                    self.nw = 1
                    while self.nw != 0:
                        self.ncmptiv += 1
                        
                        self.gact0 = self.gact
                        self.evp_actp()
                        self.gact = (self.gact0 + self.gact) / 2
                        self.molal = self.mol * self.mh2o / self.mol[11]
                        self.act = self.molal * self.gact
                        
                        for k in range(1, self.nm+1):
                            self.pai[k] = 1
                            for i in range(1, self.ncm+1):
                                self.pai[k] *= self.act[i] ** self.wmin[k, i]
                            
                            if self.pai[k] >= self.psol[k]:
                                if self.min[k] >= 0:
                                    self.lmin[k] = 1
                                elif self.min[k] < 0:
                                    self.lmin[k] = 0
                                    self.min[k] = 0
                            elif self.pai[k] < self.psol[k]:
                                if self.min[k] <= 0:
                                    self.lmin[k] = 0
                                    self.min[k] = 0
                                elif self.min[k] > 0:
                                    self.lmin[k] = 1

                        self.nminer = np.count_nonzero(self.lmin == 1)
                        self.mineraux_S = "_".join(self.mineral_S[self.lmin == 1])
                    
                    if verbose:
                        print(self.ncmptinv, self.mineraux_S)
                        
                    self.mol1 = self.mol
                    
                    self.reseq()
                    
                    self.molal = self.mol * self.mh2o / self.mol[11]
                    
                    self.nw = 0
                    for i in range(0, self.n+1):
                        if np.abs(self.mol1[i] - self.mol[i]) > self.pkmol:
                            self.nw = 1
                    
                    g = 0
                    for i in range(0, self.ncm+1):
                        if self.mol[i] > 0:
                            g += self.mol[i] * self.mu[i] + np.log(self.act[i])
                    for k in range(1, self.nm+1):
                        g += self.min[k] * self.mum[k]
                    
                    if verbose:
                        print("g = {}".format(g))
                    
                        for i in range(1, self.kinvariant):
                            if i != ii:
                                print(self.mineral_S[self.kinvar[i]])
                        
                        print()
                    
                    if g < self.gmin:
                        self.gmin = g
                        self.ksupprim = self.kinvar[ii]
                    
                    self.lmin[self.kinvar[ii]] = l
                    self.min[self.kinvar[ii]] = m
                    self.psol[self.kinvar[ii]] = p
                    
        self.mol = self.mol0
        self.tot[1:] = self.tot0[1:]
        self.lmin[1:] = self.lmin0[1:]
        self.min[1:] = self.min0[1:]
        self.minp[1:] = self.minp0[1:]
        
        self.mwev = self.mwev0
        self.nminer = self.nminer0
        self.mwev += self.mol[11] * self.increment
        self.tot[11] -= 2 * self.mol[11] * self.increment
        self.kinvariant = 0
        self.lmin[self.ksupprim] = 0
        self.min[self.ksupprim] = 0
        self.psol[self.ksupprim] = 1e+50
        self.min0[self.ksupprim] = 0
        
        if verbose:
            print("mineral definitely removed: {}"\
                  "".format(self.mineral_S[self.ksupprim]))
        
        self.loop_500(verbose, output)
        return
    
    def compact(self):
        df = pd.read_csv(self.min_file)
        df = df.loc[:, (df != 0).any(axis=0)]
        df.to_csv(self.min_file, index=False)
    
    def stop_simulation(self):
        pass
    
    def reseq(self, verbose, output):
        
        z = np.zeros((self.n+self.nminer+1, self.n+self.nminer+1))
        zz = np.zeros(self.n+self.nminer+1)
        xx = np.zeros(self.n+self.nminer+1)
        
        nconv = 1
        ncm = 15
        nt = self.n
        
        nu = 1
        while nu != 0:
            for i in range(1, nt+1):
                for j in range(1, nt+1):
                    z[i, j] = 0
                zz[i] = 0
            for i in range(1, 13):
                for j in range(1, self.n+1):
                    if self.mol[j] != 0:
                        z[i, j] = self.kmat[i, j]
            
            for i in range(13, self.n+1):
                self.kmat[i, 0] = 0
                for j in range(1, self.n+1):
                    if j != 11:
                        self.kmat[i, 0] += self.kmat[i, j]
            
            for i in range(13, self.n+1):
                p = 1
                u = 0
                for j in range(1, self.n+1):
                    if self.mol[j] != 0 and j != 11:
                        z[i, j] = self.kmat[i, j] / self.mol[j]
                        p *= self.gact[j] ** self.kmat[i, j]
                        u += self.kmat[i, j] * np.log(self.mol[j])
                    elif j == 11:
                        z[i, j] = -self.kmat[i, 0] / self.mol[j]
                        p *= self.aw ** self.kmat[i, j]
                        u -= self.kmat[i, 0] * np.log(self.mol[j])
                    elif self.mol[j] == 0:
                        z[i, j] = 0
                p *= self.mh2o ** self.kmat[i, 0]
                zz[i] = np.log(self.psc[i-12]) - np.log(p) - u
            
            l = 0
            for k in range(1, self.nm+1):
                if self.lmin[k] == 1:
                    l += 1
                    self.ica[self.n+l] = 1
                    self.wmin[k, 0] = 0
                    p = 1
                    u = 0
                    for j in range(1, ncm+1):
                        if j != 11:
                            self.wmin[k, 0] += self.wmin[k, j]
                    for j in range(1, ncm+1):
                        if j != 11 and self.mol[j] > 0:
                            z[self.n+l, j] = self.wmin[k, j] / self.mol[j]
                            p *= self.gact[j] ** self.wmin[k, j]
                            u += self.wmin[k, j] * np.log(self.mol[j])
                        elif j == 11:
                            z[self.n+l, j] = -self.wmin[k, 0] / self.mol[j]
                            p *= self.aw ** self.wmin[k, j]
                            u -= self.wmin[k, 0] * np.log(self.mol[j])
                    p *= self.mh2o ** self.wmin[k, 0]
                    zz[self.n+l] = np.log(self.psol[k]) - np.log(p) - u
                    
                    sh = 0
                    for j in range(1, ncm+1):
                        sh += self.wmin[k, j] * self.kmat[11, j]
                    for i in range(1, 11):
                        z[i, self.n+l] = self.wmin[k, i]
                    z[11, self.n+l] = sh
            
            nt = self.n + l
            for i in range(1, 11):
                u = 0
                for j in range(1, self.n+1):
                    u += self.kmat[i, j] * self.mol[j]
                for k in range(1, self.nm+1):
                    if self.lmin[k] == 1:
                        u += self.min[k] * self.wmin[k, i]
                
                zz[i] = self.tot[i] - u
            
            u = 0
            for j in range(1, self.n+1):
                u += self.kmat[11, j] * self.mol[j]
            for j in range(1, ncm+1):
                for k in range(1, self.nm+1):
                    u += self.wmin[k, j] * self.kmat[11, j] * self.min[k]
            
            zz[11] = self.tot[11] - u
            
            u = 0
            for j in range(1, self.n+1):
                u += z[12, j] * self.mol[j]
            zz[12] = self.tot[12] - u
            
            for k in range(1, 11):
                if self.tot[k] == 0:
                    self.ica[k] = 0
                    for j in range(k+1, self.n+1):
                        if self.kmat[k, j] != 0:
                            self.ica[j] = 0
            
            ni = nt
            nj = nt
            for k in range(nt, 0, -1):
                if self.ica[k] == 0:
                    for i in range(k, ni):
                        for j in range(1, nj+1):
                            z[i, j] = z[i+1, j]
                        zz[i] = zz[i+1]
                    ni -= 1
                    for j in range(k, nj):
                        for i in range(1, ni+1):
                            z[i, j] = z[i, j+1]
                    nj -= 1
            
            for k in range(1, ni+1):
                for i in range(k, ni+1):
                    if np.abs(z[i, k-1]) > self.epsilon:
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
            
            for k in range(1, nt+1):
                if self.ica[k] == 0:
                    for i in range(ni, k-1, -1):
                        xx[i+1] = xx[i]
                    xx[k] = 0
                    ni += 1
            
            nu = 1
            while nu != 0:
                for i in range(1, self.n+1):
                    if self.ica[i] == 1:
                        nu = 0
                        if self.mol[i] + xx[i] / nconv < 0:
                            nconv += 1
                            nu = 1
                            break
                if nconv >= self.max_conv:
                    if verbose:
                        print()
                        print("The equation system diverges: end of "\
                              "simulation")
                    if output:
                        with open(self.event_file, 'a') as file:
                            file.write("\n")
                            file.write("The equation system diverges: "\
                                       "end of simulation\n")
                            file.close()
                        self.compact()
                    
                    self.stop_simulation()
                    
                    return(999)
                
            for i in range(1, self.n+1):
                self.mol[i] += xx[i] / nconv
            i = self.n
            for k in range(1, self.nm+1):
                if self.lmin[k] == 1:
                    i += 1
                    self.min[k] += xx[i] / nconv
            nu = 0
            for i in range(1, nt+1):
                if self.ica[i] == 1:
                    if np.abs(xx[i]) > self.pkeq:
                        nu = 1
        return

                
    
    def reseqinv(self):
        
        nt = self.ntot-1
        
        t = np.zeros((nt+1, nt+1))
        tt = np.zeros((nt+1))
        t0 = np.zeros((nt+1, nt+1))
        tt0 = np.zeros((nt+1))
        xx = np.zeros((nt+1))
        
        if self.ninv == 0:
            hmin = 1000
            for k in range(1, self.nm+1):
                for kk in range(1, self.kinvariant+1):
                    if k == self.kinvar[kk] and self.min[k] > 0:
                        s = 0
                        for i in range(1, 16):
                            s += self.wmin[k, i] * self.kmat[nt, i]
                        if s > 0:
                            if s * self.min[k] < hmin:
                                hmin = s * self.min[k]
            self.xinv = hmin / 100
            if self.xinv <= 0:
                self.stop_simulation()
                return
        self.ninv = 1
        j = 0
        for k in range(1, self.nm+1):
            if self.lmin[k] == 1:
                for kk in range(1, self.kinvariant+1):
                    if k == self.kinvar[kk]:
                        j += 1
                        for i in range(1, nt):
                            t0[i, j] = self.wmin[k, i]
                        s = 0
                        for i in range(1, 16):
                            s += self.wmin[k, i] * self.kmat[nt, i]
                        t0[nt, j] = s
        nmin = j
        
        for i in range(1, nt):
            tt[i] = 0
            for k in range(1, self.nm+1):
                for kk in range(1, self.kinvariant+1):
                    if k == self.kinvar[kk]:
                        tt0[i] += self.min[k] * self.wmin[k, i]
        tt0[nt] = 0
        for k in range(1, self.nm+1):
            for kk in range(1, self.kinvariant+1):
                if k == self.kinvar[kk]:
                    for i in range(1, 16):
                        tt0[nt] += (self.min[k] * self.wmin[k, i] * 
                                    self.kmat[11, i])
        tt0[nt] -= self.xinv
        
        for i in range(1, nmin+1):
            for j in range(i, nmin+1):
                t[i, j] = 0
                for k in range(1, nt+1):
                    t[i, j] += t0[k, i] * t0[k, j]
                    t[j, i] = t[i, j]
        for i in range(1, nmin+1):
            tt[i] = 0
            for k in range(1, nt+1):
                tt[i] += t0[k, i] * tt0[k]
        
        for k in range(2, nmin+1):
            for i in range(k, nmin+1):
                if np.abs(t[i, k-1]) > self.epsilon:
                    u = t[k-1, k-1] / t[i, k-1]
                    for j in range(k, nmin+1):
                        t[i, j] = t[k-1, j] - t[i, j] * u
                    tt[i] = tt[k-1] - tt[i] * u
        
        xx[nmin] = tt[nmin] / t[nmin, nmin]
        for i in range(nmin-1, 0, -1):
            s = 0
            for j in range(i+1, nmin+1):
                s += t[i, j] * xx[j]
            xx[i] = (tt[i] - s) / t[i, i]
            if np.abs(xx[i]) < self.epsilon:
                xx[i] = 0
        i = 0
        for k in range(1, self.kinvariant+1):
            i += 1
            self.min[self.kinvar[k]] = xx[i]
            if self.min[self.kinvar[k]] <= 0:
                self.ninv = 0
        
        return
        
    def peritec(self):
        nt = self.ntot-1
        t = np.zeros((nt+1, nt+1))
        tt = np.zeros((nt+1))
        t0 = np.zeros((nt+1, nt+1))
        tt0 = np.zeros((nt+1))
        xx = np.zeros((nt+1))
        
        j = 0
        for k in range(1, self.nm+1):
            if self.lmin[k] == 1:
                j += 1
                for i in range(1, nt):
                    t0[i, j] = self.wmin[k, i]
                s = 0
                for i in range(1, self.ncm+1):
                    s += self.wmin[k, i] * self.kmat[11, i]
                t0[11, j] = s
        j += 1
        t0[nt, j] = 2
        nj = j
        for i in range(1, nt+1):
            tt0[i] = self.tot[i]
        
        for i in range(1, nj+1):
            for j in range(i, nj+1):
                t[i, j] = 0
                for k in range(1, nt+1):
                    t[i, j] += t0[k, i] * t0[k, j]
                    t[j, i] = t[i, j]
        for i in range(1, nj+1):
            tt[i] = 0
            for k in range(1, nt+1):
                tt[i] += t0[k, i] * tt0[k]
        
        for k in range(2, nj+1):
            for i in range(k, nj+1):
                if np.abs(t[i, k-1]) > self.epsilon:
                    u = t[k-1, k-1] / t[i, k-1]
                    for j in range(k, nj+1):
                        t[i, j] = t[k-1, j] - t[i, j] * u
                    tt[i] = tt[k-1] - tt[i] * u
        
        xx[nj] = tt[nj] / t[nj, nj]
        for i in range(nj-1, 0, -1):
            s = 0
            for j in range(i+1, nj+1):
                s += t[i, j] * xx[j]
            xx[i] = (tt[i] - s) / t[i, i]
        self.nperitec = 0
        for i in range(1, nj):
            if xx[i] < 0:
                self.nperitec = 1
        return
    
    def evp_invar(self, verbose, output):
        self.kinvariant = 0
        self.ncm = 15
        self.ninvar = self.nbmin + 3
        
        kinv = np.zeros((self.ninvar+1), dtype=int)
        minv_S = np.empty(self.ninvar+1, object)
        psminv = np.zeros((self.ninvar+1))
        winv = np.zeros((self.ninvar+1, self.ncm+1))
        minvar_S = np.empty(self.ninvar+1, object)
        psminvar = np.zeros((self.ninvar+1))
        t0 = np.zeros((self.ninvar+1, self.ninvar+1))
        t1 = np.zeros((self.ninvar+1, self.ninvar+1))
        t2 = np.zeros((self.ninvar+1, self.ninvar+1))
        t3 = np.zeros((self.ninvar+1, self.ncm+1))
        t4 = np.zeros((self.ncm+1, self.ncm+1))
        tt4 = np.zeros((self.ncm+1))
        
        
        psminv[1:4] = np.log10(self.psc[1:4])
        
        winv[1, 13] = 1
        winv[1, 12] = 1
        winv[1, 15] = -1
        winv[2, 13] = 1
        winv[2, 11] = -1
        winv[2, 14] = 1
        winv[3, 13] = 1
        winv[3, 11] = 1
        winv[3, 15] = -1
        
        n1 = 3
        for k in range(1, self.nm+1):
            if self.lmin[k] == 1:
                n1 += 1
                kinv[n1] = k
                minv_S[n1] = self.mineral_S[k]
                psminv[n1] = np.log10(self.psol[k])
                for j in range(1, self.ncm+1):
                    winv[n1, j] = self.wmin[k, j]
                winv[n1, 11], winv[n1, 15] = winv[n1, 15], winv[n1, 11]
        
        for i in range(1, n1+1):
            for j in range(i, n1+1):
                t1[i, j] = 0
                for k in range(1, self.ncm):
                    t1[i, j] += winv[i, k] * winv[j, k]
                    t1[j, i] = t1[i, j]
                    t0[i, j] = t1[i, j]
                    t0[j, i] = t0[i, j]
        for k in range(2, n1+1):
            for i in range(k, n1+1):
                if np.abs(t1[i, k-1]) > self.epsilon:
                    u = t1[k-1, k-1] / t1[i, k-1]
                    for j in range(k, n1+1):
                        t1[i, j] = t1[k-1, j] - t1[i, j] * u
                        if np.abs(t1[i, j]) < self.epsilon:
                            t1[i, j] = 0
        det = 1
        for i in range(1, n1+1):
            if np.abs(t1[i, i]) < self.epsilon:
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
                        if t2[i, k-1] != 0:
                            u = t2[k-1, k-1] / t2[i, k-1]
                            for j in range(k, n2+1):
                                t2[i, j] = t2[k-1, j] - t2[i, j] * u
                                if np.abs(t2[i, j]) < self.epsilon:
                                    t2[i, j] = 0
                det1 = 1
                for i in range(1, n2+1):
                    if np.abs(t2[i, i]) < self.epsilon:
                        det1 = 0
                        break
                if det1 == 1:
                    n3 += 1
                    self.kinvar[n3] = kinv[kk]
                    minvar_S[n3] = minv_S[kk]
                    psminvar[n3] = psminv[kk]
                    for j in range(1, self.ncm+1):
                        t3[n3, j] = winv[kk, j]
            n4 = self.ncm
            for j in range(self.ncm, 0, -1):
                u = 0
                for i in range(1, n3+1):
                    u += t3[i, j] ** 2
                if u == 0:
                    for k in range(j+1, n4+1):
                        for i in range(1, n3+1):
                            t3[i, k-1] = t3[i, k]
                    n4 -= 1
            
            for i in range(1, n4+1):
                for j in range(i, n4+1):
                    t4[i, j] = 0
                    for k in range(1, n4+1):
                        t4[i, j] += t3[k, i] * t3[k, j]
                        t4[j, i] = t4[i, j]
            for i in range(1, n4+1):
                tt4[i] = 0
                for k in range(1, n3+1):
                    tt4[i] += t3[k, i] * psminvar[k]
            
            
            for k in range(2, n4+1):
                for i in range(k, n4+1):
                    if np.abs(t4[i, k-1]) > self.epsilon:
                        u = t4[k-1, k-1] / t4[i, k-1]
                        for j in range(k, n4+1):
                            t4[i, j] = t4[k-1, j] - t4[i, j] * u
                            if np.abs(t4[i, j]) < self.epsilon:
                                t4[i, j] = 0
                        tt4[i] = tt4[k-1] - tt4[i] * u
            
            
            self.ah2o = 10 ** (tt4[n4] / t4[n4, n4])
            self.kinvariant = n3
            for i in range(self.kinvariant, 0, -1):
                if self.kinvar[i] == 0:
                    for k in range(1, self.kinvariant):
                        self.kinvar[k] = self.kinvar[k+1]
                    self.kinvariant -= 1
            if verbose:
                print("invariant system constrained by: ")
                for k in range(1, self.kinvariant+1):
                    print("    {}".format(self.mineral_S[self.kinvar[k]]))
                print()
                print("invariant aH2O = {}".format(self.ah2o))
                print("simulation aH2O = {}".format(self.aw))
            if output:
                with open(self.event_file, 'a') as file:
                    file.write("invariant system constrained by: \n")
                    for k in range(1, self.kinvariant+1):
                        file.write("    {}".format(self.mineral_S[self.kinvar[k]]))
                        file.write("\n")
                    file.write("\n")
                    file.write("invariant aH2O = {}\n".format(self.ah2o))
                    file.write("simulation aH2O = {}\n".format(self.aw))
                    file.close()
        return
        
        
    
    def evp_density(self):
        ncdens, nadens = 5, 5
        s = np.zeros((ncdens+1, nadens+1))
        cat, ic = np.zeros(ncdens+1), np.zeros(ncdens+1)
        ani, ia = np.zeros(nadens+1), np.zeros(nadens+1)
        
        for i in range(1, 9):
            if i <= ncdens:
                cat[i] = self.mol[i] / self.mol[11] * self.mh2o
            if i > ncdens:
                ani[i-ncdens] = self.mol[i] / self.mol[11] * self.mh2o
                
        ani[4] = self.mol[15] / self.mol[11] * self.mh2o
        ani[5] = self.mol[14] / self.mol[11] * self.mh2o
        
        for i in range(1, ncdens+1):
            ic[i] = self.nch[i]
        for i in range(1, 4):
            ia[i] = -self.nch[i+5]
        ia[4] = -self.nch[15]
        ia[5] = -self.nch[14]
        
        if self.ncpt == 1:
            (self.au, self.bu) = (np.zeros((ncdens+1, nadens+1)), 
                                    np.zeros((ncdens+1, nadens+1)))
    
            densite = read_file("densite")
            
            index = 0
            
            for i in range(1, 6):
                for j in range(1, 6):
                    (self.au[i, j], self.bu[i, j]) = densite[index][3:]
                    index += 1
            
        self.dens, u = 1, 0
        for j in range(1, nadens+1):
            u += ani[j] * ia[j]
        for i in range(1, ncdens+1):
            for j in range(1, nadens+1):
                s[i, j] = int((ic[i] + ia[j]) / 2) * cat[i] * ani[j] / u
                self.dens += (self.au[i, j] * s[i, j] + 
                              self.bu[i, j] * s[i, j] ** 2)
        return