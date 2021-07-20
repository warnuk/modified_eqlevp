# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:12:42 2021

@author: warnu
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import color_mapper
import file_handler

plt.style.use("ggplot")

title = 'Owens (No Dolomite)'

base = r"processed_owens_no-dol"
pco2 = 1000

temperature = 30
n = temperature-1

colormap = "mineral_colors.csv"
variable = "fc"

open_minerals = ["NORTHUPITE", "TRONA", "HALITE", "SYLVITE"]





# defining which colors are used to represent each mineral
if colormap:
    colormap = color_mapper.color_read(colormap)

# make axes
fig, axs = plt.subplots(2, 1)

ax1 = axs[0, 0]
ax2 = axs[1,0]


co = "open"
extension = "j{x}".format(x=co[0])

chem_file = r"F://EQLEVP/{base}/{pco2}ppm_output/{co}/output/trial{n}.{ext}&".format(base=base, pco2=pco2, co=co, n=n, ext=extension)
min_file = r"F://EQLEVP/{base}/{pco2}ppm_output/{co}/output/trial{n}.{ext}%".format(base=base, pco2=pco2, co=co, n=n, ext=extension)



chemdata = pd.read_csv(chem_file)

ions = ['na', 'k', 'alk', 'ca', 'mg', 'cl']

for ion in ions:
    x = np.log10(chemdata['fc'].values)
    y = chemdata[ion].values
    ax1.plot(x, y, label=ion)

# reading the output from the specified trial
mindata = pd.read_csv(min_file)

cols = mindata.columns.values.tolist()
for i in ['fc', 'aH2O', 'ionStr', 'salinity']:
    cols.remove(i)
y = {}

for mineral in cols:
    y[mineral] = mindata[mineral].values

negligible = []

for mineral in y:
    if y[mineral].sum() == 0:
        negligible.append(mineral)

for mineral in negligible:
    y.pop(mineral)


x = np.log10(mindata[variable].values)
plt.xlabel('fc')

for mineral in y:
    if colormap and mineral in colormap:
        ax2.plot(x, y[mineral], color=colormap[mineral], label=mineral)
    else:
        ax2.plot(x, y[mineral], label=mineral)



#plt.set_title("{base} {pco2}ppm {t}°C {co}".format(base=base, pco2=pco2, t=temperature, co=co))


co = "closed"
extension = "j{x}".format(x=co[0])

chem_file = r"F://EQLEVP/{base}/{pco2}ppm_output/{co}/output/trial{n}.{ext}&".format(base=base, pco2=pco2, co=co, n=n, ext=extension)
min_file = r"F://EQLEVP/{base}/{pco2}ppm_output/{co}/output/trial{n}.{ext}%".format(base=base, pco2=pco2, co=co, n=n, ext=extension)


chemdata = pd.read_csv(chem_file)

ions = ['na', 'k', 'alk', 'ca', 'mg', 'cl']

for ion in ions:
    x = np.log10(chemdata['fc'].values)
    y = chemdata[ion].values
    ax3.plot(x, y, label=ion)

# reading the output from the specified trial
mindata = pd.read_csv(min_file)

cols = mindata.columns.values.tolist()
for i in ['fc', 'aH2O', 'ionStr', 'salinity']:
    cols.remove(i)
y = {}

for mineral in cols:
    y[mineral] = mindata[mineral].values

negligible = []

for mineral in y:
    if y[mineral].sum() == 0:
        negligible.append(mineral)

for mineral in negligible:
    y.pop(mineral)


x = np.log10(mindata[variable].values)
ax4.set_xlabel('fc')

for mineral in y:
    if colormap and mineral in colormap:
        ax4.plot(x, y[mineral], color=colormap[mineral], label=mineral)
    else:
        ax4.plot(x, y[mineral], label=mineral)

ax1.set_ylabel("mmol/kg H2O")
ax1.set_yscale('log')
ax1.legend()

ax2.set_yscale('log')
ax2.set_ylabel('mmol')
ax2.legend()

ax3.set_ylabel("mmol/kg H2O")
ax3.set_yscale('log')
ax3.legend()

ax4.set_yscale('log')
ax4.set_ylabel('mmol')
ax4.legend()

ax2.set_ylim([10**-8, 10**-1])
ax4.set_ylim([10**-8, 10**-1])

ax1.set_ylim([10**-3, 10**4])
ax3.set_ylim([10**-3, 10**4])

ax1.set_title("No Backreaction")
ax3.set_title("Backreaction")

plt.suptitle('{title}, {t}°C, {pco2}ppm pCO2'.format(title=title, t=temperature, pco2=pco2), fontsize = 20)
plt.show()