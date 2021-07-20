#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Will Arnuk
"""
import os
import re
import pandas as pd
import numpy as np

def file_detect(directory, filetype):
    
    """ This function iterates through all files in a directory, returning
    a list of files that match a specified ending. """
    
    files = []
    for file in os.listdir(directory):
        if file.endswith(filetype):
            #print(os.path.join("/mydir", file))
            files.append(file)
    return(files)
    
def get_details(infile, trial):
    
    """ This function retrieves the input values (from the input file) that
    were used to run a given simulation"""
    
    data = pd.read_csv(infile)
    return(data[data['label'] == trial])      
    
def parse_outfile(file, output_path):
    
    """ This function uses regular expressions to parse through the specified
    dumpfile for concentration factor, salinity, activity of H2O, and 
    ionic strength at each step in the simulation. It compiles and saves
    a dataframe of these values as a csv file to the specified output path."""
    
    with open(file, 'r') as f:
        file_contents = f.read()
    
    fc = re.findall('concentration factor\s*=\s*(\d*.\d*)\.*', file_contents)
    salinity = re.findall('salinity\s*=\s*(\d*\.\w*.\w*)\s.*', file_contents)
    aH2O = re.findall('activity of water\s*=\s*(\d*.\d*)\.*', file_contents)
    ionStr = re.findall('ionic strength\s*=\s*(\d*\.\w*.\w*)\s.*', file_contents)
    
    data = pd.DataFrame({'fc': fc,
                         'aH2O': aH2O,
                         'ionStr': ionStr,
                         'salinity': salinity})
    data.to_csv('{path}.csv'.format(path = output_path), index=False)

def reprocess(directory, filetype):
    
    """ This function parses through all files of specified filetype (% or &) 
    in the specified directory, joining the dumpfile attributes to each file.
    
    The dumpfile attributes should be in a csv file created by the 
    output parser function. """
    
    for i in file_detect(directory, filetype):
        main_file = r'{directory}/{file}'.format(directory=directory, file=i)
        join_file = main_file.split('.')[0]+'.csv'
    
        main_df = pd.read_csv(main_file)
        
        # remove whitespace characters from headings
        cols = []
        for col in main_df.columns:
            cols.append(col.strip())
        main_df.columns = cols
        
        try:
            join_df = pd.read_csv(join_file)
        except FileNotFoundError:
            print("File '{jf}' does not exist".format(jf=join_file))
    
        for variable in join_df:
            if variable not in main_df:
                main_df[variable] = join_df[variable]
        
        main_df.to_csv(main_file, index=False)

def generate_inputs(n, ds, ph, na, k, li, ca, mg, cl, so4,
                    alk, no3, si, b, pco2=0.001, t_range=[0, 50]):
    
    label = []
    for i in range(0,n):
        label.append('trial'+str(i+1))
    label = np.array(label)
    
    pH = np.ones(shape=n) * ph
    density = np.ones(shape=n) * ds
    
    Na = np.ones(shape=n) * na
    K = np.ones(shape=n) * k
    Ca = np.ones(shape=n) * ca
    Mg = np.ones(shape=n) * mg
    Cl = np.ones(shape=n) * cl
    alk = np.ones(shape=n) * alk
    SO4 = np.ones(shape=n) * so4
    B = np.ones(shape=n) * b
    Li = np.ones(shape=n) * li
    NO3 = np.ones(shape=n) * no3
    Si = np.ones(shape=n) * si
    
    temp = np.arange(t_range[0], t_range[1], (t_range[1]-t_range[0])/n)

    pCO2 = np.ones(shape=n) * pco2
    pCO2 = np.log10(pCO2)

    x = pd.DataFrame(
        {"label": label,
         "t": temp,
         "ds": density,
         "ph": pH,
         "na": Na,
         "k": K,
         "li": Li,
         "ca": Ca,
         "mg": Mg,
         "cl": Cl,
         "so4": SO4,
         "alk": alk,
         "no3": NO3,
         "si": Si,
         "b": B,
         "pco2": pCO2})

    x.to_csv("input.dat", index=False)