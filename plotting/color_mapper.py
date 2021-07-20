#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Will Arnuk
"""
import pandas as pd

# read colormap csv into dictionary
def color_read(color_file):
    colors = pd.read_csv(color_file)
    colors.mineral = colors.mineral.apply(str.upper)
    colors = colors.set_index('mineral').T.to_dict('list')
    for i in colors:
        colors[i] = colors[i][0]
    return(colors)