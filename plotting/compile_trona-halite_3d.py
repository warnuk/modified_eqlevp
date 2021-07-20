# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:54:59 2021

@author: warnu
"""

trona = compile_xyz('TRONA', base, 'closed', 'first', z_var)
halite = compile_xyz('HALITE', base, 'closed', 'first', z_var)
np.save('3d_arrays/trona', trona)
np.save('3d_arrays/halite', halite)

halite_open = compile_xyz('HALITE', base, 'open', 'first', z_var)
np.save('3d_arrays/halite_open', halite_open)
