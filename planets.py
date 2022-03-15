# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 19:54:17 2020

@author: AlecY

Dictionaries of planets and their basic properties
"""

import numpy as np

sun={"name":    "Sun",
     "mass":    1.989e30,  # kg
     "mu":      1.327e11,  # km3 / s
     "radius":  695700.0,
     "J2":      None} # km

earth={"name":    "Earth",
       "mass":    5.972e24,    # kg
       "mu":      398600.4,    # km3 / s
       "radius":  6378.0,      # km
       "J2":      1.08262668e-3,
       
       # atmospheric parameters
       "atm":     np.array([[  63.096,2.059e-4*(10**8)],    # altitude (km)
                            [ 251.189,5.909e-11*(10**8)],   # vs density (kg/km3)
                            [1000.000,3.561e-15*(10**8)]]),
       "v_atm":   np.array([0.0,0.0,72.9211e-6])    # rotation speed of atmosphere
                                                    # relative to earth (rad/s)
      }