# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 21:30:53 2020

@author: AlecY
"""

import numpy as np
import math

def calc_atm_density(z,cb):
    """

    Parameters
    ----------
    z : altitude above sea level of the specified planetary body

    Returns
    -------
    The atmospheric density as the provided altitude

    """
    
    # first, get bounding points
    zs,rhos = find_rho_z(z,cb)
    print("alt bounds:",zs)
    print("den bounds:",rhos)
    
    # if atmospheric density is outside data range, return value of zero
    if rhos[0] == 0:
        return 0.0
    
    Hi = -(zs[1]-zs[0])/math.log(rhos[1]/rhos[0])
    
    # perform exponential interpolation
    new_rho = rhos[0]*math.exp(-(z-zs[0])/Hi)
    
    return new_rho
    
def find_rho_z(z,cb):
    """
    find the bounding altitude and atmospheric density points for interpolation
    """
    
    # if z beyond the max altitude provided
    # assume the density is zero
    if not z<cb['atm'][-1][0]:
        return [[0.0,0.0],[0.0,0.0]]

    zs   = cb['atm'][:,0]    
    rhos = cb['atm'][:,1]
    
    # find bounding points around the given altitude for interpolation
    for n in range(len(rhos)-1):
        if zs[n] < z < zs[n+1]:
            return [[zs[n],zs[n+1]],[rhos[n],rhos[n+1]]]
        
    # if still no value returned, return zeros
    return [[0.0,0.0],[0.0,0.0]]