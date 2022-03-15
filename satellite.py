# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 15:28:41 2022

@author: AlecY
"""

class Satellite:
    """
    name = satellite name
    mass = mass of satellite (kg)
    cd = coefficient of drag
    area = surface area normal to the atmospheric velocity
    """
    def __init__(self,name,mass,cd,area):
        self.name = name
        self.mass = mass # kg
        self.cd = cd
        self.area = area # km2
        self.prop = None