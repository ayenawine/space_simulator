# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:01:36 2020

@author: AlecY
"""

# import resources
import planets
from satellite import Satellite
from orbitPropagator import orbitPropagator,null_perturbations
from plotting import run_visualizer

from datetime import datetime
    
if __name__ == "__main__":
    
    # set simulation settings
    epoch = datetime.today()
    prop_time = 0.2 * (24 * 60 * 60.0) # seconds
    dt = 30.0 # seconds, time step size
    
    # initialize perturbations
    perts = null_perturbations()
    perts['J2'] = False
    perts['aero_drag'] = False
    
    # earth centric orbits
    central_body = planets.earth
    
    # set initial conditions and propagate sat
    state0 = [8000,0.01,58.0,0.0,0.0,0.0] # a (km),e,i,LAN,w,M (angles in degrees)
    sat1 = Satellite(name="Sat 1",mass=100,cd=2.2,area=1e-6)
    sat1.prop = orbitPropagator(sat1,state0,'Kep',epoch,prop_time,dt,central_body,perts)
    
    # set initial conditions and propagate sat
    state0 = [7000,0.01,58.0,0.0,0.0,0.0] # p (km),e,i,LAN,w,M (angles in degrees)
    sat2 = Satellite(name="Sat 2",mass=100,cd=2.2,area=1e-6)
    sat2.prop = orbitPropagator(sat2,state0,'Kep',epoch,prop_time,dt,central_body,perts)
    
    # add a transfer maneuver
    sat2.prop.add_maneuver(target=sat1,maneuver_type='rendezvous')
    
    sat1.prop.propagate()
    sat2.prop.propagate()
    
    # populate list of propagators
    propagators = [sat1.prop,sat2.prop]

    [viz_3d,viz_ground] = run_visualizer(propagators,central_body)