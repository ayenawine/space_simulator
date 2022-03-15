# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:01:36 2020

@author: AlecY
"""

# import resources
import planets
from orbitPropagator import orbitPropagator,null_perturbations,Satellite
from plotting import plot_body,plot_trajs,plot_ECI_states,plot_Kep_states

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
if __name__ == "__main__":
    
    # set simulation settings
    prop_time = 5 * (24 * 60 * 60.0) # seconds
    dt = 30.0 # seconds, time step size
    
    # initialize perturbations
    perts = null_perturbations()
    perts['J2'] = True
    perts['aero_drag'] = False
    
    # set initial conditions and propagate satellite 1
    r0_mag = 500 + planets.earth['radius'] # km
    v0_mag = np.sqrt(planets.earth['mu']/r0_mag)
    r0 = [r0_mag,0,0]
    v0 = [0,v0_mag,0]
    sat1 = Satellite(name="Sat 1",mass=100,cd=2.2,area=1e-6)
    satprop1 = orbitPropagator(sat1,r0,v0,prop_time,dt,planets.earth,perts)
    satprop1.propagate()
    
    # set initial conditions and propagate satellite 2
    r0_mag = 800 + planets.earth['radius'] # km
    v0_mag = np.sqrt(planets.earth['mu']/r0_mag)
    r0 = [r0_mag,0,0]
    v0 = [0,v0_mag,v0_mag*.1]
    sat2 = Satellite(name="Sat 2",mass=100,cd=2.2,area=1e-6)
    satprop2 = orbitPropagator(sat2,r0,v0,prop_time,dt,planets.earth,perts)
    satprop2.propagate()
    
    # populate list of propagators
    propagators = [satprop2]
    
    ########################################################
    ##### PLOTTING
    print("Plotting the trajectory")
    
    # initialize plot
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    plot_lim = 15000 # km
    ax.set_xlim([-plot_lim,plot_lim])
    ax.set_ylim([-plot_lim,plot_lim])
    ax.set_zlim([-plot_lim,plot_lim])

    # plot earth
    plot_body(planets.earth,fig,ax)

    # plot trajectory of states
    plot_trajs(propagators,fig,ax)
    
    # plot ECI states
    plot_ECI_states(propagators)
    
    # plot Kep states
    plot_Kep_states(propagators,units="days")
    
    fig.show()