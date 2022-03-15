# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:01:36 2020

@author: AlecY
"""

# import resources
import planets
from orbitPropagator import orbitPropagator,null_perturbations,Satellite
from plotting import plot_body,plot_trajs,plot_ECI_states,plot_Kep_states,plot_ECEF_states,plot_LatLon

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
    
    # # set initial conditions and propagate satellite 1 (Chief Sat)
    r0_mag = 621.9 + planets.earth['radius'] # km
    v0_mag = np.sqrt(planets.earth['mu']/r0_mag)
    state0 = [r0_mag,0,0,0,v0_mag,0]
    sat1 = Satellite(name="Sat 1",mass=100,cd=2.2,area=1e-6)
    sat1.prop = orbitPropagator(sat1,state0,'ECI',epoch,prop_time,dt,planets.earth,perts)
    sat1.prop.propagate()
    
    # set initial conditions and propagate satellite 2 (Deputy Sat)
    state0 = [42164,0.1,5,1,1,1] # a (km),e,i,w,LAN,M (angles in degrees)
    sat2 = Satellite(name="Sat 2",mass=100,cd=2.2,area=1e-6)
    sat2.prop = orbitPropagator(sat2,state0,'Kep',epoch,prop_time,dt,planets.earth,perts)
    sat2.prop.propagate()
    
    # populate list of propagators
    propagators = [sat1.prop,sat2.prop]
    
    ########################################################
    ##### PLOTTING
    print("Plotting the trajectory")
    plt.close('all')
    
    # initialize plot
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    plot_lim = 40000 # km
    ax.set_xlim([-plot_lim,plot_lim])
    ax.set_ylim([-plot_lim,plot_lim])
    ax.set_zlim([-plot_lim,plot_lim])

    # plot earth
    plot_body(planets.earth,fig,ax,bluemarble=False)

    # plot trajectory of states
    plot_trajs(propagators,fig,ax)
    
    # plot ECI states
    plot_ECI_states(propagators)
    
    # plot Kep states
    plot_Kep_states(propagators,units="days")
    
    # plot ECEF states
    plot_ECEF_states(propagators,units="days")
    
    # plot Lat Lon
    plot_LatLon(propagators)
    
    fig.show()