# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:01:36 2020

@author: AlecY
"""

# import resources
import planets
from orbitPropagator import orbitPropagator,null_perturbations,Satellite
from plotting import plot_body,plot_trajs
from plotting import plot_ECI_states,plot_Kep_states,plot_ECEF_states,plot_LatLon

from matplotlib import animation
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
    
    # set initial conditions and propagate iss
    # https://www.heavens-above.com/orbit.aspx?satid=25544
    p   = 11067.790
    e   =     0.83285
    i   =    87.87
    LAN =   227.89
    w   =    53.38
    v   =    92.335
    state0 = [p,e,i,LAN,w,v] # p (km),e,i,LAN,w,M (angles in degrees)
    state0 = [7000,0.01,58.0,0.0,0.0,0.0] # p (km),e,i,LAN,w,M (angles in degrees)
    sat = Satellite(name="Sat 1",mass=100,cd=2.2,area=1e-6)
    sat.prop = orbitPropagator(sat,state0,'Kep',epoch,prop_time,dt,planets.earth,perts)
    sat.prop.propagate()
    
    state0 = [8000,0.01,58.0,0.0,0.0,0.0] # p (km),e,i,LAN,w,M (angles in degrees)
    sat2 = Satellite(name="Sat 2",mass=100,cd=2.2,area=1e-6)
    sat2.prop = orbitPropagator(sat2,state0,'Kep',epoch,prop_time,dt,planets.earth,perts)
    sat2.prop.propagate()
    
    # populate list of propagators
    propagators = [sat.prop,sat2.prop]
    
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
    animation = plot_trajs(propagators,fig,ax,animated=False)
    
    # plot ECI states
    plot_ECI_states(propagators)
    
    # plot Kep states
    plot_Kep_states(propagators,units="days")
    
    # plot ECEF states
    plot_ECEF_states(propagators,units="days")
    
    # plot Lat Lon
    plot_LatLon(propagators)
    
    fig.show()