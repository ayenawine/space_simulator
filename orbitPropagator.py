# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 20:21:36 2020

@author: AlecY
"""

import numpy as np
from scipy.integrate import ode
from datetime import datetime

import planets
from utils import calc_atm_density
from state_conversions import KeptoECI,ECItoECEF,ECEFtoLatLonAlt,ECItoKep

def null_perturbations():
    perturbations = {'J2':          False,
                     'aero_drag':   False,
                     'lunar_grav':  False,
                     }
    return perturbations

class orbitPropagator:
    """
    sat = Satellite object class
    r0 = initial position vector (km)
    v0 = initial velocity vector (km/s)
    prop_time = propagation time (s)
    dt = delta time, time step (s)
    cb = central body (earth,sun,etc)
    name = propagator/satellite name
    perturbations = dictionary of perturbation forces
    """
    def __init__(self,sat,state0,state0type,epoch,prop_time,dt,cb=planets.earth,perturbations=None):
        self.sat = sat
        self.state0 = state0
        self.state0type = state0type
        self.prop_time = prop_time
        self.dt = dt
        self.cb = cb
        self.name = sat.name
        self.perturbations = perturbations
        self.epoch = epoch
        
        # get j2000 reference time in seconds
        j2000 = datetime(2000,1,1,0,0,0,0)
        self.epochJ2000 = (self.epoch-j2000).total_seconds()
        
        # initialize state arrays
        self.n_steps      = int(np.ceil(self.prop_time/self.dt))
        self.times        = np.zeros((self.n_steps,1)) # array of times
        self.times_J2000  = np.zeros((self.n_steps,1)) # array of times in J2000 reference frame
        
        self.statesECI       = np.zeros((self.n_steps,6)) # array of states for each time step
        self.statesKep       = np.zeros((self.n_steps,6)) # array of states in Kepler Orbital Elements
        self.statesECEF      = np.zeros((self.n_steps,6)) # array of states in ECEF
        self.statesLatLonAlt = np.zeros((self.n_steps,3)) # array of Lat Lons
    
    def propagate(self):
    
        print("Propagating {}".format(self.name))
        
        # set initial conditions in ECI (convert if provided in kepler elements)
        if self.state0type == 'Kep':
            # use Kep state to get first ECI state
            self.statesKep[0]  = np.array(self.state0)
            self.statesECI[0]  = KeptoECI(self.statesKep[0],self.cb['mu'])
        elif self.state0type == 'ECI':
            # use provided ECI state and also get first kepler state
            self.statesECI[0]  = np.array(self.state0)
            self.statesKep[0]  = ECItoKep(self.statesECI[0],self.cb['mu'])
        else:
            raise Exception("Provided initial state type does not match and supported types.")
        
        print("initial ECI state:",self.statesECI[0])
        print("initial Kep state:",self.statesKep[0])
        
        # get other initial states and move on
        self.step = 0
        self.statesECEF[self.step]      = ECItoECEF(self.statesECI[self.step],self.epoch,self.times[self.step])
        self.statesLatLonAlt[self.step] = ECEFtoLatLonAlt(self.statesECEF[self.step],self.cb['radius'])
        self.step = 1
    
        # initialize the ode solver
        self.solver = ode(self.diff_eq)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.statesECI[0],0) # initial state at time zero
    
        # propagate the orbit using the ode solver for each time step
        while self.step < self.n_steps and self.solver.successful():
            
            # solve next sate using the ode solver
            self.solver.integrate(self.solver.t+self.dt)
            self.times[self.step] = self.solver.t
            self.statesECI[self.step] = self.solver.y
            
            # get the other states
            self.statesKep[self.step]       = ECItoKep(self.statesECI[self.step],self.cb['mu'])
            self.statesECEF[self.step]      = ECItoECEF(self.statesECI[self.step],self.epoch,self.times[self.step])
            self.statesLatLonAlt[self.step] = ECEFtoLatLonAlt(self.statesECEF[self.step],self.cb['radius'])
            
            # debug
            # print("time: ",self.times[self.step])
            # print("ECI state:",self.statesECI[self.step])
            # print("Kep state:",self.statesKep[self.step])
            
            # proceed to next time step
            self.step += 1
        
    def diff_eq(self,t,state):
    
        # unpack state
        rx,ry,rz,vx,vy,vz = state
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])
    
        # calculate magnitude of the radius
        r_norm = np.linalg.norm(r)
    
        # diff_eq equation
        # calculate accelleration vector
        a = -r*self.cb['mu'] / r_norm**3 # equation 2.15
    
        # add J2 acceleration
        if self.perturbations['J2']:
            z2 = r[2]**2
            r2 = r_norm**2
            tx = r[0]/r_norm*(5*z2/r2-1)
            ty = r[1]/r_norm*(5*z2/r2-1)
            tz = r[2]/r_norm*(5*z2/r2-3)
            
            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/r_norm**4*np.array([tx,ty,tz])
    
            a += a_j2
        
        # add aerodynamic drag
        if self.perturbations['aero_drag']:
            # calculate altitude and air density
            z = r_norm - self.cb['radius']    # altitude
            rho = calc_atm_density(z,self.cb) # atmospheric density at altitude
            
            print("alt:",z,"density:",rho)
            
            # calculate velocity of spacecraft relative to the rotating atmosphere
            v_rel = v - np.cross(self.cb['v_atm'],r)
            print("v_sat:",v)
            print("v_air:",np.cross(self.cb['v_atm'],r))
            print("ang vel:",self.cb['v_atm'])
            
            a_aero = -0.5*((self.sat.cd*self.sat.area)/self.sat.mass)*rho*v_rel*np.linalg.norm(v_rel)
            
            print("a_aero:",a_aero)
            print("a:     ",a)
            
            a += a_aero
            
        ax, ay, az = a
        return [vx,vy,vz,ax,ay,az]
    
    def add_maneuver(target=None,maneuver_type=None):
        print("add maneuver in progess")        