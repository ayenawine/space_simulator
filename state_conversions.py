# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 11:46:30 2022

@author: AlecY
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
from datetime import datetime
import spiceypy as spice # get pck here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
# import kernal data
spice.furnsh("C:\\Users\\AlecY\\Documents\\python\\trajectory_optimization\\resources\\solar_system_kernal.bpc")
j2000 = datetime(2000,1,1,0,0,0,0)

def KeptoECI(stateKep,mu):
    # source: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
    
    # unpack state
    # ==================================
    # p = semiparameter (km)
    # e = eccentricity
    # i = inclination
    # w = argument of periapsis (deg)
    # LAN = longitude of ascending nodes
    # M = mean anomaly (deg)
    # u = argument of latitude
    # lam_true = true longitude
    [p,e,i,LAN,w,v] = stateKep
    v   = np.deg2rad(v)
    LAN = np.deg2rad(LAN)
    i   = np.deg2rad(i)
    w   = np.deg2rad(w)
    
    # check edge cases
    if e <= 0.001 and i <= 5.0:
        w = 0.0
        LAN = 0.0
    elif e <= 0.001:
        w = 0.0
    elif e > 0.001:
        LAN = 0.0
    
    # NOTE: using true anomaly input instead to follow book
    # # solve for the eccentric anomaly using keplers equation & Newton–Raphson method
    # diff = 1.0
    # Eold = np.deg2rad(M)
    # while diff > 0.000001:
    #     E = Eold - (Eold-e*np.sin(Eold)-np.deg2rad(M)) / (1 - e*np.cos(Eold))
    #     diff = abs(E - Eold)
    #     Eold = E
    
    # # additional calcs: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors
    # nu = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E/2))
    # h = np.sqrt(self.cb['mu'] * a * (1-e**2))
    
    # # solve for true anomaly
    # x = np.sqrt(1+e)*np.sin(E/2)
    # y = np.sqrt(1-e)*np.cos(E/2)
    # v = np.arctan2(x,y) # true anomaly
    
    # # calculate radial distance to central body
    # r = a * ( 1-e*np.cos(E) )
    
    # # get position and velocity vector in the orbital frame
    # o_pos = r * np.array([np.cos(v),np.sin(v),0]);
    # o_vel = (np.sqrt(self.cb['mu']*a)/r) * np.array([-np.sin(E),np.sqrt(1-e**2)*np.cos(E),0]);
    
    r_PQW = np.array([ (p*np.cos(v))/(1+e*np.cos(v)) , (p*np.sin(v))/(1+e*np.cos(v)) , 0 ])
    v_PQW = np.array([ -np.sin(v) , e+np.cos(v) , 0 ]) * np.sqrt(mu/p)
    
    # rotate to the intertial body centric frame (using rot matrices)
    # Rz1 = R.from_euler('z', -LAN, degrees=True).as_matrix()
    # Rx  = R.from_euler('x', -i,   degrees=True).as_matrix()
    # Rz2 = R.from_euler('z', -w,   degrees=True).as_matrix()
    # rot = Rz1*Rx*Rz2
    X11 =  np.cos(LAN)*np.cos(w)-np.sin(LAN)*np.sin(w)*np.cos(i)
    X12 = -np.cos(LAN)*np.sin(w)-np.sin(LAN)*np.cos(w)*np.cos(i)
    X13 =  np.sin(LAN)*np.sin(i)
    X21 =  np.sin(LAN)*np.cos(w)+np.cos(LAN)*np.sin(w)*np.cos(i)
    X22 = -np.sin(LAN)*np.sin(w)+np.cos(LAN)*np.cos(w)*np.cos(i)
    X23 = -np.cos(LAN)*np.sin(i)
    X31 =  np.sin(w)*np.sin(i)
    X32 =  np.cos(w)*np.sin(i)
    X33 =  np.cos(i)
    rot =  np.array([[X11, X12, X13],
                    [X21, X22, X23],
                    [X31, X32, X33],
                   ])
    r_pos = np.dot(rot,r_PQW) # matrix multiplication is necessary here for some reason
    r_vel = np.dot(rot,v_PQW)
    
    # join the ECI state
    ECIstate = np.concatenate((r_pos,r_vel),axis=0)
    
    # get ECI state using spice instead for now
    # https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/conics_c.html
    '''
    RP      Perifocal distance.
    ECC     Eccentricity.
    INC     Inclination.
    LNODE   Longitude of the ascending node.
    ARGP    Argument of periapse.
    M0      Mean anomaly at epoch.
    T0      Epoch.
    MU      Gravitational parameter.
    '''
    #stateKep2 = [a*(1-e),e,i,w,LAN,M,self.epochJ2000,self.cb['mu']]
    #ECIstate = spice.conics(stateKep2,self.epochJ2000)
    
    return ECIstate

def ECItoKep(stateECI,mu):
    # source: http://control.asu.edu/Classes/MAE462/462Lecture07.pdf
    
    # unpack state
    rx,ry,rz,vx,vy,vz = stateECI
    
    r_vec = np.array([rx,ry,rz])
    r_norm = np.linalg.norm(r_vec)
    
    v_vec = np.array([vx,vy,vz])
    v_norm = np.linalg.norm(v_vec)

    h_vec = np.cross(r_vec,v_vec)
    h_norm = np.linalg.norm(h_vec)
    
    # normal vector pointing towards the right ascention
    n_vec = np.cross(np.array([0,0,1]),h_vec)
    n_norm = np.linalg.norm(n_vec)
    
    # eccentricity vector
    e_vec = np.cross(v_vec,h_vec)/mu - r_vec/r_norm
    e_norm = np.linalg.norm(e_vec)
    
    # true anomaly
    if np.inner(r_vec,v_vec) >= 0:
        v = np.arccos(np.inner(e_vec,r_vec)/(e_norm*r_norm))
    else:
        v = 2*np.pi - np.arccos(np.inner(e_vec,r_vec)/(e_norm*r_norm))
    
    # calculate orbital elements
    e = e_norm                                          # eccentricity
    E = 2*np.arctan(np.tan(v/2)/np.sqrt((1+e)/(1-e)))   # eccentric anomaly
    i = np.arccos(h_vec[2]/h_norm)                      # inclination
    
    # longitude of ascending node
    if n_vec[1] >= 0:
        Omega = np.arccos(n_vec[0]/n_norm)
    else:
        Omega = 2*np.pi - np.arccos(n_vec[0]/n_norm)
        
    # argument of periapse
    if e_vec[2] >= 0:
        omega = np.arccos(np.inner(n_vec,e_vec)/(n_norm*e_norm))
    else:
        omega = 2*np.pi - np.arccos(np.inner(n_vec,e_vec)/(n_norm*e_norm))
        
    M = E - e*np.sin(E)                      # mean anomaly
    a = 1 / (2/r_norm - (v_norm**2)/mu)      # semi major axis
    p = a*(1-e**2)                           # semi parameter (works for parabolic orbits)
    
    stateKep = [p,e,np.rad2deg(i),np.rad2deg(Omega),np.rad2deg(omega),np.rad2deg(v)]
    
    '''
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/oscltx_c.html
    Input
    
    
    Output
    RP      Perifocal distance. a(1-e)
    ECC     Eccentricity.
    INC     Inclination.
    LNODE   Longitude of the ascending node.
    ARGP    Argument of periapsis.
    M0      Mean anomaly at epoch.
    T0      Epoch.
    MU      Gravitational parameter.
    NU      True anomaly at epoch.
    A       Semi-major axis. A is set to zero if
            it is not computable.
    TAU     Orbital period. Applicable only for
            elliptical orbits. Set to zero otherwise.
    '''
    # t = self.epochJ2000+self.step*self.dt
    # print(self.step*self.dt)
    # rp,e,i,LAN,w,M,epoch,mu,v,a,n = spice.oscltx(stateECI,t,self.cb['mu'])
    # a = rp / (1-e)
    # stateKep = [a,e,np.rad2deg(i),np.rad2deg(w),np.rad2deg(LAN),np.rad2deg(M)]
    
    return stateKep

def ECItoECEF(stateECI,epoch,time,frame='J2000'):
    # convert from earth inertial to earth fixed frame
    # source: https://www.youtube.com/watch?v=wqhKcV-BubE&ab_channel=AlfonsoGonzalez-Astrodynamics%26SEPodcast
    
    # get times in j2000 frame
    epoch_j2000 = (epoch-j2000).total_seconds() # time of beginning epoch in j2000
    currTimeJ2000 = time + epoch_j2000 # current time in j2000 reference
    
    # split into position and velocity states
    ECIpos = stateECI[0:3];
    ECIvel = stateECI[3:6];
    
    # calculate ECEF
    Rot_matx = spice.pxform(frame,'ITRF93',currTimeJ2000)
    ECEFpos = np.dot(Rot_matx,ECIpos)
    ECEFvel = np.dot(Rot_matx,ECIvel)
    
    # return full state
    stateECEF = np.concatenate([ECEFpos, ECEFvel], -1)
    
    return stateECEF
    
def ECEFtoLatLonAlt(stateECEF,r):
    '''
    calculates all lat lon states for given ECEF state
    R = body radius
    '''
    
    r_norm, lat, lon = spice.reclat(stateECEF[0:3])
    stateLatLonAlt = [np.degrees(lat), np.degrees(lon), r_norm-r]
    
    return stateLatLonAlt