# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 22:09:40 2020

@author: AlecY
"""

import numpy as np
import PIL

from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import cm

plt.style.use('dark_background')

def run_visualizer(propagators,central_body):
    
    print("Running visualizer")
    plt.close('all')

    # plot trajectory of states
    viz_3D = plot_trajs(central_body,propagators,animated=True)
    
    # plot Lat Lon
    viz_ground = plot_LatLon(propagators,animated=False)
    
    # plot ECI states
    plot_ECI_states(propagators)
    
    # plot Kep states
    plot_Kep_states(propagators,units="days")
    
    # plot ECEF states
    plot_ECEF_states(propagators,units="days")
    
    return viz_3D,viz_ground

def plot_body(ax,body,bluemarble=False):
    
    if bluemarble:
        # NOTE: blue marble will be fixed in ECI, which is not correct
        
        # define sphere mesh & smoothness
        N=20 # 20=standard
        stride=5
        sampling=8 # larger=lower resolution image
        
        # load bluemarble with PIL
        bm = PIL.Image.open('./resources/blue_marble.jpg')
        # it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept 
        bm = np.array( bm.resize([int(d/sampling) for d in bm.size]) )/256.
        
        # coordinates of the image - don't know if this is entirely accurate, but probably close
        img_lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
        img_lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180
        
        x = (np.outer(np.cos(img_lons),np.cos(img_lats)).T)*body['radius']
        y = (np.outer(np.sin(img_lons),np.cos(img_lats)).T)*body['radius']
        z = (np.outer(np.ones(np.size(img_lons)),np.sin(img_lats)).T)*body['radius']
        ax.plot_surface(x, y, z, rstride=stride, cstride=stride, facecolors = bm)
        
    else:
        # define sphere mesh & smoothness
        N=20 # 20=standard
        stride=1
        
        # get base sphere coordinates (in ECI frame)
        u = np.linspace(0, 2 * np.pi, N)
        v = np.linspace(0, np.pi, N)
        x = np.outer(np.cos(u), np.sin(v))*body['radius']
        y = np.outer(np.sin(u), np.sin(v))*body['radius']
        z = np.outer(np.ones(np.size(u)), np.cos(v))*body['radius']
        my_colors = cm.jet(z/np.amax(z))
        ax.plot_surface(x, y, z, linewidth=0.0, cstride=stride, rstride=stride, zorder=0, facecolors=my_colors )
        
    # draw ECI frame
    len = body['radius']*2
    ax.quiver(0,0,0,0,0,len,label="Z (ECI)",color='r',zorder=10)
    ax.quiver(0,0,0,0,len,0,label="Y (ECI)",color='g',zorder=10)
    ax.quiver(0,0,0,len,0,0,label="X (ECI)",color='w',zorder=10)

def plot_trajs(central_body,propagators,animated=False):
    
    # set up figure
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    text_ax = fig.add_axes([0.0,0.95,0.1,0.05]) # time text axis
    text_ax.axis("off")
    
    plot_lim = 40000 # km
    ax.set_xlim([-plot_lim,plot_lim])
    ax.set_ylim([-plot_lim,plot_lim])
    ax.set_zlim([-plot_lim,plot_lim])
    
    # plot earth
    plot_body(ax,central_body,bluemarble=False)
    
    if not animated:
        trajectory_display_size = 1.0
        
        for propagator in propagators:
            x = propagator.statesECI[:,0]
            y = propagator.statesECI[:,1]
            z = propagator.statesECI[:,2]
            ax.plot3D(x,y,z,label=propagator.name,zorder=10,linewidth=trajectory_display_size,linestyle='-')
        
        anim = None
    else:
        time = text_ax.text(0.5,0.5, str(0), ha="left", va="top")
        
        points = [None] * len(propagators)
        lines  = [None] * len(propagators)
        
        # plot first point and full trajectory as a line
        for idx,propagator in enumerate(propagators):
            x = propagator.statesECI[:,0]
            y = propagator.statesECI[:,1]
            z = propagator.statesECI[:,2]
            
            lines[idx],  = ax.plot(x, y, z, label=propagator.name)
            points[idx], = ax.plot([x[0]], [y[0]], [z[0]], 'o')
        
        # Updating function, to be repeatedly called by the animation
        def update_point(n):
            for idx,propagator in enumerate(propagators):
                #print('updating point:',n)
                x = propagator.statesECI[:,0]
                y = propagator.statesECI[:,1]
                z = propagator.statesECI[:,2]
                
                points[idx].set_data(np.array([x[n], y[n]]))
                points[idx].set_3d_properties(z[n], 'z')
            
            time.set_text("Sim Time: {:.3f} days".format(propagators[0].times[n][0]/86400.0))
            
            return points,time
    
        # create animation
        interval_time = 10 # ms    
        anim = animation.FuncAnimation(fig, update_point, frames=len(propagators[0].statesECI),
                                       interval=interval_time)
        
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    
    return anim
        
def plot_ECI_states(propagators):
    
    fig, ax = plt.subplots(3)
    
    for propagator in propagators:
        ax[0].plot(propagator.times,propagator.statesECI[:,0],label=propagator.name)
        ax[1].plot(propagator.times,propagator.statesECI[:,1],label=propagator.name)
        ax[2].plot(propagator.times,propagator.statesECI[:,2],label=propagator.name)
    
    fig.suptitle('ECI States', fontsize=14)
    ax[0].set_ylabel('X',rotation=0.0)
    ax[1].set_ylabel('Y',rotation=0.0)
    ax[2].set_ylabel('Z',rotation=0.0)
    
    ax[2].set_xlabel('time (sec)',rotation=0.0,ha='center')
    
    # only take the labels from the last axis
    handles, labels = ax[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    
def plot_Kep_states(propagators,units="days"):
    
    if units=="days":
        scale_factor=24 * 60 * 60.0
    else:
        scale_factor=1
    
    fig, ax = plt.subplots(6)
    
    for propagator in propagators:
        ax[0].plot(propagator.times/scale_factor,propagator.statesKep[:,0],label=propagator.name)
        ax[1].plot(propagator.times/scale_factor,propagator.statesKep[:,1],label=propagator.name)
        ax[2].plot(propagator.times/scale_factor,propagator.statesKep[:,2],label=propagator.name)
        ax[3].plot(propagator.times/scale_factor,propagator.statesKep[:,3],label=propagator.name)
        ax[4].plot(propagator.times/scale_factor,propagator.statesKep[:,4],label=propagator.name)
        ax[5].plot(propagator.times/scale_factor,propagator.statesKep[:,5],label=propagator.name)
    
    fig.suptitle('Kepler States', fontsize=14)
    ax[0].set_ylabel('semi-parameter',          rotation=15.0,ha='right')
    ax[1].set_ylabel('eccentricity',            rotation=15.0,ha='right')
    ax[2].set_ylabel('inclination (deg)',       rotation=15.0,ha='right')
    ax[3].set_ylabel('Lon of Asc Node (deg)',   rotation=15.0,ha='right')
    ax[4].set_ylabel('arg of peri (deg)',       rotation=15.0,ha='right')
    ax[5].set_ylabel('mean anomaly (deg)',      rotation=15.0,ha='right')
    
    ax[5].set_xlabel('time (%s)'.format(units),rotation=0.0,ha='center')
    
    # formatting
    for i in range(5):
        ax[i].set_xlim(left=0)
        ax[i].ticklabel_format(useOffset=False)
    
    ax[3].set_ylim([0-10,360+10])
    ax[4].set_ylim([0-10,360+10])
    ax[5].set_ylim([0,360])
        
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax[2].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    # only take the labels from the last axis
    handles, labels = ax[5].get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    
def plot_ECEF_states(propagators,units="days"):
    
    if units=="days":
        scale_factor=24 * 60 * 60.0
    else:
        scale_factor=1
    
    fig, ax = plt.subplots(3)
    
    for propagator in propagators:
        ax[0].plot(propagator.times/scale_factor,propagator.statesECEF[:,0],label=propagator.name)
        ax[1].plot(propagator.times/scale_factor,propagator.statesECEF[:,1],label=propagator.name)
        ax[2].plot(propagator.times/scale_factor,propagator.statesECEF[:,2],label=propagator.name)
    
    fig.suptitle('ECEF States', fontsize=14)
    ax[0].set_ylabel('X',rotation=0.0)
    ax[1].set_ylabel('Y',rotation=0.0)
    ax[2].set_ylabel('Z',rotation=0.0)
    
    ax[2].set_xlabel('time (%s)'.format(units),rotation=0.0,ha='center')
    
    # only take the labels from the last axis
    handles, labels = ax[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    
def plot_LatLon(propagators,animated=False):
    
    fig,ax = plt.subplots(1)
    points = [None] * len(propagators)
    anim = None
    
    # plot ground track for each propagator
    for idx,propagator in enumerate(propagators):
            
        lons = propagator.statesLatLonAlt[:,0]
        lats = propagator.statesLatLonAlt[:,1]
        
        # separate list into multiple lists
        # this removes line artifacts when crossing 180 to -180 boundary
        split_idxs = []
        for i in range(0, len(lons)-1):
            if lons[i] > 0.0 and lons[i + 1] < 0.0:
                split_idxs.append(i)
        
        # plot each pass separately
        for idx,split_idx in enumerate(split_idxs):
            if idx == 0:
                lats_pass = lats[0:split_idx+1]
                lons_pass = lons[0:split_idx+1]
                p = ax.plot(lons_pass,lats_pass,label=propagator.name)
            else:
                lats_pass = lats[split_idxs[idx-1]+1:split_idx]
                lons_pass = lons[split_idxs[idx-1]+1:split_idx]
                aColor = p[0].get_color()
                ax.plot(lons_pass,lats_pass,color=aColor)
    
    # plot animation of ground position over time
    if animated:
        
        for idx,propagator in enumerate(propagators):
            points[idx], = ax.plot([lons[0]], [lats[0]], 'o')
        
        # Updating function, to be repeatedly called by the animation
        def update_point(n):
            for idx,propagator in enumerate(propagators):
                #print('updating point:',n)
                lons = propagator.statesLatLonAlt[:,0]
                lats = propagator.statesLatLonAlt[:,1]
                
                points[idx].set_data(np.array([lons[n], lats[n]]))
            
            return points
        
        # create animation
        interval_time = 10 # ms    
        anim = animation.FuncAnimation(fig, update_point, frames=len(propagators[0].statesLatLonAlt),
                                       interval=interval_time)
    
    # plot coastlines
    # coast_coords = np.genfromtxt('./resources/coastlines.csv',delimiter=',')
    # ax.plot(coast_coords[:,0],coast_coords[:,1],'mo',markersize=0.3)
    
    # plot earth background image
    img = plt.imread('./resources/blue_marble.jpg')
    ax.imshow(img, extent=[-180, 180, -90, 90])
    
    #fig.suptitle('Lat Lon', fontsize=14)
    
    ax.set_ylabel('Lat (deg)',rotation=0.0)
    ax.set_xlabel('Lon (deg)',rotation=0.0,ha='center')
    
    ax.set_ylim([-90,90])
    ax.set_xlim([-180,180])
    
    # only take the labels from the last axis
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    
    return anim