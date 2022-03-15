# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 22:32:39 2022

@author: AlecY
"""


from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax = fig.add_subplot(projection='3d')

# create the parametric curve
t=np.arange(0, 2*np.pi, 2*np.pi/100)
x=np.cos(t)
y=np.sin(t)
z=t/(2.*np.pi)

# create the first plot
point, = ax.plot([x[0]], [y[0]], [z[0]], 'o')
line, = ax.plot(x, y, z, label='parametric curve')
ax.legend()
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])

# second option - move the point position at every frame
def update_point(n, x, y, z, point):
    point.set_data(np.array([x[n], y[n]]))
    point.set_3d_properties(z[n], 'z')
    
    print('updating point:',n)
    
    return point

ani=animation.FuncAnimation(fig, update_point, 99, fargs=(x, y, z, point))

plt.show()