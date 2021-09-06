# Simulation of Heat Equation suing Finite Difference Method

import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

# Plate Dimension and time initialisation
size = 100
max_time = 500
dx = 1                  # Division difference
K = 40
dt = (dx**2)/(4*K)      # Stability Condition

Thot,Tcool = 100,0

# Initialisation of Array
U = np.empty((max_time,size,size))
U[:,:,:] = Tcool

# Initial Temperatures
T_right = 100
T_left = 100
T_top = 100
T_bottom = 100

# Initial Boundary Conditions
U[:,25:75,25:75] = Thot
"""
U[:,size-1,:] = T_bottom
U[:,:,size-1] = T_right
U[:,:,0] = T_left
U[:,0,:] = T_top
"""
# Calculating the numerical solution
for k in range(max_time-1):
    for i in range(size-1):
        for j in range(size-1):
            Laplace = U[k,i+1,j] + U[k,i-1,j] + U[k,i,j+1] + U[k,i,j-1] - 4*U[k,i,j]
            U[k+1,i,j] = U[k,i,j] + (K*dt)/(dx**2)*Laplace

fig,ax = plt.subplots()
title = plt.title("Temperature")
vector = ax.imshow(U[0],cmap=plt.cm.jet,vmax=Thot, vmin=0.0,interpolation="bilinear")
bar = plt.colorbar(vector)
bar.set_label("Temperature in degree Celsius")

ax.set_aspect("equal")

def animate(k):
    plt.title(f"Temperature at t = {k*dt:.3f} unit time")
    vector.set_array(U[k])
    return vector,
anim = FuncAnimation(fig,animate,interval=40,frames=max_time,repeat=True)
anim.save("Heat_Equation.mp4",writer="ffmpeg")
plt.show()
