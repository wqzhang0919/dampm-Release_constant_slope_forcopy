#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 16:14:04 2021

@author: slab_guest
"""




#############################################
#POSTPROCESSING FOR DEPTH-AVERAGED AND FAILURE
#############################################


# INFOS, PARAMETERS AND USEFUL FUNCTIONS

from math import pi, cos, sin, tan, sqrt

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
figsize = (18, 12)


# name = "test_complex_terrain_prime"
name = "Complex_terrain_elastic_slope35_phi30" 
# name = "Comparaison_2d_3d_theta45"

#path = "build/dumps/" + name
#path = "../Simulation/comparaison_Bertil/Phi_21/" + name
path = '../Simulation/Complex_terrain/' + name

info = np.loadtxt(path + "/info.txt")
Np = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame

# Acharger dan le truc d'info
loading_time = 0.25
#PST_speed = 3.5
PST_speed = float(info[12])


dx = float(info[2])
rho   = float(info[3])
g = float(info[4])
phi = float(info[5])
Cohesion_WL = float(info[6])
#h = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[0, 6]
h = 0.5
L1 = float(info[7])
L2 = float(info[8])
E = float(info[9])
nu = float(info[10])
K_WL = float(info[11])

NloopX = round(L1/dx)
NloopY = round(L2/dx)


frames = np.arange(0, end_frame+1) # last argument not included
time = np.linspace(0, final_time, end_frame+1)
PST_length = PST_speed * (time - loading_time)
mask_loading = time >= loading_time

# We keep the values of the averaged cohesion for all particles and the angle theta of the slope 
# at each frame
mean_cohesion = np.zeros(end_frame+1)


# we look at a cut in the slab of the same coordinate x

height = np.zeros([NloopY, end_frame+1])
displacement = np.zeros([NloopY, end_frame+1])
cohesion = np.zeros([NloopY, end_frame+1])
Y = np.zeros([NloopY, end_frame+1])
tau_yz = np.zeros([NloopY, end_frame+1])
tau_xz = np.zeros([NloopY, end_frame+1])

sigma_xx = np.zeros([NloopY, end_frame+1])
sigma_yy = np.zeros([NloopY, end_frame+1])
sigma_xy = np.zeros([NloopY, end_frame+1])

velocity_x = np.zeros([NloopY, end_frame+1])
velocity_y = np.zeros([NloopY, end_frame+1])

Initial_position_y = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 0]



for f in frames:
    particles_data = np.loadtxt(path + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    
    displacement[:, f] = - particles_data[:, 4]
    cohesion[:, f] = particles_data[:, 5]
    
    tau_yz[:, f] = particles_data[:,6]
    tau_xz[:, f] = particles_data[:, 7]
    
    
    sigma_yy[:, f] = particles_data[:, 8]
    sigma_xx[:, f] = particles_data[:, 9]
    sigma_xy[:, f] = particles_data[:, 10]
    
    height[:, f] = particles_data[:, 11]
    velocity_x[:, f] = particles_data[:, 1]
    velocity_y[:, f] = particles_data[:, 2]

    Y[:, f] = particles_data[:, 0]

# def coord_z(y):
#     #for y >= 0
#     return (80/pi) * (np.sin(pi/40*(y - 20)) + 1)

def coord_z(y, theta):
    if y < 0:
        return 0
    elif y <= 10:
        return tan(theta)/20 * y**2
    elif y<= 30:
        return tan(theta)*(y - 5)
    else:
        return tan(theta)*(-(y-30)**2/20 + (y-30) + 25)
    
def coord_z_np(Y, theta):
    Z = np.zeros(len(Y))
    for k in range(len(Y)):
        Z[k] = coord_z(Y[k], theta)
    return Z
#%%
# theta_init = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[0, 8]
# theta_init = 30*pi/
theta_init = 35*pi/180
tau_p = rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL
tau_g = rho * h * g * sin(theta_init)
critical_displacement = (rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL)/K_WL

#for benchmark num 2
Lambda = sqrt( E * h/((1 - nu*nu) * K_WL) )
#critical_length = Lambda * (tau_p/tau_g - 1) * (tan(theta_init)/(tan(theta_init) - tan(phi)))
G = E/(2*(1 + nu))
C_s = sqrt( G / rho )
Wave_speed = sqrt(E/rho)

def trace_cohesion(f):
    plt.figure(figsize=figsize)
    plt.plot(Y[:, f], cohesion[:, f])
    plt.show()

#%%

frame = 480


plt.figure(figsize=figsize)
plt.plot(Y[:, frame], coord_z_np(Y[:, frame], theta_init), lw=0.1, marker=3, color="b")
plt.xlabel("position Y")
plt.ylabel("position Z")
plt.show()



#%%
f = 480

fig, (ax1, ax2) = plt.subplots(2, figsize=figsize)
plt.title("Crack arrest à la fin de la simulation. en haut : geometrie de la pente ; en bas : cohesion. Theta = 35")
ax1.set_xlim(0, 40)
ax2.set_xlim(0, 40)
plt.xlabel("y")
plt.setp(ax1, ylabel="z")
plt.setp(ax2, ylabel="hauteur h")
ax1.grid()
ax2.grid()
# ax3.grid()
# Y1 = Y[:, f]
ax1.plot(Y[:, f], coord_z_np(Y[:, f], theta_init), lw=0.1, marker=3, color="b")
# ax1.plot(Y[:, f], -tau_yz[:, f], lw=1, color='b')
# ax1.plot(Y[:, f], tau_p*np.ones(len(Y[:, f])), 'r')
ax2.plot(Y[:, f], cohesion[:, f], lw=1, color='g')
#ax2.plot(Y[:, f], displacement[:, f], lw=1, color='g')
#ax2.plot(Y[:, f], critical_displacement*np.ones(len(Y[:, f])), 'r')
#ax3.plot(Y[:, f], velocity_y[:, f], lw=1, color='g')
plt.show()


#%%

M = len(Y[:, 0])//2

def find_front():
    front = np.zeros(end_frame+1)
    Y0 = Initial_position_y[:]
    for f in frames:
        # diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        diff_cohesion = abs( cohesion[M:-1, f] - cohesion[M+1:, f])
        n = np.argmax(diff_cohesion)
        front[f] = Y0[n]
    return front

def speed_front():
    front = find_front()
    mask = (front > dx) & (front < L2-dx)
    new_front = front[mask]
    time_front = time[mask]
    speed_front = np.zeros(len(new_front))
    for k in range(1, len(new_front)):
        delta_y = new_front[k] - new_front[k-1]
        delta_t = time_front[k] - time_front[k-1]
        speed_front[k] = (delta_y / delta_t) * (1/C_s)
        
    plt.figure(figsize=figsize)
    plt.plot(new_front, speed_front, "r")
    plt.plot(new_front, speed_front, "b")
    plt.plot(new_front, (Wave_speed/C_s)*np.ones(len(new_front)), "g")
    plt.xlabel("position Y")
    plt.ylabel("vitesse de la rupture")
    plt.title("dx / lambda = " + str(dx/Lambda))
    plt.show()



 #%%

import matplotlib.animation as animation

def data_gen():
    f = data_gen.f
    cnt = 0
    while cnt < 720:
        cnt+=1
        f += 1
        y1 = Y[:, f]
        y2 = sigma_yy[:, f]
        z = coord_z(Y[:, f])
        # adapted the data generator to yield both sin and cos
        yield z, y1, y2

data_gen.f = 0

# create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2,1)

# intialize two line objects (one in each axes)
line1, = ax1.plot([], [], lw=2, marker="o", color="b")
line2, = ax2.plot([], [], lw=2, color='r')
line = [line1, line2]

# the same axes initalizations as before (just now we do it for both of them)
ax1.set_xlim(0, 40)
ax2.set_xlim(0, 40)

ax1.set_ylim(0, 15)
ax2.set_ylim(-200, 500)

ax1.grid()
ax2.grid()
 
# for ax in [ax1, ax2]:
#     ax.set_ylim(-1.1, 1.1)
#     ax.set_xlim(0, 5)
#     ax.grid()

# initialize the data arrays 
xdata, y1data, y2data = [], [], []
def run(data):
    # update the data
    t, y1, y2 = data
    xdata.append(t)
    y1data.append(y1)
    y2data.append(y2)

    # axis limits checking. Same as before, just for both axes
    for ax in [ax1, ax2]:
        xmin, xmax = ax.get_xlim()
        if t >= xmax:
            ax.set_xlim(xmin, 2*xmax)
            ax.figure.canvas.draw()

    # update the data of both line objects
    line[0].set_data(xdata, y1data)
    line[1].set_data(xdata, y2data)

    return line

ani = animation.FuncAnimation(fig, run, data_gen, blit=True, interval=10,
    repeat=False)
plt.show()   
