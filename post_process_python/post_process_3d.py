#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:24:42 2021

@author: slab_guest
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
figsize = (18, 12)
from math import *


name = "elastic_slope45_light"
# name = "MCC_p0_10000_Light_dx_40"
# name = "elastic_Light_dx_40"

#path = "build/dumps/" + name
#path = "../Simulation/comparaison_Bertil/Phi_21/" + name
path = '../Simulation/Comparaison_Bertil_2D/' + name

info = np.loadtxt(path + "/info.txt")
Np = np.loadtxt(path + "/out_part_vertical_frame_0.csv", delimiter=",", skiprows=1).shape[0]

end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame

# Acharger dan le truc d'info
loading_time = 0.25
#PST_speed = 3.5
PST_speed = float(info[12])

theta_init = float(info[13])*pi/180
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


# displacement = np.zeros([NloopY, end_frame+1])
cohesion_v = np.zeros([2*NloopY, end_frame+1])
cohesion_h = np.zeros([2*NloopY, end_frame+1])
Y_v = np.zeros([2*NloopY, end_frame+1])
X_h = np.zeros([2*NloopY, end_frame+1])
Y_h = np.zeros([2*NloopY, end_frame+1])
X_v = np.zeros([2*NloopY, end_frame+1])
# tau_yz = np.zeros([NloopY, end_frame+1])
# tau_xz = np.zeros([NloopY, end_frame+1])

# sigma_xx = np.zeros([NloopY, end_frame+1])
# sigma_yy = np.zeros([NloopY, end_frame+1])
# sigma_xy = np.zeros([NloopY, end_frame+1])

# velocity_x = np.zeros([NloopY, end_frame+1])
# velocity_y = np.zeros([NloopY, end_frame+1])

Initial_position_y = np.loadtxt(path + "/out_part_vertical_frame_0.csv", delimiter=",", skiprows=1)[:, 1]
Initial_position_x = np.loadtxt(path + "/out_part_horizontal_frame_0.csv", delimiter=",", skiprows=1)[:, 0]



for f in frames:
    particles_data_v = np.loadtxt(path + "/out_part_vertical_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    particles_data_h = np.loadtxt(path + "/out_part_horizontal_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    
    cohesion_v[:, f] = particles_data_v[:, 6]
    cohesion_h[:, f] = particles_data_h[:, 6]

    Y_v[:, f] = particles_data_v[:, 1]
    X_v[:, f] = particles_data_v[:, 0]
    Y_h[:, f] = particles_data_h[:, 1]
    X_h[:, f] = particles_data_h[:, 0]

    # tau_yz[:, f] = particles_data[:,6]
    # tau_xz[:, f] = particles_data[:, 7]
    
    
    # sigma_yy[:, f] = particles_data[:, 8]
    # sigma_xx[:, f] = particles_data[:, 9]
    # sigma_xy[:, f] = particles_data[:, 10]
    
    
    # velocity_x[:, f] = particles_data[:, 1]
    # velocity_y[:, f] = particles_data[:, 2]



#%%
# theta_init = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[0, 8]
# theta_init = 30*pi/180
tau_p = rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL
tau_g = rho * h * g * sin(theta_init)
critical_displacement = (rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL)/K_WL

#for benchmark num 2
Lambda = sqrt( E * h/((1 - nu*nu) * K_WL) )
critical_length = Lambda * (tau_p/tau_g - 1) * (tan(theta_init)/(tan(theta_init) - tan(phi)))
G = E/(2*(1 + nu))
C_s = sqrt( G / rho )
Wave_speed = sqrt(E/rho)

def trace_cohesion(f):
    plt.figure(figsize=figsize)
    plt.plot(Initial_position_y, cohesion[:, f])
    plt.show()

My = len(Y_v[:, 0])//2
Mx = len(X_h[:, 0])//2

def find_front_vertical():
    front = np.zeros(end_frame+1)
    Y0 = Initial_position_y[:]
    for f in frames:
        # diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        diff_cohesion = abs( cohesion_v[My:-1, f] - cohesion_v[My+1:, f])
        n = np.argmax(diff_cohesion)
        front[f] = Y0[n]
    return front

def find_front_horizontal():
    front = np.zeros(end_frame+1)
    X0 = Initial_position_x[:]
    for f in frames:
        # diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        diff_cohesion = abs( cohesion_h[Mx:-1, f] - cohesion_h[Mx+1:, f])
        n = np.argmax(diff_cohesion)
        front[f] = X0[n]
    return front    

def speed_front(front):
    mask = (front > dx) & (front < L2-dx)
    new_front = front[mask]
    time_front = time[mask]
    speed_front = np.zeros(len(new_front))
    for k in range(1, len(new_front)):
        delta_y = new_front[k] - new_front[k-1]
        delta_t = time_front[k] - time_front[k-1]
        speed_front[k] = (delta_y / delta_t) * (1/C_s)
    return speed_front, new_front
    
velocity_y, front_y = speed_front(find_front_vertical())
velocity_x, front_x = speed_front(find_front_horizontal())



# fig, (ax1, ax2) = plt.subplots(2, figsize=figsize)
# plt.title("Crack speed, Theta = 45")
# ax1.set_xlim(0, 20)
# ax2.set_xlim(0, 20)
# plt.xlabel("position X")
# plt.setp(ax1, ylabel="vertical speed")
# plt.setp(ax2, ylabel="horizontal speed")
# ax1.grid()
# ax2.grid()
# ax1.plot(front_y, velocity_y, color="b")
# ax2.plot(front_x, velocity_x, color="b")
# ax1.plot(front_y, (Wave_speed/C_s)*np.ones(len(front_y)), "g")
# ax2.plot(front_x, np.ones(len(front_x)), "r")
# plt.show()


fig, ax1 = plt.subplots(1, figsize=figsize)
# plt.title("Crack speed, Theta = 45")
ax1.set_xlim(0, 20)
plt.xlabel("position X", fontsize=12)
plt.ylabel("crack speed propagation normalized", fontsize=12)
# ax1.grid()
ax1.plot(front_y, velocity_y, color="r", label="down-slope")
ax1.plot(front_x, velocity_x, color="g", label="cross-slope")
plt.legend(fontsize=12)
plt.show()


# plt.figure(figsize=figsize)
# plt.plot(new_front[mask1], speed_front[mask1], "r")
# plt.plot(new_front[~mask1], speed_front[~mask1], "b")
# plt.plot(new_front, (Wave_speed/C_s)*np.ones(len(new_front)), "g")
# plt.xlabel("position Y")
# plt.ylabel("vitesse de la rupture")
# plt.title("dx / lambda = " + str(dx/Lambda))
# plt.show()

    

#%%

# frame = 230


# plt.figure(figsize=figsize)
# plt.plot(Y[:, frame], sigma_yy[:, frame])
# plt.xlabel("position Y")
# plt.ylabel("contrainte sigma_yy")
# plt.show()

# plt.figure(figsize=figsize)
# plt.plot(Y[:, frame], -tau_yz[:, frame])
# plt.xlabel("position Y")
# plt.ylabel("contrainte tau_yz")
# plt.show()


    
