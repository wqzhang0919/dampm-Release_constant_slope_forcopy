#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 12:33:23 2021

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

name3 = "Complex_terrain_elastic_slope40_phi30"
name2 = "Complex_terrain_elastic_slope38_phi30"
name1 = "Complex_terrain_elastic_slope35_phi30" 

path1 = '../Simulation/Complex_terrain/' + name1
path2 = '../Simulation/Complex_terrain/' + name2
path3 = '../Simulation/Complex_terrain/' + name3


info = np.loadtxt(path1 + "/info.txt")
Np = np.loadtxt(path1 + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

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


#Chqrgement des donnees

frames = np.arange(0, end_frame+1) # last argument not included
time = np.linspace(0, final_time, end_frame+1)
PST_length = PST_speed * (time - loading_time)
mask_loading = time >= loading_time


# we look at a cut in the slab of the same coordinate x

cohesion1 = np.zeros([NloopY, end_frame+1])
Y1 = np.zeros([NloopY, end_frame+1])
cohesion2 = np.zeros([NloopY, end_frame+1])
Y2 = np.zeros([NloopY, end_frame+1])
cohesion3 = np.zeros([NloopY, end_frame+1])
Y3 = np.zeros([NloopY, end_frame+1])


Initial_position_y = np.loadtxt(path1 + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 0]



for f in frames:
    particles_data1 = np.loadtxt(path1 + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    particles_data2 = np.loadtxt(path2 + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    particles_data3 = np.loadtxt(path3 + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    
    cohesion1[:, f] = particles_data1[:, 5]
    cohesion2[:, f] = particles_data2[:, 5]
    cohesion3[:, f] = particles_data3[:, 5]
    
    Y1[:, f] = particles_data1[:, 0]
    Y2[:, f] = particles_data2[:, 0]
    Y3[:, f] = particles_data3[:, 0]


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


theta_init = 35*pi/180
bascule = 30 + 10*(-tan(phi)/tan(theta_init) + 1)
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
    plt.plot(Initial_position_y, cohesion1[:, f])
    plt.show()


M = len(Initial_position_y)//2

def find_front(i):
    if i==1:
        cohesion = cohesion1
    elif i==2:
        cohesion = cohesion2
    elif i==3:
        cohesion = cohesion3
    front = np.zeros(end_frame+1)
    Y0 = Initial_position_y[:]
    for f in frames:
        # diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        diff_cohesion = abs( cohesion[M:-1, f] - cohesion[M+1:, f])
        n = np.argmax(diff_cohesion)
        front[f] = Y0[n]
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
    
    return new_front, speed_front

front1 = find_front(1)
front2 = find_front(2)
front3 = find_front(3)

new_front1, speed_front1 = speed_front(front1)
new_front2, speed_front2 = speed_front(front2)
new_front3, speed_front3 = speed_front(front3)

plt.figure(figsize=figsize)
plt.plot(new_front1 + 20, speed_front1, "b", label="theta=35")
plt.plot(new_front2 + 20, speed_front2, "r", label="theta=38")
plt.plot(new_front3 + 20, speed_front3, "g", label="theta=40")
plt.plot(new_front3 + 20, (Wave_speed/C_s)*np.ones(len(new_front3)), color="black", label="C_s")

plt.grid()
plt.legend()
plt.xlabel("position X")
plt.ylabel("Crack speed")
plt.title("Different crack speed profile as a function of the slope angle")
plt.vlines(x=bascule, ymin=0, ymax=2, color="grey", line="--", label="slope angle = phi")

plt.show()
