#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 15:31:25 2021

@author: slab_guest
"""

#########################################
#Comparasion critical_legnth pseudo 2d-3d
#########################################

import matplotlib.pyplot as plt
from math import cos, sin, sqrt, tan, pi
import numpy as np

figsize = (9, 6)

name = "tensile_length.csv"
path = "../Simulation/Benchmark3/"

info = np.loadtxt(path + "p0_8000_theta_30" + "/info.txt")

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

data = np.loadtxt(path + name, delimiter=",", skiprows=1)
Theta = data[:, 0]
Super_critical_length = data[:, 1]
critical_length = data[:, 2]

Lambda = sqrt( E * h/((1 - nu*nu) * K_WL) )

#critical_displacement = (rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL)/K_WL

#for benchmark num 2
G = E/(2*(1 + nu))
C_s = sqrt( G / rho )
Wave_speed = sqrt(E/rho)
sigma_t = 4363

def compute_tensile_length(theta):
    a_critic = sigma_t/(rho*g*(np.sin(theta) - np.cos(theta)*tan(phi)))
    return a_critic

angles = np.linspace(35, 60, 500)
theoric = compute_tensile_length(angles*pi/180)

plt.figure(figsize=figsize)
plt.plot(Theta, Super_critical_length, 'bo', label="quasi-static", markerfacecolor='none', linewidth=5)
plt.plot(Theta, critical_length, 'rs', label="dynamic", linewidth=5)
plt.plot(angles, theoric, 'k--', label="shear model")
plt.xlabel("slope angle $\Psi$ ($^\circ$)", fontsize=12)
plt.ylabel("crack length (m)", fontsize=12)
plt.legend(fontsize=12)
plt.show()