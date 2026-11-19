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

figsize = (18, 12)

name = "critical_length.csv"
path = "../Simulation/comparaison_Bertil/Phi_21/"

info = np.loadtxt(path + "Comparaison_2d_3d_theta30" + "/info.txt")

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


def compute_critical_length(theta):
    tau_p = rho*g*h*np.cos(theta*pi/180)*tan(phi) + Cohesion_WL
    tau_g = rho * h * g * np.sin(theta*pi/180)
    a_critic = Lambda * (tau_p/tau_g - 1) * (np.tan(theta*pi/180)/(np.tan(theta*pi/180) - tan(phi)))
    return a_critic

# angles = np.linspace(24, 55, 500)*pi/180
angles = np.linspace(24, 55, 500)
theoric = compute_critical_length(angles)

plt.figure(figsize=figsize)
plt.plot(Theta[:-1]*180/pi, Super_critical_length[:-1], 'bo', label="super Critical length", markerfacecolor='none', linewidth=3)
plt.plot(Theta[:-1]*180/pi, critical_length[:-1], 'rs', label="critical length", linewidth=3)
plt.plot(angles, theoric, 'k--', label="theoric")
plt.xlabel("slope angle $\Psi$ ($^\circ$)", fontsize=12)
plt.ylabel("crack length (m)", fontsize=12)
plt.legend(fontsize=12)
plt.show()