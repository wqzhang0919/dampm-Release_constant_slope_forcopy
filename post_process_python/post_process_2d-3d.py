#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 10:36:19 2021

@author: slab_guest
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#############################################
#POSTPROCESSING FOR DEPTH-AVERAGED AND 3D COMPARAISON
#############################################


# INFOS, PARAMETERS AND USEFUL FUNCTIONS

from math import pi, cos, sin, tan, sqrt

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
figsize = (18, 12)


# name = "test_complex_terrain_prime"
# name = "p0_5000_theta_45"
# name = "p0_8000_theta_40"
name = "PST_MCC_p0_15000_M12_dx1_light"
# name = "Test_RMA_p0_8000_xi_1_compressible_prime"

#path = "build/dumps/" + name
# path = "../Simulation/Benchmark3/" + name
path = '../Simulation/Comparaison_Bertil_2D/' + name

info = np.loadtxt(path + "/info.txt")
Np = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

# end_frame = int(info[0])
end_frame= 293
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


displacement = np.zeros([NloopY, end_frame+1])
p0 = np.zeros([NloopY, end_frame+1])
cohesion = np.zeros([NloopY, end_frame+1])
Y = np.zeros([NloopY, end_frame+1])
tau_yz = np.zeros([NloopY, end_frame+1])
tau_xz = np.zeros([NloopY, end_frame+1])

P = np.zeros([NloopY, end_frame+1])
P_init = np.zeros([NloopY, end_frame+1])
Q = np.zeros([NloopY, end_frame+1])
rho_p = np.zeros([NloopY, end_frame+1])
area = np.zeros([NloopY, end_frame+1])
mass = np.zeros([NloopY, end_frame+1])

sigma_xx = np.zeros([NloopY, end_frame+1])
sigma_yy = np.zeros([NloopY, end_frame+1])
sigma_xy = np.zeros([NloopY, end_frame+1])
height = np.zeros([NloopY, end_frame+1])

velocity_x = np.zeros([NloopY, end_frame+1])
velocity_y = np.zeros([NloopY, end_frame+1])
Eps = np.zeros([NloopY, end_frame+1])

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
    p0[:, f] = particles_data[:, 12]
    velocity_x[:, f] = particles_data[:, 1]
    velocity_y[:, f] = particles_data[:, 2]

    Y[:, f] = particles_data[:, 0]
    P[:, f] = particles_data[:, 13]
    P_init[:, f] = particles_data[:, 14]
    Q[:, f] = particles_data[:, 15]
    
    Eps[:, f] = particles_data[:, 16]
    
    area[:, f] = particles_data[:, 17]
    rho_p[:, f] = particles_data[:, 18]
    mass[:, f] = area[:, f]*rho_p[:, f]*height[:, f]



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

def trace_displacement(p):
    plt.figure(figsize=figsize)
    plt.plot(time, displacement[p, :])
    plt.xlabel("temps")
    plt.ylabel("deplacement u_y")
    plt.show()
    plt.show()

def trace_tau(p):
    plt.figure(figsize=figsize)
    plt.plot(time, -tau_yz[p, :], "b")
    plt.xlabel("temps")
    plt.ylabel("tau_zy")
    #plt.plot(time, tau_b, "b")
    plt.show()

    
def tau_displacement(p):
    plt.figure(figsize=figsize)
    mask1 = displacement[p, :] > critical_displacement
    mask2 = displacement[p, :] <= critical_displacement
    plt.plot(displacement[p, mask2], -tau_yz[p, mask2], "b+")
    plt.plot(displacement[p, mask1], -tau_yz[p, mask1], "r+")
    plt.xlabel("displacement")
    plt.ylabel("tau_yz")
    plt.show()
    
    
    
def plot_displacement(f):
    
    fig, ax = plt.subplots(figsize=figsize)
    
    L = PST_length[f]
    #cohesion_f = cohesion[:, f]
    mean_disp =  displacement[:, f]
    Y = Initial_position_y[:]
    mask1 = Y <= L
    mask2 = mean_disp <= critical_displacement
    cri = critical_displacement * np.ones(len(Y))
    plt.plot(Y[mask2], mean_disp[mask2], "b+")
    plt.plot(Y[~mask2], mean_disp[~mask2], "r+")
    plt.plot(Y[mask1], cri[mask1], "g")
    #plt.ylim([-1e-5, 5e-4])

    plt.show()   

def find_speed_cohesion(f1, f2):
    t1 = f1 * frame_dt
    t2 = f2 * frame_dt
    diff_cohesion1 = abs( cohesion[:, f1] - np.concatenate( (cohesion[1:, f1], [0]) ))
    diff_cohesion2 = abs( cohesion[:, f2] - np.concatenate( (cohesion[1:, f2], [0]) ))
    n1 = np.argmax(diff_cohesion1)
    n2 = np.argmax(diff_cohesion2)
    deltaY = L1*(n2 - n1)/60
    deltaT = t2 - t1
    return deltaY/deltaT

def find_front():
    front = np.zeros(end_frame+1)
    Y0 = Initial_position_y[:]
    for f in frames:
        # diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        diff_cohesion = abs( cohesion[:-1, f] - cohesion[1:, f])
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
    
    mask1 = new_front < critical_length
    plt.figure(figsize=figsize)
    plt.plot(new_front[mask1], speed_front[mask1], "r")
    plt.plot(new_front[~mask1], speed_front[~mask1], "b")
    plt.plot(new_front, (Wave_speed/C_s)*np.ones(len(new_front)), "g")
    plt.xlabel("position Y")
    plt.ylabel("vitesse de la rupture")
    plt.title("dx / lambda = " + str(dx/Lambda))
    plt.show()

def find_critical_length(threshold):
    #on considere aue l'on a trouver le longueur critiue lorsdaue le front est superieur a un seuil + longueur PST
    front = find_front()
    diff = (front[mask_loading] - PST_length[mask_loading]) - threshold
    frame_simulation = frames[mask_loading]
    frame_aux = frame_simulation[diff > 0]
    crit_length = PST_length[frame_aux[0]]
    return crit_length
    

#%%

frame = 160


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


# plt.figure(figsize=figsize)
# plt.plot(Y[:, frame], height[:, frame])
# plt.xlabel("position Y")
# plt.ylabel("hauteur h")
# plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[:, frame], p0[:, frame])
plt.xlabel("position Y")
plt.ylabel("p0")
plt.show()

# plt.figure(figsize=figsize)
# plt.plot(Y[particle_middle, frame], area[particle_middle, frame]*height[particle_middle, frame])
# plt.xlabel("position Y")
# plt.ylabel("volume")
# plt.show()


# plt.figure(figsize=figsize)
# plt.plot(Y[particle_middle, frame], sigma_xx[particle_middle, frame])
# plt.xlabel("position Y")
# plt.ylabel("contrainte sigma_xx")
# plt.show()

# plt.figure(figsize=figsize)
# plt.plot(Y[particle_middle, frame], sigma_xy[particle_middle, frame])
# plt.xlabel("position Y")
# plt.show()
# plt.ylabel("contrainte sigma_xy")

# plt.figure(figsize=figsize)
# plt.plot(Y[particle_middle, frame], velocity_y[particle_middle, frame])
# plt.xlabel("position Y")
# plt.ylabel("velocity v_y")
# plt.show()

# plt.figure(figsize=figsize)
# plt.plot(Y[particle_middle, frame], velocity_x[particle_middle, frame])
# plt.xlabel("position Y")
# plt.ylabel("velocity v_x")
# plt.show()


#%%
fig = plt.figure(figsize=figsize)
a1 = fig.add_axes([0,0,1,1])
# a1.set_ylim(-critical_displacement, 20*critical_displacement)
a1.set_ylim(0, 3000)
plt.grid()
for f in frames:
    data = p0[:, f]
    if f == 0:
        line, = plt.plot(Initial_position_y, data)
    else:
        line.set_ydata(data)
    plt.pause(0.01) # pause avec duree en secondes
    plt.xlabel("position Y")
    plt.ylabel("diplacement at frame " + str(f))
plt.show()

#%%

f = 286


fig, (ax1, ax2) = plt.subplots(2, figsize=figsize)
plt.title("Crack arrest à la fin de la simulation. en haut : geometrie de la pente ; en bas : cohesion. Theta = 35")
ax1.set_xlim(0, 20)
ax2.set_xlim(0, 20)
plt.xlabel("y")
plt.setp(ax1, ylabel="z")
plt.setp(ax2, ylabel="hauteur h")
ax1.grid()
ax2.grid()
# ax3.grid()
# Y1 = Y[:, f]
ax1.plot(Y[:, f], p0[:, f], lw=1, color='b')
# ax1.plot(Y[:, f], tau_p*np.ones(len(Y[:, f])), 'r')
ax2.plot(Y[:, f], mass[:, f], lw=1, color='g')
#ax2.plot(Y[:, f], displacement[:, f], lw=1, color='g')
#ax2.plot(Y[:, f], critical_displacement*np.ones(len(Y[:, f])), 'r')
#ax3.plot(Y[:, f], velocity_y[:, f], lw=1, color='g')
plt.show()


    
