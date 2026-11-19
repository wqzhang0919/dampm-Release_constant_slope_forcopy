#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 16:04:21 2021

@author: slab_guest
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#############################################
#POSTPROCESSING FOR BENCHMARK n 2
#############################################

from math import pi, cos, sin, tan, sqrt

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
figsize = (18, 12)


# name = "test new version of code_Non_Parallel"
<<<<<<< HEAD
name = "Parallel_test_L20_h1_T6_fifth"
name1 = "Non-Parallel_test_L20_h1_T6_second"

# path = "build/dumps/" + name
path = "build/dumps/" + name
path1 = "build/dumps/" + name1

# name = "parallel_24_prime"
# name1 = "Non_parallel_prime"

# # path = "build/dumps/" + name
# path = "/media/berclaz/HDD4/larisempm_louis/Test_parallel/" + name
# path1 = "/media/berclaz/HDD4/larisempm_louis/Test_parallel/" + name1


info = np.loadtxt(path1 + "/info.txt")
Np = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame

# Acharger dan le truc d'info

loading_time = 1
PST_speed = 1

# loading_time = 0.25
# PST_speed = 4



dx = float(info[2])
rho   = float(info[3])
g = float(info[4])
phi = float(info[5])
Cohesion_WL = float(info[6])
h = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[0, 6]
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


# We keep the values of the averaged cohesion for all particles and the angle theta of the slope 

# at each frame
mean_cohesion = np.zeros(end_frame+1)
theta = np.zeros(end_frame+1)

# we look at a cut in the slab of the same coordinate x
particle_middle = [4*NloopY*round(NloopX/2) + 2*j for j in range(2*NloopY)]
point_middle = 4*NloopY*round(NloopX/2) + 2*round(NloopY/2)+1
point_front = 4*NloopY*round(NloopX/2) + 4*NloopY

displacement = np.zeros([Np, end_frame+1])
displacement1 = np.zeros([Np, end_frame+1])

cohesion = np.zeros([2*NloopY, end_frame+1])
cohesion1 = np.zeros([2*NloopY, end_frame+1])

X = np.zeros([Np, end_frame+1])
Y = np.zeros([Np, end_frame+1])
Y1 = np.zeros([Np, end_frame+1])
tau_yz = np.zeros([Np, end_frame+1])
tau_yz1 = np.zeros([Np, end_frame+1])
mean_tau_yz = np.zeros(end_frame+1)
tau_yz_front = np.zeros(end_frame+1)
tau_xz = np.zeros([Np, end_frame+1])

area = np.zeros([Np, end_frame+1])

displacement_middle = np.zeros(end_frame+1)
mean_displacement = np.zeros(end_frame+1)
displacement_front = np.zeros(end_frame+1)
sigma_xx = np.zeros([Np, end_frame+1])
sigma_yy = np.zeros([Np, end_frame+1])
sigma_xy = np.zeros([Np, end_frame+1])

sigma_xx1 = np.zeros([Np, end_frame+1])
sigma_yy1 = np.zeros([Np, end_frame+1])
sigma_xy1 = np.zeros([Np, end_frame+1])

height = np.zeros([Np, end_frame+1])

velocity_x = np.zeros([Np, end_frame+1])
velocity_y = np.zeros([Np, end_frame+1])

velocity_x1 = np.zeros([Np, end_frame+1])
velocity_y1 = np.zeros([Np, end_frame+1])

Initial_position_x = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 0]
Initial_position_y = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 1]



for f in frames:
    particles_data = np.loadtxt(path + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    particles_data1 = np.loadtxt(path1 + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)

    # mean_cohesion[f] = np.mean( particles_data[:, 9] )
    # theta[f] = np.mean( particles_data[:, 8] )

    displacement[:, f] = - particles_data[:, 5]
    cohesion[:, f] = particles_data[particle_middle, 9]

    displacement1[:, f] = - particles_data1[:, 5]
    cohesion1[:, f] = particles_data1[particle_middle, 9]

    tau_yz[:, f] = particles_data[:, 10]
    tau_yz1[:, f] = particles_data1[:, 10]

    # tau_xz[:, f] = particles_data[:, 11]

    # mean_tau_yz[f] = np.mean( particles_data[:, 10] )
    # tau_yz_front[f] = particles_data[point_front, 10]

    # mean_displacement[f] = np.mean( -particles_data[point_middle, 5])

    # displacement_middle[f] = -particles_data[point_middle, 5]
    # displacement_front[f] = -particles_data[point_front, 5]

    sigma_yy[:, f] = particles_data[:, 12]
    sigma_xx[:, f] = particles_data[:, 13]
    sigma_xy[:, f] = particles_data[:, 14]

    sigma_yy1[:, f] = particles_data1[:, 12]
    sigma_xx1[:, f] = particles_data1[:, 13]
    sigma_xy1[:, f] = particles_data1[:, 14]

     
    # area[:, f] = particles_data[:, 15]
    # height[:, f] = particles_data[:, 6]
    

    # area[:, f] = particles_data[:, 15]
    # height[:, f] = particles_data[:, 6]


    velocity_x[:, f] = particles_data[:, 2]
    velocity_y[:, f] = particles_data[:, 3]

    velocity_x1[:, f] = particles_data1[:, 2]
    velocity_y1[:, f] = particles_data1[:, 3]

    # #for the verfication
    # X[:, f] = particles_data[:, 0]
    Y[:, f] = particles_data[:, 1]
    Y1[:, f] = particles_data1[:, 1]

#%%

frame = 200


Err_Y = (1/Np*end_frame)*np.sum(np.abs(Y-Y1))
Err_tau_yz = (1/Np*end_frame)*np.sum(np.abs(tau_yz - tau_yz1))

print(Err_Y)
print(Err_tau_yz)


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], -tau_yz[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], -tau_yz1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.ylabel("contrainte tau_yz")
plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_yy[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], sigma_yy1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.ylabel("contrainte sigma_yy")
plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_xx[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], sigma_xx1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.ylabel("contrainte sigma_xx")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_xy[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], sigma_xy1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.show()
plt.ylabel("contrainte sigma_xy")

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], velocity_y[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], velocity_y1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.ylabel("velocity v_y")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], velocity_x[particle_middle, frame])
plt.plot(Y1[particle_middle, frame], velocity_x1[particle_middle, frame], "g")
plt.xlabel("position Y")
plt.ylabel("velocity v_x")
plt.show()


#%%

tau_b = rho * h * g * (1/final_time) * np.multiply(time, np.sin(theta))


theta_init = 35*pi/180

tau_p = rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL
tau_g = rho * h * g * sin(theta_init)
critical_displacement = (rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL)/K_WL

#for benchmark num 2
Lambda = sqrt( E * h/((1 - nu*nu) * K_WL) )
critical_length = Lambda * (tau_p/tau_g - 1) * (tan(theta_init)/(tan(theta_init) - tan(phi)))
G = E/(2*(1 + nu))
C_s = sqrt( G / rho )


def trace_cohesion(f):
    plt.figure(figsize=figsize)
    plt.plot(cohesion[:, f])
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

def trace_mean_tau():
    plt.figure(figsize=figsize)
    plt.plot(time, -mean_tau_yz, "r")
    plt.plot(time, tau_b, "b")
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

def tau_displacement_front():
    plt.figure(figsize=figsize)
    mask1 = displacement_front > critical_displacement
    mask2 = displacement_front <= critical_displacement
    plt.plot(displacement_front[mask2], -tau_yz_front[mask2], "b+")
    plt.plot(displacement_front[mask1], -tau_yz_front[mask1], "r+")
    plt.show()

def mean_tau_displacement():
    plt.figure(figsize=figsize)
    mask1 = mean_displacement > critical_displacement
    mask2 = mean_displacement <= critical_displacement
    plt.plot(mean_displacement[mask2], -mean_tau_yz[mask2], "b+")
    plt.plot(mean_displacement[mask1], -mean_tau_yz[mask1], "r+")
    plt.show()    


def plot_displacement(f):

    fig, ax = plt.subplots(figsize=figsize)

    L = PST_length[f]
    #cohesion_f = cohesion[:, f]
    mean_disp =  displacement[particle_middle, f]
    Y = Initial_position_y[particle_middle]
    mask1 = Y <= L
    mask2 = mean_disp <= critical_displacement
    cri = critical_displacement * np.ones(len(particle_middle))
    plt.plot(Y[mask2], mean_disp[mask2], "b+")
    plt.plot(Y[~mask2], mean_disp[~mask2], "r+")
    plt.plot(Y[mask1], cri[mask1], "g")
    #plt.ylim([-1e-5, 5e-4])

    plt.show()   


#%%

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
    Y0 = Initial_position_y[particle_middle]
    for f in frames:
        diff_cohesion = abs( cohesion[:, f] - np.concatenate( (cohesion[1:, f], [0]) ))
        n = np.argmax(diff_cohesion)
        front[f] = Y0[n]
    return front

def speed_front():
    front = find_front()
    mask = (front > dx/4) & (front < L1-dx/4)
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
    plt.xlabel("position Y")
    plt.ylabel("vitesse de la rupture")
    plt.title("dx / lambda = " + str(dx/Lambda))
    plt.show()



#%%

frame = 200

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_yy[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("contrainte sigma_yy")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], -tau_yz[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("contrainte tau_yz")
plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], tau_xz[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("contrainte tau_xz")
plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], area[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("area")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], area[particle_middle, frame]*height[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("volume")
plt.show()


plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_xx[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("contrainte sigma_xx")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], sigma_xy[particle_middle, frame])
plt.xlabel("position Y")
plt.show()
plt.ylabel("contrainte sigma_xy")

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], velocity_y[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("velocity v_y")
plt.show()

plt.figure(figsize=figsize)
plt.plot(Y[particle_middle, frame], velocity_x[particle_middle, frame])
plt.xlabel("position Y")
plt.ylabel("velocity v_x")
plt.show()

tau_displacement(point_middle)
tau_displacement(point_front)
tau_displacement(7)

#%%
# fig = plt.figure(figsize=figsize)
# a1 = fig.add_axes([0,0,1,1])
# a1.set_ylim(-500, 300)
# for f in frames:
#     data = tau_yz[particle_middle, f]
#     if f == 0:
#         line, = plt.plot(Y[particle_middle, f], data)
#     else:
#         line.set_ydata(data)
#     plt.pause(0.01) # pause avec duree en secondes
#     plt.xlabel("position Y")
#     plt.ylabel("stress tau_yz at frame " + str(f))
# plt.show()

#%%

fig, ax = plt.subplots(figsize=figsize)
plt.plot(Y[particle_middle, frame], -tau_yz[particle_middle, frame])
plt.ylim([0, 1300])
#plt.xlim([2.5, 5])
plt.xlabel("position Y")
plt.ylabel("contrainte tau_yz")
plt.show()


