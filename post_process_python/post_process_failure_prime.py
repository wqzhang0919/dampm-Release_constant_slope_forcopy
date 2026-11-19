##################
#POSTPROCESSING FOR BENCHMARK n 2
#############################################

from math import pi, cos, sin, tan, sqrt

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
figsize = (18, 12)


# name = "test_complex_terrain_prime"
name = "Real_complex_terrain_fifth"
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

# We keep the values of the averaged cohesion for all particles and the angle theta of the slope 
# at each frame
mean_cohesion = np.zeros(end_frame+1)
theta = np.zeros(end_frame+1)

# we look at a cut in the slab of the same coordinate x
particle_middle = [4*NloopY*round(NloopX/2) + 2*j for j in range(2*NloopY)]
point_middle = 4*NloopY*round(NloopX/2) + 2*round(NloopY/2)+1
point_front = 4*NloopY*round(NloopX/2) + 4*NloopY
displacement = np.zeros([Np, end_frame+1])
cohesion = np.zeros([2*NloopY, end_frame+1])
X = np.zeros([Np, end_frame+1])
Y = np.zeros([Np, end_frame+1])
tau_yz = np.zeros([Np, end_frame+1])
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
height = np.zeros([Np, end_frame+1])

velocity_x = np.zeros([Np, end_frame+1])
velocity_y = np.zeros([Np, end_frame+1])

Initial_position_x = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 0]
Initial_position_y = np.loadtxt(path + "/out_part_frame_0.csv", delimiter=",", skiprows=1)[:, 1]



for f in frames:
    particles_data = np.loadtxt(path + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    
    mean_cohesion[f] = np.mean( particles_data[:, 9] )
    theta[f] = np.mean( particles_data[:, 8] )
    
    displacement[:, f] = - particles_data[:, 5]
    cohesion[:, f] = particles_data[particle_middle, 9]
    
    tau_yz[:, f] = particles_data[:, 10]
    tau_xz[:, f] = particles_data[:, 11]
    
    mean_tau_yz[f] = np.mean( particles_data[:, 10] )
    tau_yz_front[f] = particles_data[point_front, 10]

    mean_displacement[f] = np.mean( -particles_data[point_middle, 5]) 
    displacement_middle[f] = -particles_data[point_middle, 5]
    displacement_front[f] = -particles_data[point_front, 5]

    sigma_yy[:, f] = particles_data[:, 12]
    sigma_xx[:, f] = particles_data[:, 13]
    sigma_xy[:, f] = particles_data[:, 14]
    
    area[:, f] = particles_data[:, 15]
    height[:, f] = particles_data[:, 6]
    
    velocity_x[:, f] = particles_data[:, 2]
    velocity_y[:, f] = particles_data[:, 3]

    #for the verfication
    X[:, f] = particles_data[:, 0]
    Y[:, f] = particles_data[:, 1]
    
    
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
theta_init = 45*pi/180
tau_p = 1450

#critical_displacement = (rho*g*h*cos(theta_init)*tan(phi) + Cohesion_WL)/K_WL

#for benchmark num 2
Lambda = sqrt( E * h/((1 - nu*nu) * K_WL) )
#critical_length = Lambda * (tau_p/tau_g - 1) * (tan(theta_init)/(tan(theta_init) - tan(phi)))
G = E/(2*(1 + nu))
C_s = sqrt( G / rho )
Wave_speed = sqrt(E/rho)

def trace_cohesion(f):
    plt.figure(figsize=figsize)
    plt.plot(cohesion[:, f])
    plt.show()
    
#%%
f = 291

Y1 = Y[particle_middle, f]
S = sigma_yy[particle_middle, f]
T = -tau_yz[particle_middle, f]
H = height[particle_middle, f]
A = area[particle_middle, f]

fig, (ax1, ax2) = plt.subplots(2, figsize=figsize)
ax1.set_xlim(0, 40)
ax2.set_xlim(0, 40)
#ax1.set_ylim(0, 30)
#ax2.set_ylim(-2000, 1000)
ax1.grid()
ax2.grid()
# Y1 = Y[:, f]
# ax1.plot(Y[:, f], coord_z_np(Y[:, f], theta_init), lw=0.1, marker=3, color="b")
ax1.plot(Y1, S, lw=1, color='b')
ax2.plot(Y1, T, lw=1, color='g')
plt.show()