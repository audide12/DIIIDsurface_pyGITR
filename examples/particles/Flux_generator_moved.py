#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:50:33 2022

@author: audide
"""


import pyGITR
from pyGITR.Particles import *
ParticleFile='particleConf2.nc'
p = ParticleDistribution()
p_C = ParticleDistribution()
p_H = ParticleDistribution()



# Attributes of particles x,y,z,v,x,vy,vz are stored in the dictionary p.Particles
# Display the list of attributes of particles.

p.SetAttr('Np', 10)


# Set positions of particles
p.SetAttr(['z','y'],'Uniform',xmin=-0.05,xmax=0.05) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p.SetAttr('x',-0.01) # set all values of x to -0.01

# Set velocities of particles
p.SetAttr(['vx'],'Gaussian',beta=1.6)
p.SetAttr(['vy','vz'],'Gaussian')

# Rescale velocity by characteristic velocity
vpara = 1#1e4
vperp = 1#1e5
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)


#p.ShowAttr()



# Parameters for the model
alpha_c = 0.2
sigma = 1
Delta_t = 0.1
n_atom = 10**23
Delta_implant = 1
beta_C = 0.1 



Time = np.linspace(0,10,100)

Gamma_C_ero = np.zeros(Time.size)
Gamma_W_ero = np.zeros(Time.size)
time = 0;
C_C = np.zeros(Time.size)
C_W = np.zeros(Time.size)
C_W[0] = 1
C_C[0] = 0


p_C.SetAttr('Np', Time.size)
p_H.SetAttr('Np', Time.size)

# Set positions of particles
p_C.SetAttr(['z','y'],'Uniform',xmin=-0.05,xmax=0.05) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p_C.SetAttr('x',-0.01) # set all values of x to -0.01

p_H.SetAttr(['z','y'],'Uniform',xmin=-0.05,xmax=0.05) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p_H.SetAttr('x',-0.01) # set all values of x to -0.01



# Set velocities of particles
p_C.SetAttr(['vx'],'Gaussian',sigma = 1.825e4,beta=3.16e14)
p_C.SetAttr(['vy','vz'],'Gaussian',sigma = 1.825e4,beta=3.16e14)

p_H.SetAttr(['vx'],'Gaussian',sigma = 6.324e4, beta = 9.12e13)
p_H.SetAttr(['vy','vz'],'Gaussian',sigma = 6.324e4, beta = 9.12e13)

#p_H.SetAttr(['vx'],'Gaussian',sigma = 1, beta=0)
#p_H.SetAttr(['vy','vz'],'Gaussian',sigma = 1, beta=0)


Gamma_W_redep = 0
Gamma_C_redep = 0




for i in range(Time.size-1):
    
    sr_object = Sputtering_and_reflection()
    
    #Gamma_C_incident = Gamma_C_redep + alpha_c*(1/(sigma*Delta_t))
    Gamma_C_incident = alpha_c*(1/(sigma*Delta_t))
    Gamma_H_incident = (1/(sigma*Delta_t))
    
    Energy_H = 0.5 * (p_H.Particles['vx'][time]**2+p_H.Particles['vy'][time]**2+p_H.Particles['vz'][time]**2) * (1.0 * 1e-8)   # in eV
    Energy_C = 0.5 * (p_C.Particles['vx'][time]**2 + p_C.Particles['vy'][time]**2 + p_C.Particles['vz'][time]**2 ) * (12 * 1e-8)   # in eV
    
    
    Gamma_W_ero[time] = sr_object.Calculate_PhysicalSputteringParameters('H','W',Energy_H)*C_W[time]*Gamma_H_incident + sr_object.Calculate_PhysicalSputteringParameters('C','W',Energy_C)*C_W[time]*Gamma_C_incident 
    Gamma_C_ero[time] = sr_object.Calculate_PhysicalSputteringParameters('H','C',Energy_H)*C_C[time]*Gamma_H_incident + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_C)*C_C[time]*Gamma_C_incident 
    
    Gamma_C_dep = (1- sr_object.Calculate_ReflectionCoefficients('C','W',Energy_C))*C_W[time]*Gamma_C_incident + (1- sr_object.Calculate_ReflectionCoefficients('C','C',Energy_C))*C_C[time]*Gamma_C_incident 
    Gamma_C_redep = beta_C*Gamma_C_ero[time]
    Gamma_C_ero[time] = (1-beta_C)*Gamma_C_ero[time]
    
    
    
    Gamma_C_bulk = 0
    Gamma_W_bulk = 0
    
    #Updating C_C and C_W
    
    Gamma_C_net = Gamma_C_dep - Gamma_C_ero[time]
    Gamma_W_net = - Gamma_W_ero[time]
    
    if (Gamma_C_net+Gamma_W_net) > 0: #deposition regime
        Gamma_C_bulk = C_C[time]*(Gamma_C_net + Gamma_W_net)
        Gamma_W_bulk = C_W[time]*(Gamma_C_net + Gamma_W_net)
    elif (Gamma_C_net+Gamma_W_net) < 0:#erosion regime
        Gamma_C_bulk = 0
        Gamma_W_bulk = (Gamma_C_net + Gamma_W_net)
        
        
    C_C[time+1] = C_C[time] + Delta_t*(Gamma_C_net - Gamma_C_bulk)/(n_atom*Delta_implant)
    C_W[time+1] = C_W[time] + Delta_t*(Gamma_W_net - Gamma_W_bulk)/(n_atom*Delta_implant)
    
    
    time = time + 1
    


    
    
plt.plot(Time,C_W)    























