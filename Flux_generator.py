#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:50:33 2022

@author: audide
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generation of particles distribution for GITR.
@author: guterl
"""
import pyGITR
from pyGITR.Particles import *
ParticleFile='particleConf2.nc'
p = ParticleDistribution()
p_C = ParticleDistribution()
p_H = ParticleDistribution()



# Attributes of particles x,y,z,v,x,vy,vz are stored in the dictionary p.Particles
# Display the list of attributes of particles.


# First, set numbers of particles. This is not affecting existing attributes
# but only the generation of those attributes
p.SetAttr('Np', 10)



# Set distribution for attributes of particles
# First, show list of available distribution pdf
#p.ShowAvailablePdfs()

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


p.ShowAttr()


alpha_c = 0.2
sigma = 1
Delta_t = 0.001
n_atom = 10**23
Delta_implant = 1


Time = np.linspace(0,10,10000)

Gamma_C_ero = np.zeros(Time.size)
Gamma_W_ero = np.zeros(Time.size)
time = 0;
C_C = np.zeros(Time.size)
C_W = np.zeros(Time.size)
C_W[0] = 1

Gamma_W_redep = 0
Gamma_C_redep = 0


for i in p.Particles['vx']:
    Gamma_C_incident = Gamma_C_redep + alpha_c*(1/(sigma*Delta_t))
    Gamma_H_incident = Gamma_H_redep
    
    Gamma_W_ero[time] = Sputtering*C_W[time]*Gamma_H_incident + Sputtering*C_W[time]*Gamma_C_incident 
    Gamma_C_ero[time] = Sputtering*C_C[time]*Gamma_H_incident + Sputtering*C_C[time]*Gamma_C_incident 
    
    Gamma_C_dep = (1-R)*C_W[time]*Gamma_C_incident + (1-R)*C_C[time]*Gamma_C_incident 
    Gamma_C_redep = beta_C*Gamma_C_ero[time]
    Gamma_C_ero[time] = (1-beta_C)*Gamma_C_ero[time]
    
    
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
    


    
    
    























