#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:39:30 2022

@author: de
"""


import electronvolt_num as units
import io, libconf
import numpy as np
import pyGITR
import netCDF4
from netCDF4 import Dataset
import os
import math


FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit = np.array[PositionData['surfaceHit']]


surface_vx = np.array(PositionData['vx'])
surface_vy = np.array(PositionData['vy'])
surface_vz = np.array(PositionData['vz'])


Energy_particles = np.array(0.5*184*1.66e-27*(surface_vx**2 + surface_vy**2 + surface_vz**2)/1.602e-19)

#%%
# Setting up of particles by mesh element.

#%%
filename="gitrGeom.cfg"

with io.open(filename) as f:
    config = libconf.load(f)

x1 = np.array(config.geom.x1)
x2 = np.array(config.geom.x2)
x3 = np.array(config.geom.x3)
y1 = np.array(config.geom.y1)
y2 = np.array(config.geom.y2)
y3 = np.array(config.geom.y3)
z1 = np.array(config.geom.z1)
z2 = np.array(config.geom.z2)
z3 = np.array(config.geom.z3)
area = np.array(config.geom.area)
surf = np.array(config.geom.surface)
Z = np.array(config.geom.Z)
materialSurfaceInidces = np.nonzero(Z)
surfIndArray = np.asarray(materialSurfaceInidces)
a = np.array(config.geom.a)
b = np.array(config.geom.b)
c = np.array(config.geom.c)
d = np.array(config.geom.d)
plane_norm = np.array(config.geom.plane_norm)
#%%

#Reading the surface features

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")


C_C = SurfaceConcentrationData['surface_concentration'][:,:]
Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)


# for looping over all the surface elements and tracking concentrations and carbon erosion
new_entry_C = np.zeros(len(x1))
Gamma_C_ero_global = np.zeros(len(x1))

#Calculate Flux_C and obtain Flux_H

Flux_H = 1.43222e+19 #Obtained from SOLPS simulation (Zack)
beta_C = 0  # zero redeposition
Delta_t = 0.1 # in seconds
Delta_t_gitr = 1e-7
Delta_implant = 0.2 # enter parameter value and units

#Calculate Edist_H

# The following arrays will keep track of entries to be made in the next GITR run. 
vx_array = np.zeros(0)
vy_array = np.zeros(0)
vz_array = np.zeros(0)

x_array = np.zeros(0)
y_array = np.zeros(0)
z_array = np.zeros(0)


amu = 12 #for carbon
n_atom = 1e24 # average number density

# to loop over only the non-zero entries for Edist
specified = np.nonzero(Edist)
nP_global = 0 #tracks total number of eroded particles

for (x,y,z) in zip(*specified):
    
    Energy_C = Edist[x,y,z]*gridE[y] # check this that the Energy is in eV
    Angle_C = gridA[z] # check this
    
    sr_object = Sputtering_and_reflection()
    
    Flux_C = Edist[x,y,z]/(Delta_t*area[x])  # check this
    #Gamma_W_ero = sr_object.Calculate_PhysicalSputteringParameters('H','W',Energy_H)*C_W[time]*Gamma_H_incident + sr_object.Calculate_PhysicalSputteringParameters('C','W',Energy_C)*C_W[time]*Gamma_C_incident 
    Gamma_C_ero = sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_C)*C_C[x,counter-1]*Flux_C
    
    Gamma_C_dep = (1- sr_object.Calculate_ReflectionCoefficients('C','C',Energy_C))*C_C[x,counter-1]*Flux_C
    Gamma_C_redep = beta_C*Gamma_C_ero
    Gamma_C_ero = (1-beta_C)*Gamma_C_ero    
    
    Gamma_C_bulk = 0
    Gamma_W_bulk = 0
    
    #Updating C_C
    
    Gamma_C_net = Gamma_C_dep - Gamma_C_ero
    Gamma_W_net = 0 # for the present case
    
    if (Gamma_C_net+Gamma_W_net) > 0: #deposition regime
        Gamma_C_bulk = C_C[x,counter-1]*(Gamma_C_net + Gamma_W_net)
        #Gamma_W_bulk = C_W[time]*(Gamma_C_net + Gamma_W_net)
    elif (Gamma_C_net+Gamma_W_net) < 0:#erosion regime
        Gamma_C_bulk = 0
        #Gamma_W_bulk = (Gamma_C_net + Gamma_W_net)
        
        
    new_entry_C[x] = new_entry_C[x] + (Gamma_C_net - Gamma_C_bulk)    
    Gamma_C_ero_global[x] = Gamma_C_ero_global[x] + Gamma_C_ero 
    
    no_of_C = math.ceil(Gamma_C_ero*Delta_t*area[x])
    
    #Let us say it is one for each 
    nP_global = nP_global + 1
    
    p_C = ParticleDistribution()
    p_C.SetAttr('Np', 1)
    
    p_C.SetAttr('x','Uniform',xmin=min(x1[x],x2[x]),xmax=max(x1[x],x2[x])) 
    p_C.SetAttr('y','Uniform',xmin=min(y1[x],y2[x]),xmax=max(y1[x],y2[x]))
    p_C.SetAttr('z','Uniform',xmin=min(z1[x],z2[x]),xmax=max(z1[x],z2[x]))
    
    p_C.SetAttr('vx',0.0)
    p_C.SetAttr('vy',0.0)
    
    p_C.SetAttr('vz','Thomson')
    E_C = p_C.Particles['vz']
    
    E_C = np.sqrt(2*E_C*units.eV/(amu*units.mp))
    p_C.SetAttr('vz',E_C)
        
    vx_array = np.append(vx_array,p_C.Particles['vx'])
    vy_array = np.append(vy_array,p_C.Particles['vy'])
    vz_array = np.append(vz_array,p_C.Particles['vz'])
    
    x_array = np.append(x_array,p_C.Particles['x'])
    y_array = np.append(y_array,p_C.Particles['y'])
    z_array = np.append(z_array,p_C.Particles['z'])
    
    

#C_C[x][counter] = C_C[x][counter-1] + Delta_t*new_entry_C[x]/(n_atom*Delta_implant)
print(C_C.shape)
C_C_new = (C_C[:,1]).reshape(len(x1),1)
C_C = np.concatenate((C_C,C_C_new),axis=1)
print(C_C.shape)

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t)

#nP_global = np.multiply(Gamma_C_ero,area)
#nP_global = nP_global*Delta_t
#N_P = np.sum(nP_global) # total number of carbon atoms genererated during erosion

p_C = ParticleDistribution()

p_C.SetAttr('Np', nP_global)

# Set positions of particles
p_C.SetAttr('x',x_array)
p_C.SetAttr('y',y_array)
p_C.SetAttr('z',z_array)

# Set velocities of particles
p_C.SetAttr('vx',vx_array)
p_C.SetAttr('vy',vy_array)
p_C.SetAttr('vz',vz_array)

# Writing the particle source file for the next GITR run.
p_C.WriteParticleFile('particleConf.nc')

#Writing the surface features with time

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_new.nc', 'w', format='NETCDF4')
s_number_dim = ncFile.createDimension('surface_dim', len(x1)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))
s_concentration = ncFile.createVariable('surface_concentration',np.float64,('surface_dim','time_dim'))


s_number[:] = np.linspace(1,len(x1),len(x1))
s_time[:] = Surface_time
s_concentration[:,:] = C_C

ncFile.close()



#run the next simulation
    
    