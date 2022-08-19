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

#Calculate Flux_C and obtain Flux_H

Flux_H = 1.43222e+19 # Obtained from SOLPS simulation (Zack)
beta_C = 0  # zero redeposition
Delta_t = 0.1 # in seconds
Delta_t_gitr = 1e-7
Delta_implant = 1e-3 # enter parameter value and units
amu = 12 #for carbon
n_atom = 1e24 # average number density
weight_gitr = Delta_t/Delta_t_gitr
#weight_gitr = 1

#Calculate Edist_H

#%%
# Reading position files

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit = np.array(PositionData['surfaceHit'])


surface_vx = np.array(PositionData['vx'])
surface_vy = np.array(PositionData['vy'])
surface_vz = np.array(PositionData['vz'])


Energy_particles = np.array(0.5*amu*1.66e-27*(surface_vx**2 + surface_vy**2 + surface_vz**2)/1.602e-19) # make sure that this energy is in eV


count = 0
for i in surfacehit:
    if i != -1:
        count+=1
print(count)

#print(surface_vx)

#Angle_particles = add this functionality


#%%
#Reading geometry files

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

# This section is to initialize the C_C file

C_C = np.full((len(x1),1),0.4)

Surface_time = np.full((1,1),0.0)

Surface_number = np.array(range(len(x1)))

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C.nc', 'w', format='NETCDF4')
s_number_dim = ncFile.createDimension('surface_dim', len(x1)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))
s_concentration = ncFile.createVariable('surface_concentration',np.float64,('surface_dim','time_dim'))


s_number[:] = np.linspace(1,len(x1),len(x1))
s_time[:] = Surface_time
s_concentration[:,:] = C_C

ncFile.close()


#%%

#Reading the surface features

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")


C_C = SurfaceConcentrationData['surface_concentration'][:,:]
Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)


# for looping over all the surface elements and tracking concentrations and carbon erosion

Gamma_C_ero_global = np.zeros(len(x1))

Gamma_C_dep_global = np.zeros(len(x1))



#%%
# Calculation of erosion fluxes

for i in range(len(Energy_particles)):
    if surfacehit[i] != -1:
        surface_index = int(surfacehit[i])
        sr_object = Sputtering_and_reflection()
        
        Flux_C = weight_gitr/(Delta_t_gitr*area[surface_index])  # check this
        #Gamma_W_ero = sr_object.Calculate_PhysicalSputteringParameters('H','W',Energy_H)*C_W[time]*Gamma_H_incident + sr_object.Calculate_PhysicalSputteringParameters('C','W',Energy_C)*C_W[time]*Gamma_C_incident 
        Gamma_C_ero = sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles[i])*C_C[surface_index,-1]*Flux_C
        #print(Flux_C)
        #print(sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles[i]))
        
        
        Gamma_C_dep = (1- sr_object.Calculate_ReflectionCoefficients('C','C',Energy_particles[i]))*C_C[surface_index,-1]*Flux_C
        
        
        Gamma_C_reflected = (sr_object.Calculate_ReflectionCoefficients('C','C',Energy_particles[i]))*C_C[surface_index,-1]*Flux_C
        print(sr_object.Calculate_ReflectionCoefficients('C','C',Energy_particles[i]))
        
        
        Gamma_C_redep = beta_C*Gamma_C_ero
        Gamma_C_dep = Gamma_C_dep + Gamma_C_redep
        Gamma_C_ero = (1-beta_C)*Gamma_C_ero    
        
        # print("-----")
        # print("Incident Flux",1.0,"   Erosion flux: ",(Gamma_C_ero)/Flux_C)
        # print("Incident Flux",1.0,"   Deposition flux: ",(Gamma_C_dep)/Flux_C)
        # print("Incident Flux",1.0,"   Reflected Flux: ",(Gamma_C_reflected)/Flux_C)
        # print("Incident Flux",1.0,"   Erosion+Deposition+Reflection: ",(Gamma_C_dep+Gamma_C_ero+Gamma_C_reflected)/Flux_C)
        # print("-----")
            
        Gamma_C_ero_global[surface_index] = Gamma_C_ero_global[surface_index] + Gamma_C_ero 
        
        Gamma_C_dep_global[surface_index] = Gamma_C_dep_global[surface_index] + Gamma_C_dep 
        


#%%

# The following arrays will keep track of entries to be made in the next GITR run. 

vx_array = np.zeros(0)
vy_array = np.zeros(0)
vz_array = np.zeros(0)

x_array = np.zeros(0)
y_array = np.zeros(0)
z_array = np.zeros(0)

nP_global = 0 #tracks total number of eroded particles

new_entry_C = np.zeros((len(x1),1)) 
Gamma_C_net_global = np.zeros((len(x1),1)) 
Gamma_C_bulk_global = np.zeros((len(x1),1)) # without the concentration term

for surface_index in range(len(x1)):
    
    Gamma_C_bulk = 0
    Gamma_W_bulk = 0   
    
    Gamma_C_net_global[surface_index] = Gamma_C_dep_global[surface_index] - Gamma_C_ero_global[surface_index]
    
    Gamma_W_net = 0 # for the present case
    
    if (Gamma_C_net_global[surface_index] + Gamma_W_net) > 0: #deposition regime
        print("Deposition")
        Gamma_C_bulk_global[surface_index] = (Gamma_C_net_global[surface_index] + Gamma_W_net)        
        #Gamma_W_bulk = C_W[time]*(Gamma_C_net + Gamma_W_net)
    elif (Gamma_C_net_global[surface_index] + Gamma_W_net) < 0:#erosion regime
        print("erosion")
        Gamma_C_bulk_global[surface_index] = 0
        #Gamma_W_bulk = (Gamma_C_net + Gamma_W_net)
    
    #no_of_C = int(math.ceil(Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/weight_gitr))
    no_of_C = int(math.floor(Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/weight_gitr))
    no_of_C_frac = (Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/weight_gitr)%no_of_C
    
    if (np.random.uniform(low = 0.0, high = 1.0) < no_of_C_frac):
        no_of_C += 1
        
    #print(no_of_C)

    nP_global = nP_global + no_of_C
    
    if no_of_C > 0:
        
        p_C = ParticleDistribution()
        p_C.SetAttr('Np', no_of_C)
        
        # if min(x1[surface_index],x2[surface_index])<max(x1[surface_index],x2[surface_index]):
        #     p_C.SetAttr('x','Uniform',xmin=min(x1[surface_index],x2[surface_index]),xmax=max(x1[surface_index],x2[surface_index]))
        # else:
        #     p_C.SetAttr('x',np.full((no_of_C,1),x1[surface_index]))
            
        # if min(y1[surface_index],y2[surface_index])<max(y1[surface_index],y2[surface_index]):
        #     p_C.SetAttr('y','Uniform',xmin=min(y1[surface_index],y2[surface_index]),xmax=max(y1[surface_index],y2[surface_index]))
        # else:
        #     p_C.SetAttr('y',np.full((no_of_C,1),y1[surface_index]))
         
        # if min(z1[surface_index],z2[surface_index])<max(z1[surface_index],z2[surface_index]):
        #     p_C.SetAttr('z','Uniform',xmin=min(z1[surface_index],z2[surface_index]),xmax=max(z1[surface_index],z2[surface_index]))
        # else:
        #     p_C.SetAttr('z',np.full((no_of_C,1),z1[surface_index]))
        
        centroid_x = (x1[surface_index] + x2[surface_index] +x3[surface_index])/2
        centroid_y = (y1[surface_index] + y2[surface_index] +y3[surface_index])/2
        centroid_z = (z1[surface_index] + z2[surface_index] +z3[surface_index])/2
        
        p_C.SetAttr('x',np.full((no_of_C,1),centroid_x))
        p_C.SetAttr('y',np.full((no_of_C,1),centroid_y))
        p_C.SetAttr('z',np.full((no_of_C,1),centroid_z))
        
        p_C.SetAttr('vx',np.full((no_of_C,1),0.0))
        p_C.SetAttr('vy',np.full((no_of_C,1),0.0))
        
        p_C.SetAttr('vz','Thomson')
        E_C = p_C.Particles['vz'][:]
        
        E_C = np.sqrt(2*E_C*units.eV/(amu*units.mp))   #check these
        #print(E_C)
        p_C.SetAttr('vz',E_C)
            
        vx_array = np.append(vx_array,p_C.Particles['vx'])
        vy_array = np.append(vy_array,p_C.Particles['vy'])
        vz_array = np.append(vz_array,p_C.Particles['vz'])
        
        x_array = np.append(x_array,p_C.Particles['x'])
        y_array = np.append(y_array,p_C.Particles['y'])
        z_array = np.append(z_array,p_C.Particles['z'])
        
#%%
# Evolution of C_C

Time = Delta_t
Time_steps = 1e6
Delta_Time = Delta_t/Time_steps

old_C = np.reshape(C_C[:,-1],(8,1))
new_entry_C[:] =  old_C[:] + Delta_Time*(Gamma_C_net_global[:] - old_C[:]*Gamma_C_bulk_global[:])/(Delta_implant*n_atom)
for t in range(2,int(Time_steps+1)):
    new_entry_C[:] = new_entry_C[:] + Delta_Time*(Gamma_C_net_global[:] - new_entry_C[:]*Gamma_C_bulk_global[:])/(Delta_implant*n_atom)

#%%


C_C = np.concatenate((C_C,new_entry_C),axis=1)

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t)

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

#%%
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




#run the next simulation
    
    