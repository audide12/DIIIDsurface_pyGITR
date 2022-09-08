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

#Flux_H = 1.43222e+19 # Obtained from SOLPS simulation (Zack)

Flux_H = 1.0e20
alpha_c = 0.02       # Carbon concentration in the background plasma
Flux_C = alpha_c*Flux_H

Delta_t = 0.1 # in seconds
Delta_t_gitr = 1e-7
Delta_implant = 1e-5 # enter parameter value and units
amu_C = 12 #for carbon
amu_W = 74 #for tungsten

n_atom = 6e22 # average number density
weight_gitr = Delta_t/Delta_t_gitr
#weight_gitr = 1



Sputtering_yield_H_to_C = 0.005  # max 0.005
Sputtering_yield_C_to_C = 0.22   # 0.22 treated almost constant

Sputtering_yield_H_to_W = 0.002  # max 0.002 
Sputtering_yield_C_to_W = 0.5  # 0.5 treated almost constant

Reflection_yield_C_to_C = 0.005    # max 0.005 min 0.001  steady-state :  0.9
Reflection_yield_C_to_W = 0.75  # max 0.95 min 0.67    steady-state :  0.005

#Calculate Edist_H

#%%
# Reading position files of Carbon

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output_C/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit_C = np.array(PositionData['surfaceHit'])


surface_vx_C = np.array(PositionData['vx'])
surface_vy_C = np.array(PositionData['vy'])
surface_vz_C = np.array(PositionData['vz'])


Energy_particles_C = np.array(0.5*amu_C*1.66e-27*(surface_vx_C**2 + surface_vy_C**2 + surface_vz_C**2)/1.602e-19) # make sure that this energy is in eV


count_C = 0
count_C_6 = 0

for i in surfacehit_C:
    if i != -1:
        count_C+=1
    if i == 6:
        count_C_6+=1
print(count_C)
print(count_C_6)

#weight_flux_C = token_flux/count_C # this is the delta Gamma in the pdf

#print(surface_vx)

#Angle_particles = add this functionality

#%%

# Reading position files of Tungsten

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output_W/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit_W = np.array(PositionData['surfaceHit'])


surface_vx_W = np.array(PositionData['vx'])
surface_vy_W = np.array(PositionData['vy'])
surface_vz_W = np.array(PositionData['vz'])


Energy_particles_W = np.array(0.5*amu_W*1.66e-27*(surface_vx_W**2 + surface_vy_W**2 + surface_vz_W**2)/1.602e-19) # make sure that this energy is in eV


count_W = 0
for i in surfacehit_W:
    if i != -1:
        count_W+=1
print(count_W)

#weight_flux_W = token_flux/count_W

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

in_direction = np.array(config.geom.inDir)
plane_norm = np.array(config.geom.plane_norm)

#%%

# This section is to initialize the surface_evolution netcdf file

C_C = np.full((len(x1),1),0.0)
C_W = np.full((len(x1),1),0.0)

for k in range(len(x1)):
    if k<4:
        C_C[k] = 1.0
    else:
        C_W[k] = 1.0



N_GITR = 10000 # number of GITR particles
initial_token_flux = 1.0e19  # tunable parameter
Flux_proportionality_C = 0
Flux_proportionality_W = 0

for k in range(len(x1)):
    Flux_proportionality_C = Flux_proportionality_C + initial_token_flux*area[k]*Delta_t_gitr/N_GITR
    Flux_proportionality_W = Flux_proportionality_W + initial_token_flux*area[k]*Delta_t_gitr/N_GITR
    

Surface_time = np.full((1,1),0.0)

Surface_number = np.array(range(len(x1)))

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_W.nc', 'w', format='NETCDF4')
s_number_dim = ncFile.createDimension('surface_dim', len(x1)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))
s_concentration_C = ncFile.createVariable('surface_concentration_C',np.float64,('surface_dim','time_dim'))
s_concentration_W = ncFile.createVariable('surface_concentration_W',np.float64,('surface_dim','time_dim'))
flux_proportionality_C = ncFile.createVariable('Flux_Conversion_C',np.float64,('time_dim'))
flux_proportionality_W = ncFile.createVariable('Flux_Conversion_W',np.float64,('time_dim'))


s_number[:] = np.linspace(1,len(x1),len(x1))
s_time[:] = Surface_time
s_concentration_C[:,:] = C_C
s_concentration_W[:,:] = C_W
flux_proportionality_C[:] = Flux_proportionality_C
flux_proportionality_W[:] = Flux_proportionality_W

ncFile.close()


#%%

#Reading the surface features from the surface evolution netcdf file

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_W.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")


C_C = SurfaceConcentrationData['surface_concentration_C'][:,:]
C_W = SurfaceConcentrationData['surface_concentration_W'][:,:]
Flux_proportionality_C = SurfaceConcentrationData['Flux_Conversion_C'][:]
Flux_proportionality_W = SurfaceConcentrationData['Flux_Conversion_W'][:]

Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]

counter = len(Surface_time)

#%%
# Calculation of erosion and deposition fluxes for Carbon and Tungsten for each GITRb particle

Gamma_C_redep = np.zeros((len(x1),1))
Y_CW_Gamma_C_redep = np.zeros((len(x1),1))
Y_CC_Gamma_C_redep = np.zeros((len(x1),1))
Y_WC_Gamma_C_redep = np.zeros((len(x1),1))



for i in range(len(Energy_particles_C)):
    if surfacehit_C[i] != -1:
        surface_index = int(surfacehit_C[i])
        sr_object = Sputtering_and_reflection()
        
        Flux_C_local = Flux_proportionality_C[-1]/(Delta_t_gitr*area[surface_index])
        
        Gamma_C_redep[surface_index] = Gamma_C_redep[surface_index] + Flux_C_local  # check this
        Y_CW_Gamma_C_redep[surface_index] = Y_CW_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','W',Energy_particles_C[i])*Flux_C_local
        Y_CC_Gamma_C_redep[surface_index] = Y_CC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i])*Flux_C_local
        Y_WC_Gamma_C_redep[surface_index] = Y_WC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('W','C',Energy_particles_C[i])*Flux_C_local

Gamma_W_redep = np.zeros((len(x1),1))
Y_WW_Gamma_W_redep = np.zeros((len(x1),1))
        
for i in range(len(Energy_particles_W)):
    if surfacehit_W[i] != -1:
        surface_index = int(surfacehit_W[i])
        sr_object = Sputtering_and_reflection()
        
        Flux_W_local = Flux_proportionality_W[-1]/(Delta_t_gitr*area[surface_index])
        
        Gamma_W_redep[surface_index] = Gamma_W_redep[surface_index] + Flux_W_local  # check this
        Y_WW_Gamma_W_redep[surface_index] = Y_WW_Gamma_W_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('W','W',Energy_particles_W[i])*Flux_W_local
    
             
        
chi_W_ero =  Y_WW_Gamma_W_redep + Y_CW_Gamma_C_redep + Sputtering_yield_H_to_W*Flux_H + Sputtering_yield_C_to_W*Flux_C    

chi_C_ero_1 =  Y_CC_Gamma_C_redep + Sputtering_yield_H_to_C*Flux_H + Sputtering_yield_C_to_C*Flux_C   
chi_C_ero_2 =  Y_WC_Gamma_C_redep

chi_C_dep_1 = np.zeros((len(x1),1)) + (1-Reflection_yield_C_to_C)*Flux_C
chi_C_dep_2 = np.zeros((len(x1),1)) + (1-Reflection_yield_C_to_W)*Flux_C
   
Gamma_C_ero_global = np.reshape(C_C[:,-1],(len(x1),1)) * chi_C_ero_1 + np.reshape(C_W[:,-1],(len(x1),1)) * chi_C_ero_2
Gamma_C_dep_global = np.reshape(C_C[:,-1],(len(x1),1)) * chi_C_dep_1 + np.reshape(C_W[:,-1],(len(x1),1)) * chi_C_dep_2 + Gamma_C_redep
Gamma_W_ero_global = np.reshape(C_W[:,-1],(len(x1),1)) * chi_W_ero
Gamma_W_dep_global = Gamma_W_redep

#%%

# The following arrays will keep track of entries to be made in the next GITR run. 

vx_C_array = np.zeros(0)
vy_C_array = np.zeros(0)
vz_C_array = np.zeros(0)

x_C_array = np.zeros(0)
y_C_array = np.zeros(0)
z_C_array = np.zeros(0)

vx_W_array = np.zeros(0)
vy_W_array = np.zeros(0)
vz_W_array = np.zeros(0)

x_W_array = np.zeros(0)
y_W_array = np.zeros(0)
z_W_array = np.zeros(0)

nP_C_global = 0 #tracks total number of eroded particles
nP_W_global = 0 #tracks total number of eroded particles

prop_W = 0
prop_C = 0

for k in range(len(x1)):
    prop_W = prop_W + Gamma_W_ero_global[k]*area[k]*Delta_t_gitr
    prop_C = prop_C + Gamma_C_ero_global[k]*area[k]*Delta_t_gitr

prop_W = prop_W/N_GITR
prop_C = prop_C/N_GITR

#%%    
for surface_index in range(len(x1)):
    
    
    #no_of_C = int(math.floor(Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_C))
    no_of_C = round(np.array(Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_C).item())
    # no_of_C_frac = (Gamma_C_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_C)%no_of_C
    
    # if (np.random.uniform(low = 0.0, high = 1.0) < no_of_C_frac):
    #     no_of_C += 1
        
    #print(no_of_C)

    nP_C_global = nP_C_global + no_of_C
    
    
    if no_of_C > 0:
        
        
        # print("No of C: ",no_of_C)
        # print("Surface: ",surface_index)
        
        
        p_C = ParticleDistribution()
        p_C.SetAttr('Np', no_of_C)
        
        centroid_x = (x1[surface_index] + x2[surface_index] +x3[surface_index])/3
        centroid_y = (y1[surface_index] + y2[surface_index] +y3[surface_index])/3
        centroid_z = (z1[surface_index] + z2[surface_index] +z3[surface_index])/3
        
        p_C.SetAttr('x',np.full((no_of_C,1),centroid_x))
        p_C.SetAttr('y',np.full((no_of_C,1),centroid_y))
        p_C.SetAttr('z',np.full((no_of_C,1),centroid_z))
        
        p_C.SetAttr('vx',np.full((no_of_C,1),0.0))
        p_C.SetAttr('vy',np.full((no_of_C,1),0.0))
        
        p_C.SetAttr('vz','Thomson')
        E_C = p_C.Particles['vz'][:]
        
        E_C = np.sqrt(2*E_C*units.eV/(amu_C*units.mp))   #check these
        #print(E_C)
        p_C.SetAttr('vz',E_C)
        
        if (surface_index<4):
            sign = 1
        else:
            sign = -1
            
        p_C.ScaleAttr('vz',sign)    
            
        vx_C_array = np.append(vx_C_array,p_C.Particles['vx'])
        vy_C_array = np.append(vy_C_array,p_C.Particles['vy'])
        vz_C_array = np.append(vz_C_array,p_C.Particles['vz'])
        
        #print("Surface Index: ", surface_index,"   and vz:  ", sign*p_C.Particles['vz'])
        
        x_C_array = np.append(x_C_array,p_C.Particles['x'])
        y_C_array = np.append(y_C_array,p_C.Particles['y'])
        z_C_array = np.append(z_C_array,p_C.Particles['z'])
        
        
    #no_of_W = int(math.floor(Gamma_W_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_W))    
    no_of_W = round(np.array(Gamma_W_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_W).item())
    # no_of_W_frac = (Gamma_W_ero_global[surface_index]*Delta_t_gitr*area[surface_index]/prop_W)%no_of_W
    
    # if (np.random.uniform(low = 0.0, high = 1.0) < no_of_W_frac):
    #     no_of_W += 1
        
    #print(no_of_C)

    nP_W_global = nP_W_global + no_of_W
    
    if no_of_W > 0:
        
        
        # print("No of C: ",no_of_C)
        # print("Surface: ",surface_index)
        
        
        p_W = ParticleDistribution()
        p_W.SetAttr('Np', no_of_W)
        
        centroid_x = (x1[surface_index] + x2[surface_index] +x3[surface_index])/3
        centroid_y = (y1[surface_index] + y2[surface_index] +y3[surface_index])/3
        centroid_z = (z1[surface_index] + z2[surface_index] +z3[surface_index])/3
        
        p_W.SetAttr('x',np.full((no_of_W,1),centroid_x))
        p_W.SetAttr('y',np.full((no_of_W,1),centroid_y))
        p_W.SetAttr('z',np.full((no_of_W,1),centroid_z))
        
        p_W.SetAttr('vx',np.full((no_of_W,1),0.0))
        p_W.SetAttr('vy',np.full((no_of_W,1),0.0))
        
        p_W.SetAttr('vz','Thomson')
        E_W = p_W.Particles['vz'][:]
        
        E_W = np.sqrt(2*E_W*units.eV/(amu_W*units.mp))   #check these
        #print(E_C)
        p_W.SetAttr('vz',E_W)
                
        if (surface_index<4):
            sign = 1
        else:
            sign = -1
            
        p_W.ScaleAttr('vz',sign)      
            
        vx_W_array = np.append(vx_W_array,p_W.Particles['vx'])
        vy_W_array = np.append(vy_W_array,p_W.Particles['vy'])
        vz_W_array = np.append(vz_W_array,p_W.Particles['vz'])
        
        x_W_array = np.append(x_W_array,p_W.Particles['x'])
        y_W_array = np.append(y_W_array,p_W.Particles['y'])
        z_W_array = np.append(z_W_array,p_W.Particles['z'])
        
        
#%%
# Writing the particle list for the next GITR run  
p_C = ParticleDistribution()



p_C.SetAttr('Np', nP_C_global)

# Set positions of particles
p_C.SetAttr('x',x_C_array)
p_C.SetAttr('y',y_C_array)
p_C.SetAttr('z',z_C_array)

# Set velocities of particles
p_C.SetAttr('vx',vx_C_array)
p_C.SetAttr('vy',vy_C_array)
p_C.SetAttr('vz',vz_C_array)

# Writing the particle source file for the next GITR run.
p_C.WriteParticleFile('particleConf_C.nc')


# the same for Tungsten
p_W = ParticleDistribution()

p_W.SetAttr('Np', nP_W_global)

# Set positions of particles
p_W.SetAttr('x',x_W_array)
p_W.SetAttr('y',y_W_array)
p_W.SetAttr('z',z_W_array)

# Set velocities of particles
p_W.SetAttr('vx',vx_W_array)
p_W.SetAttr('vy',vy_W_array)
p_W.SetAttr('vz',vz_W_array)

# Writing the particle source file for the next GITR run.
p_W.WriteParticleFile('particleConf_W.nc')

#run the next simulation        
        
#%%


# The actual surface model differential equation

# Evolution of C_C and C_W

Time = Delta_t
Time_steps = 1e4
Delta_Time = Delta_t/Time_steps

new_entry_C = np.reshape(C_C[:,-1],(len(x1),1))

new_entry_W = np.reshape(C_W[:,-1],(len(x1),1))
#new_entry_C = np.zeros(len(x1))

#new_entry_C[:] =  old_C[:] + Delta_Time*(Gamma_C_net_global[:] - old_C[:]*Gamma_C_bulk_global[:])/(Delta_implant*n_atom)

        
for t in range(1,int(Time_steps)):
    
    Gamma_W_ero = new_entry_W*chi_W_ero
    
    Gamma_C_ero = new_entry_C*chi_C_ero_1 + new_entry_W*chi_C_ero_2
    
    Gamma_C_dep = new_entry_C*chi_C_dep_1 + new_entry_W*chi_C_dep_2 + Gamma_C_redep 
    
    Gamma_W_dep = Gamma_W_redep
    
    
    # determining erosion or deposition
    Gamma_C_net = Gamma_C_dep - Gamma_C_ero
    
    Gamma_W_net = -Gamma_W_ero
    
    Gamma_C_bulk = np.zeros((len(x1),1))
    Gamma_W_bulk = np.zeros((len(x1),1))
    
    #print(Gamma_C_net)
    
    for surface_index in range(len(x1)):
        if (Gamma_C_net[surface_index] + Gamma_W_net[surface_index]) > 0: # deposition regime
            #print("deposition")
            Gamma_C_bulk[surface_index] = new_entry_C[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
            Gamma_W_bulk[surface_index] = new_entry_W[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
        
        else:  #  erosion regime
            #print("erosion")
            Gamma_C_bulk[surface_index] = 0
            Gamma_W_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
    
    #print(t)
    new_entry_C = new_entry_C + Delta_Time*(Gamma_C_net - Gamma_C_bulk)/(Delta_implant*n_atom)
    
    new_entry_W = new_entry_W + Delta_Time*(Gamma_W_net - Gamma_W_bulk)/(Delta_implant*n_atom)
        

#%%


C_C = np.concatenate((C_C,new_entry_C),axis=1)
C_W = np.concatenate((C_W,new_entry_W),axis=1)

Flux_proportionality_C = np.append(Flux_proportionality_C,(1/prop_C))
Flux_proportionality_W = np.append(Flux_proportionality_W,(1/prop_W))

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t)

#Writing the surface features with time

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_W_new.nc', 'w', format='NETCDF4')
s_number_dim = ncFile.createDimension('surface_dim', len(x1)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))
s_concentration_C = ncFile.createVariable('surface_concentration_C',np.float64,('surface_dim','time_dim'))
s_concentration_W = ncFile.createVariable('surface_concentration_W',np.float64,('surface_dim','time_dim'))
flux_proportionality_C = ncFile.createVariable('Flux_Conversion_C',np.float64,('time_dim'))
flux_proportionality_W = ncFile.createVariable('Flux_Conversion_W',np.float64,('time_dim'))


s_number[:] = np.linspace(1,len(x1),len(x1))
s_time[:] = Surface_time
s_concentration_C[:,:] = C_C
s_concentration_W[:,:] = C_W
flux_proportionality_C[:] = Flux_proportionality_C
flux_proportionality_W[:] = Flux_proportionality_W

ncFile.close()

#%%

    
    