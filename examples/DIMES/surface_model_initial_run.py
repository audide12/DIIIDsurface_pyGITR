#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:12:31 2022

@author: de
"""


import electronvolt as units
import io, libconf
import numpy as np
import pyGITR
from pyGITR.Particles import *
import netCDF4
from netCDF4 import Dataset
import os
import math

#from pyGITR.process import *
#from pyGITR.process_functions import *
from surface_model_functions import *
from pyGITR.Physical_Sputtering import *
from pyGITR.make_particleSource import *



#%%
# Reading position files of Carbon

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output_C/positions.nc'
PositionData = netCDF4.Dataset(FileNameHistory)

surfacehit_C = np.array(PositionData['surfaceHit'])
surface_vx_C = np.array(PositionData['vx'])
surface_vy_C = np.array(PositionData['vy'])
surface_vz_C = np.array(PositionData['vz'])

Energy_particles_C = np.array(0.5*amu_C*1.66e-27*(surface_vx_C**2 + surface_vy_C**2 + surface_vz_C**2)/1.602e-19) # make sure that this energy is in eV

count_C = 0

for i in surfacehit_C:
    if i != -1:
        count_C+=1
    
print(count_C,"have hit a mesh element (not necessarily a surface)")

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
print(count_W,"have hit a mesh element (not necessarily a surface)")


#%%

#Reading geometry files

GeomFile = "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/gitrGeom.cfg"
x1,x2,x3,y1,y2,y3,z1,z2,z3,area,surf,Z,a,b,c,d,in_direction,plane_norm = getGeom(GeomFile)

# Initialize the surface_evolution netcdf file
# Only care about surfaces
Zs = []
Surfaces = []
idx = np.arange(0,len(surf))
for surface,z,i in zip(surf,Z,idx):
    if surface!=0:
        Zs.append(z)
        Surfaces.append(i)
Zs = np.unique(Zs)
print(Zs,"make up the", len(Surfaces),"surface mesh elements")
# print(Surfaces)

#%%

# Initiallize all surfaces with concentrations equal to the Z at that surface

Zs = np.append(Zs,74)# Add tungsten to the mix.

Concentration = {}
for z in Zs:
    Concentration[z] = np.full((len(Surfaces),1), 0.0)
    # step through the surfaces, set to 1.0 for those surfaces with Z
    for k,surface in enumerate(Surfaces):
        if Z[surface] == z:
            Concentration[z][k] = 1.0


# Create initial ncFile
os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_W.nc")
makeInitNC(len(Surfaces),area,Concentration)


#%%

#Reading the surface features from the surface evolution netcdf file
FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_W.nc'
SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")

# Record concentrations of all surface elements and their initial Z
Flux_proportionality = {}
for z in Zs:
    Concentration[z] = SurfaceConcentrationData['surface_concentration_{}'.format(z)][:,:]
    Flux_proportionality[z] = SurfaceConcentrationData['Flux_Conversion_{}'.format(z)][:]

Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)
#%%
# Calculation of erosion and deposition fluxes for Carbon and Tungsten for each GITRb particle

Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CW_Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CC_Gamma_C_redep = np.zeros((len(Surfaces),1))


for i in range(len(Energy_particles_C)):
    if surfacehit_C[i] != -1:
        surface_index = int(surfacehit_C[i])
        sr_object = Sputtering_and_reflection()

        for j in Surfaces:
            if j == surface_index:
                Flux_C_local = Flux_proportionality[6][-1]/(Delta_t_gitr*area[surface_index])
                
                Gamma_C_redep[surface_index] = Gamma_C_redep[surface_index] + Flux_C_local  # check this
                Y_CW_Gamma_C_redep[surface_index] = Y_CW_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','W',Energy_particles_C[i])*Flux_C_local
                Y_CC_Gamma_C_redep[surface_index] = Y_CC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i])*Flux_C_local
                

Gamma_W_redep = np.zeros((len(Surfaces),1))
Y_WW_Gamma_W_redep = np.zeros((len(Surfaces),1))
Y_WC_Gamma_W_redep = np.zeros((len(Surfaces),1))
        
for i in range(len(Energy_particles_W)):
    if surfacehit_W[i] != -1:
        surface_index = int(surfacehit_W[i])
        sr_object = Sputtering_and_reflection()
        
        for j in Surfaces:
            if j == surface_index:
                Flux_W_local = Flux_proportionality[74][-1]/(Delta_t_gitr*area[surface_index])
                
                Gamma_W_redep[surface_index] = Gamma_W_redep[surface_index] + Flux_W_local  # check this
                Y_WW_Gamma_W_redep[surface_index] = Y_WW_Gamma_W_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('W','W',Energy_particles_W[i])*Flux_W_local
                Y_WC_Gamma_W_redep[surface_index] = Y_WC_Gamma_W_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('W','C',Energy_particles_W[i])*Flux_W_local
                

        
        
chi_W_ero =  Y_WW_Gamma_W_redep + Y_CW_Gamma_C_redep + Sputtering_yield_H_to_W*Flux_H + Sputtering_yield_C_to_W*Flux_C    

chi_C_ero_1 =  Y_CC_Gamma_C_redep + Sputtering_yield_H_to_C*Flux_H + Sputtering_yield_C_to_C*Flux_C   
chi_C_ero_2 =  Y_WC_Gamma_W_redep

chi_C_dep_1 = np.zeros((len(Surfaces),1)) + (1-Reflection_yield_C_to_C)*Flux_C
chi_C_dep_2 = np.zeros((len(Surfaces),1)) + (1-Reflection_yield_C_to_W)*Flux_C
   
Gamma_C_ero_global = np.reshape(Concentration[6][:,-1],(len(Surfaces),1)) * chi_C_ero_1 + np.reshape(Concentration[74][:,-1],(len(Surfaces),1)) * chi_C_ero_2
Gamma_C_dep_global = np.reshape(Concentration[6][:,-1],(len(Surfaces),1)) * chi_C_dep_1 + np.reshape(Concentration[74][:,-1],(len(Surfaces),1)) * chi_C_dep_2 + Gamma_C_redep
Gamma_W_ero_global = np.reshape(Concentration[74][:,-1],(len(Surfaces),1)) * chi_W_ero
Gamma_W_dep_global = Gamma_W_redep

#%%

# The following arrays will keep track of entries to be made in the next GITR run. 

nP_C_global = 0 #tracks total number of eroded particles
nP_W_global = 0 #tracks total number of eroded particles

prop_W = 0
prop_C = 0

for surface_index in range(len(Surfaces)):
    prop_W = prop_W + Gamma_W_ero_global[surface_index]*area[surface_index]*Delta_t_gitr
    prop_C = prop_C + Gamma_C_ero_global[surface_index]*area[surface_index]*Delta_t_gitr

prop_W = prop_W/N_GITR
prop_C = prop_C/N_GITR


particleSourceDict_C = {}
for i,surface in enumerate(Surfaces):
    if Z[i] == 6:
        num_particles = np.array(Gamma_C_ero_global[i]*Delta_t_gitr*area[surface]/prop_C).item()
        if num_particles!=0: 
            print("Surface:",surface,"C particles:",num_particles)
            particleSourceDict_C[surface] = round(num_particles)
            
particleSourceDict_W = {}
for i,surface in enumerate(Surfaces):
    if Z[i] == 74:
        num_particles = np.array(Gamma_W_ero_global[i]*Delta_t_gitr*area[surface]/prop_W).item()
        if num_particles!=0: 
            print("Surface:",surface,"W particles:",num_particles)
            particleSourceDict_W[surface] = round(num_particles)            

makeParticleSource(particleSourceDict_C, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_C.nc")
#makeParticleSource(particleSourceDict_W, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_W.nc")


#%%

# Estimating the total time evolution for the surface model

last_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))
last_entry_W = np.reshape(Concentration[74][:,-1],(len(Surfaces),1))

Gamma_W_ero = last_entry_W*chi_W_ero
Gamma_C_ero = last_entry_C*chi_C_ero_1 + last_entry_W*chi_C_ero_2
Gamma_C_dep = last_entry_C*chi_C_dep_1 + last_entry_W*chi_C_dep_2 + Gamma_C_redep 
Gamma_W_dep = Gamma_W_redep

Gamma_C_net = Gamma_C_dep - Gamma_C_ero

Gamma_W_net = -Gamma_W_ero

Gamma_C_bulk = np.zeros((len(Surfaces),1))
Gamma_W_bulk = np.zeros((len(Surfaces),1))

#print(Gamma_C_net)

for surface_index in range(len(Surfaces)):
    if (Gamma_C_net[surface_index] + Gamma_W_net[surface_index]) > 0: # deposition regime
        #print("deposition")
        Gamma_C_bulk[surface_index] = last_entry_C[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
        Gamma_W_bulk[surface_index] = last_entry_W[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
    
    else:  #  erosion regime
        #print("erosion")
        Gamma_C_bulk[surface_index] = 0
        Gamma_W_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_W_net[surface_index])


RHS_C = Gamma_C_net - Gamma_C_bulk
RHS_W = Gamma_W_net - Gamma_W_bulk

Stopping_criteria = 0.1 # for C_C and C_W
        
Delta_t_surface_estimate_C = (Delta_implant*n_atom*Stopping_criteria)/RHS_C

Delta_t_surface_estimate_W = (Delta_implant*n_atom*Stopping_criteria)/RHS_W

Delta_t_surface = min(np.amin(Delta_t_surface_estimate_C),np.amin(Delta_t_surface_estimate_C))        


#%%
# The actual surface model differential equation
# Evolution of C_C and C_W
# Stopping criterion implemented
# Delta_t is a constant for the surface model

Time = Delta_t_surface
Time_steps = 1e4
Delta_Time = Delta_t/Time_steps
Delta_t_Stopping = 0
Stopping_criteria = 0.1 # for C_C and C_W

new_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))
new_entry_W = np.reshape(Concentration[74][:,-1],(len(Surfaces),1))
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
    
    Gamma_C_bulk = np.zeros((len(Surfaces),1))
    Gamma_W_bulk = np.zeros((len(Surfaces),1))
    
    #print(Gamma_C_net)
    
    for surface_index in range(len(Surfaces)):
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
    
    Delta_t_Stopping += Delta_Time
        
    if (np.abs(new_entry_C-last_entry_C)>Stopping_criteria).any() or (np.abs(new_entry_W-last_entry_W)>Stopping_criteria).any():
        print(Delta_t_Stopping," Delta_t_Stopping ", t)
        break
        

#%%
# Appending time to all the surface characteristics

Concentration[6] = np.concatenate((Concentration[6],new_entry_C),axis=1)
Concentration[74] = np.concatenate((Concentration[74],new_entry_W),axis=1)

Flux_proportionality[6] = np.append(Flux_proportionality[6],(1/prop_C))
Flux_proportionality[74] = np.append(Flux_proportionality[74],(1/prop_W))

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t_Stopping)


#%%
#Writing the surface features with time

os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_W.nc")

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C_W.nc', 'w', format='NETCDF4')
s_number_dim = ncFile.createDimension('surface_dim', len(Surfaces)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))


s_concentration = {}
flux_proportionality = {}
for z in Concentration.keys():
    s_concentration[z] = ncFile.createVariable('surface_concentration_{}'.format(z), np.float64, ('surface_dim','time_dim'))
    flux_proportionality[z] = ncFile.createVariable('Flux_Conversion_{}'.format(z),np.float64,('time_dim'))
    print(z)


s_number[:] = np.linspace(1,len(Surfaces),len(Surfaces))
s_time[:] = Surface_time


for z in Zs:
   s_concentration[z][:,:] = Concentration[z]
   flux_proportionality[z][:] = Flux_proportionality[z]

ncFile.close()

#%%

    
    