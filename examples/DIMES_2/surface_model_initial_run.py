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


from pyGITR.process import *
from pyGITR.process_functions import *
from surface_model_functions import *

from pyGITR.process import *
from pyGITR.process_functions import *
from surface_model_functions import *

from pyGITR.Physical_Sputtering import *
from pyGITR.make_particleSource import *



#%%
# Reading position files of Carbon

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/output_C_1/positions.nc'

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


FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/output_Si_1/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit_Si = np.array(PositionData['surfaceHit'])
surface_vx_Si = np.array(PositionData['vx'])
surface_vy_Si = np.array(PositionData['vy'])
surface_vz_Si = np.array(PositionData['vz'])


Energy_particles_Si = np.array(0.5*amu_Si*1.66e-27*(surface_vx_Si**2 + surface_vy_Si**2 + surface_vz_Si**2)/1.602e-19) # make sure that this energy is in eV

count_Si = 0
for i in surfacehit_Si:
    if i != -1:
        count_Si+=1
print(count_Si,"have hit a mesh element (not necessarily a surface)")



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


Concentration = {}
for z in Zs:
    Concentration[z] = np.full((len(Surfaces),1), 0.0)
    # step through the surfaces, set to 1.0 for those surfaces with Z
    for k,surface in enumerate(Surfaces):
        if Z[surface] == z:
            Concentration[z][k] = 1.0


# Create initial ncFile

#os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_Si.nc")

os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_W.nc")

makeInitNC(len(Surfaces),area,Concentration)


#%%

#Reading the surface features from the surface evolution netcdf file

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_Si.nc'

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
Y_CSiC_Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CSi_Gamma_C_redep = np.zeros((len(Surfaces),1)) # wasn't needed
Y_CC_Gamma_C_redep = np.zeros((len(Surfaces),1))


for i in range(len(Energy_particles_C)):
    if surfacehit_C[i] != -1:

        surface_index = int(surfacehit_C[i])
        sr_object = Sputtering_and_reflection()

        for j in Surfaces:
            if j == surface_index:

                #print("yes")
                Flux_C_local = Flux_proportionality[6][-1]/(Delta_t_gitr*area[surface_index])  # we start with initial weight 1 (uniform)
                #print(Flux_C_local)
                
                Gamma_C_redep[surface_index] = Gamma_C_redep[surface_index] + Flux_C_local  # check this
                Y_CSiC_Gamma_C_redep[surface_index] = Y_CSiC_Gamma_C_redep[surface_index] + 0.5*sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i])*Flux_C_local
                Y_CC_Gamma_C_redep[surface_index] = Y_CC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i])*Flux_C_local
                

Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_SiSi_Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_SiSiC_Gamma_Si_redep = np.zeros((len(Surfaces),1))
        
        
for i in range(len(Energy_particles_Si)):
    if surfacehit_Si[i] != -1:
        surface_index = int(surfacehit_Si[i])

        sr_object = Sputtering_and_reflection()
        
        for j in Surfaces:
            if j == surface_index:

                Flux_Si_local = Flux_proportionality[14][-1]/(Delta_t_gitr*area[surface_index]) # we start with initial weight 1 (uniform)
                
                Gamma_Si_redep[surface_index] = Gamma_Si_redep[surface_index] + Flux_Si_local  # check this
                Y_SiSi_Gamma_Si_redep[surface_index] = Y_SiSi_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','Si',Energy_particles_Si[i])*Flux_Si_local
                Y_SiSiC_Gamma_Si_redep[surface_index] = Y_SiSiC_Gamma_Si_redep[surface_index] + 0.5*sr_object.Calculate_PhysicalSputteringParameters('Si','Si',Energy_particles_Si[i])*Flux_Si_local
                

        
beta_eroC1 = np.zeros((len(Surfaces),1)) + Sputtering_yield_H_to_C*Flux_H
beta_eroC2 = np.zeros((len(Surfaces),1)) + Sputtering_yield_C_to_C*Flux_C
beta_C_dep = np.zeros((len(Surfaces),1)) + (1-Reflection_yield_C_to_C)*Flux_C    
beta_SiC1 =  np.zeros((len(Surfaces),1)) + Sputtering_yield_H_to_SiC*Flux_H + Sputtering_yield_C_to_SiC*Flux_C
beta_eroSi1 = np.zeros((len(Surfaces),1)) + Sputtering_yield_H_to_Si*Flux_H


beta_SiC = beta_SiC1 + Y_CSiC_Gamma_C_redep + Y_SiSiC_Gamma_Si_redep
beta_eroSi2 = Y_SiSi_Gamma_Si_redep
beta_eroC3 = Y_CC_Gamma_C_redep
   
Gamma_SiC_ero_global = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))*beta_SiC

#print(sum(Gamma_SiC_ero_global))

Gamma_Si_ero_global = Gamma_SiC_ero_global + np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*beta_eroSi1 + np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*beta_eroSi2
Gamma_Si_ero_exclusive = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*beta_eroSi1 + np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*beta_eroSi2
Gamma_C_ero_global = Gamma_SiC_ero_global + np.reshape(Concentration[6][:,-1],(len(Surfaces),1))*(beta_eroC1 + beta_eroC2 + beta_eroC3)
Gamma_C_ero_exclusive = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))*(beta_eroC1 + beta_eroC2 + beta_eroC3)
Gamma_C_dep_global = Gamma_C_redep + np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*Flux_C + np.reshape(Concentration[6][:,-1],(len(Surfaces),1))*beta_C_dep +  np.reshape(Concentration[20][:,-1],(len(Surfaces),1))*beta_C_dep     
Gamma_Si_dep_global = Gamma_Si_redep




#%%

# The following arrays will keep track of entries to be made in the next GITR run. 

nP_C_global = 0 #tracks total number of eroded particles
nP_Si_global = 0 #tracks total number of eroded particles

prop_Si = 0
prop_C = 0
prop_SiC = 1

for surface_index in range(len(Surfaces)):
    prop_Si = prop_Si + Gamma_Si_ero_global[surface_index]*area[surface_index]*Delta_t_gitr
    prop_C = prop_C + Gamma_C_ero_global[surface_index]*area[surface_index]*Delta_t_gitr

prop_Si = prop_Si/N_GITR
prop_C = prop_C/N_GITR


particleSourceDict_C = {}
for i,surface in enumerate(Surfaces):

    num_particles = np.array(Gamma_C_ero_global[i]*Delta_t_gitr*area[surface]/prop_C).item()
    if num_particles!=0: 
        print("Surface:",surface,"C particles:",num_particles)
        particleSourceDict_C[surface] = round(num_particles)
        nP_C_global+=round(num_particles)
            
particleSourceDict_Si = {}
for i,surface in enumerate(Surfaces):
    num_particles = np.array(Gamma_Si_ero_global[i]*Delta_t_gitr*area[surface]/prop_Si).item()
    if num_particles!=0: 
        print("Surface:",surface,"Si particles:",num_particles)
        particleSourceDict_Si[surface] = round(num_particles)
        nP_Si_global+=round(num_particles)            


makeParticleSource(particleSourceDict_C, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_C.nc","/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/weights_C.nc")

if nP_Si_global>0:
    
    makeParticleSource(particleSourceDict_Si, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_Si.nc","/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/weights_Si.nc")
            
          
#%%
# Estimating the total time evolution for the surface model

last_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))
last_entry_Si = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))
last_entry_SiC = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))


Gamma_C_net = Gamma_C_dep_global - Gamma_C_ero_exclusive

Gamma_Si_net = Gamma_Si_dep_global - Gamma_Si_ero_exclusive

Gamma_SiC_net = (-1)*Gamma_SiC_ero_global


Gamma_C_bulk = np.zeros((len(Surfaces),1))
Gamma_Si_bulk = np.zeros((len(Surfaces),1))
Gamma_SiC_bulk = np.zeros((len(Surfaces),1))

#print(Gamma_C_net)

for surface_index in range(len(Surfaces)):

    if (Gamma_C_net[surface_index] + Gamma_Si_net[surface_index] + Gamma_SiC_net[surface_index]) > 0: # deposition regime
        print("deposition")
        Gamma_C_bulk[surface_index] = last_entry_C[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
        Gamma_Si_bulk[surface_index] = last_entry_Si[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
        Gamma_SiC_bulk[surface_index] = last_entry_SiC[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
    
    else:  #  erosion regime
        print("erosion")
        Gamma_C_bulk[surface_index] = 0.0
        Gamma_Si_bulk[surface_index] = 0.0
        Gamma_SiC_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
    

RHS_C = Gamma_C_net - Gamma_C_bulk
RHS_Si = Gamma_Si_net - Gamma_Si_bulk
RHS_SiC = Gamma_SiC_net - Gamma_SiC_bulk

Stopping_criteria = 0.1 # for C_C and C_W

RHS_C   = np.abs(RHS_C) 
RHS_Si  = np.abs(RHS_Si)
RHS_SiC = np.abs(RHS_SiC)
       
Delta_t_surface_estimate_C = (Delta_implant*n_atom*Stopping_criteria)/RHS_C

Delta_t_surface_estimate_Si = (Delta_implant*n_atom*Stopping_criteria)/RHS_Si

Delta_t_surface_estimate_SiC = (Delta_implant*n_atom*Stopping_criteria)/RHS_SiC

Delta_t_surface = min(np.amin(Delta_t_surface_estimate_C),np.amin(Delta_t_surface_estimate_Si),np.amin(Delta_t_surface_estimate_SiC))        



#%%
# The actual surface model differential equation
# Evolution of C_C and C_Si
# Stopping criterion implemented
# Delta_t is a constant for the surface model

Time = Delta_t_surface
Time_steps = 1e4
Delta_Time = Delta_t/Time_steps
Delta_t_Stopping = 0

Stopping_criteria = 0.2 # for C_C and C_Si

new_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))   # populating it with last entry
new_entry_Si = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))
new_entry_SiC = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))


Stopping_criteria = 0.1 # for C_C and C_Si

        
for t in range(1,int(Time_steps)):
    

    
    Gamma_SiC_ero = beta_SiC * new_entry_SiC
    Gamma_Si_ero = Gamma_SiC_ero + beta_eroSi1*new_entry_Si + beta_eroSi2*new_entry_Si
    Gamma_Si_ero_exclusive = beta_eroSi1*new_entry_Si + beta_eroSi2*new_entry_Si
    
    Gamma_C_ero = Gamma_SiC_ero + beta_eroC1*new_entry_C + beta_eroC2*new_entry_C + beta_eroC3*new_entry_C
    Gamma_C_ero_exclusive = beta_eroC1*new_entry_C + beta_eroC2*new_entry_C + beta_eroC3*new_entry_C
    
    Gamma_C_dep = new_entry_Si*Flux_C + new_entry_C*beta_C_dep + Gamma_C_redep + new_entry_SiC*beta_C_dep
    Gamma_Si_dep = Gamma_Si_redep
    
    
    # determining erosion or deposition
    Gamma_C_net = Gamma_C_dep - Gamma_C_ero_exclusive
    
    Gamma_Si_net = Gamma_Si_dep - Gamma_Si_ero_exclusive
    
    Gamma_SiC_net = - Gamma_SiC_ero
    
    Gamma_C_bulk = np.zeros((len(Surfaces),1))
    Gamma_Si_bulk = np.zeros((len(Surfaces),1))
    Gamma_SiC_bulk = np.zeros((len(Surfaces),1))

    
    for surface_index in range(len(Surfaces)):

        if (Gamma_C_net[surface_index] + Gamma_Si_net[surface_index] + Gamma_SiC_net[surface_index]) > 0: # deposition regime
            #print("deposition")
            Gamma_C_bulk[surface_index] = new_entry_C[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
            Gamma_Si_bulk[surface_index] = new_entry_Si[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
            Gamma_SiC_bulk[surface_index] = new_entry_SiC[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
        
        else:  #  erosion regime
            #print("erosion")
            Gamma_C_bulk[surface_index] = 0.0
            Gamma_Si_bulk[surface_index] = 0.0
            Gamma_SiC_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
               
    
    #print(t)
    new_entry_C = new_entry_C + Delta_Time*(Gamma_C_net - Gamma_C_bulk)/(Delta_implant*n_atom)
    

    new_entry_Si = new_entry_Si + Delta_Time*(Gamma_Si_net - Gamma_Si_bulk)/(Delta_implant*n_atom)
    
    new_entry_SiC = new_entry_SiC + Delta_Time*(Gamma_SiC_net - Gamma_SiC_bulk)/(Delta_implant*n_atom)
    
    Delta_t_Stopping += Delta_Time
        
    if (np.abs(new_entry_C-last_entry_C)>Stopping_criteria).any() or (np.abs(new_entry_Si-last_entry_Si)>Stopping_criteria).any() or (np.abs(new_entry_SiC-last_entry_SiC)>Stopping_criteria).any():
        print(Delta_t_Stopping," Delta_t_Stopping ", t)
        break
        

#%%
# Appending time to all the surface characteristics

Concentration[6] = np.concatenate((Concentration[6],new_entry_C),axis=1)
Concentration[14] = np.concatenate((Concentration[14],new_entry_Si),axis=1)
Concentration[20] = np.concatenate((Concentration[20],new_entry_SiC),axis=1)


Flux_proportionality[6] = np.append(Flux_proportionality[6],prop_C)
Flux_proportionality[14] = np.append(Flux_proportionality[14],prop_Si)
Flux_proportionality[20] = np.append(Flux_proportionality[20],prop_SiC)

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t_Stopping)


#%%
#Writing the surface features with time


#Writing the surface features with time


os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_Si.nc")

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/surface_evolution_C_Si.nc', 'w', format='NETCDF4')


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

    
    