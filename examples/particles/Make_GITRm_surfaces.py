#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:15:20 2023

@author: de
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import Dataset

nE = 50
nA = 40
nS = 2   # 1 is for C and 2 is for W

energies = np.logspace(1,3,nE)
angles = np.linspace(0,89.5,nA)
species = np.array([1,2])  # 1 is for C and 2 is for W

interactions = [('C','C'),('C','W'),('W','C'),('W','W')]

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/surface_model_GITRm_C_W.nc', 'w', format='NETCDF4')

angle_dim = ncFile.createDimension('angle_dim', nA) # angle dimension
energy_dim = ncFile.createDimension('energy_dim', nE) # energy dimension
projectile_dim = ncFile.createDimension('projectile_dim', nS) # species dimension 
target_dim = ncFile.createDimension('target_dim', nS) # species dimension 

angle_number = ncFile.createVariable('Angles', np.float32, ('angle_dim',))
energy_number = ncFile.createVariable('Energies', np.float32, ('energy_dim',))
projectile_number = ncFile.createVariable('Projectiles', np.float32, ('projectile_dim',))  # 1 is for C and 2 is for W
target_number = ncFile.createVariable('Targets', np.float32, ('target_dim',))  # 1 is for C and 2 is for W


Species_Sputtering = ncFile.createVariable('Physical_Sputtering', np.float64, ('projectile_dim','target_dim','angle_dim','energy_dim'))    
   
angle_number[:] = angles
energy_number[:] = energies
projectile_number[:] = species
target_number[:] = species

Sputtering_yield = np.zeros((nS,nS,nA,nE)) 

for target_idx,target in enumerate(species):
    for projectile_idx,projectile in enumerate(species):
        for angle_idx,angle in enumerate(angles):
            for energy_idx,energy in enumerate(energies):
                sr_object = Sputtering_and_reflection()

                #print(angle, "    ", energy)
                # yld = sputtering_yield(incident, substrate, energy, angle, num_samples)
                if target == 1:
                    substrate = 'C'
                else:
                    substrate = 'W'
                if projectile == 1:
                    incident = 'C'
                else:
                    incident = 'W'    
                
                yld = sr_object.Calculate_PhysicalSputteringParameters(incident,substrate,energy,angle)
                
                if (np.isnan(yld)):
                    print("energy:  ",energy, "   and angle    ",angle,"   gives yield   ",yld)
                #print(yld)    
                    
                Sputtering_yield[projectile_idx,target_idx,angle_idx,energy_idx] = yld

Species_Sputtering[:,:,:,:] = Sputtering_yield             

ncFile.close()         

#%%
FileName='/Users/de/Research/DIIIDsurface_pyGITR/examples/surface_model_GITRm_C_W.nc'

s_model_data = Dataset(FileName, "r", format="NETCDF4")