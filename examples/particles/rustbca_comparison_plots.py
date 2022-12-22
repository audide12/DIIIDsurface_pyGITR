#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 09:17:22 2022

@author: audide
"""

#%%
# Reading the rustbca netcdf files

import netCDF4
from netCDF4 import Dataset

FileName='/Users/audide/Research/DIIIDsurface_pyGITR/examples/particles/Rustbca_Dataset_Energy.nc'
Rustbca_energy = Dataset(FileName, "r", format="NETCDF4")

FileName='/Users/audide/Research/DIIIDsurface_pyGITR/examples/particles/Rustbca_Dataset_Angle.nc'
Rustbca_angle = Dataset(FileName, "r", format="NETCDF4")

energy_rustbca = np.array(Rustbca_energy['energy'])

angle_rustbca = np.array(Rustbca_angle['angle'])


H_W_Rustbca  = np.array(Rustbca_energy['Sputtering_Hydrogen_Tungsten'])
H_Si_Rustbca = np.array(Rustbca_energy['Sputtering_Hydrogen_Silicon'])
H_C_Rustbca  = np.array(Rustbca_energy['Sputtering_Hydrogen_Carbon'])


H_W_angle_Rustbca  = np.array(Rustbca_angle['Sputtering_Hydrogen_Tungsten_Angle'])



#%%
# Plotting comparison plots

import matplotlib.pyplot as plt       



plt.figure()       
       
       

Energies = np.logspace(1, 3, 100)
Sputtering = Physical_Sputtering_Reflection_Plots.Sputtering_yields('H', 'C',Energies)


plt.plot(Energies,Sputtering.real,'ro',label='Surface Model')

plt.scatter(energy_rustbca,H_C_Rustbca,label='RustBCA')

# plt.xlim(0.1,1e5)
# plt.ylim(1e-5,2)

plt.xlabel('E_incident(eV)')
plt.ylabel('Sputering Yield (Y)')
plt.title(r'$H \rightarrow C$')
plt.legend()



