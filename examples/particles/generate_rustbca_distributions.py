#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:05:14 2022

@author: audide
"""

from libRustBCA.pybca import *
from scripts.materials import *
import numpy as np
import time
import netCDF4
from netCDF4 import Dataset


tungsten['Eb'] = 3.0

num_samples = 100000

energies = np.logspace(1, 3, 100)


Sputtering_Hydrogen_Tungsten = [sputtering_yield(hydrogen, tungsten, energy, 0.0, num_samples) for energy in energies] # for normal incidence ?

Sputtering_Hydrogen_Carbon = [sputtering_yield(hydrogen, carbon, energy, 0.0, num_samples) for energy in energies] # for normal incidence ?

Sputtering_Hydrogen_Silicon = [sputtering_yield(hydrogen, silicon, energy, 0.0, num_samples) for energy in energies] # for normal incidence ? 


energy_constant = 100
Angle_space = np.linspace(0,90,num=100) # deg

Sputtering_Hydrogen_Tungsten_Angle = [sputtering_yield(hydrogen, tungsten, energy_constant, angle, num_samples) for angle in Angle_space] # for oblique incidence


# Writing in a netcdf file for energy

ncFile = netCDF4.Dataset('Rustbca_Dataset_Energy.nc', 'w', format='NETCDF4')

energy_dim = ncFile.createDimension('energy_dim', len(energies)) # energy dimension

energy_variable = ncFile.createVariable('energy', np.float32, ('energy_dim',))

Sputtering_H_W_variable = ncFile.createVariable('Sputtering_Hydrogen_Tungsten',np.float64,('energy_dim'))

Sputtering_H_C_variable = ncFile.createVariable('Sputtering_Hydrogen_Carbon',np.float64,('energy_dim'))

Sputtering_H_Si_variable = ncFile.createVariable('Sputtering_Hydrogen_Silicon',np.float64,('energy_dim'))

energy_variable[:] = energies
Sputtering_H_W_variable[:] = Sputtering_Hydrogen_Tungsten
Sputtering_H_C_variable[:] = Sputtering_Hydrogen_Carbon
Sputtering_H_Si_variable[:] = Sputtering_Hydrogen_Silicon

ncFile.close()

# Writing in a netcdf file for energy

ncFile = netCDF4.Dataset('Rustbca_Dataset_Angle.nc', 'w', format='NETCDF4')

angle_dim = ncFile.createDimension('angle_dim', len(Angle_space)) # energy dimension

angle_variable = ncFile.createVariable('angle', np.float32, ('angle_dim',))

Sputtering_H_W_variable = ncFile.createVariable('Sputtering_Hydrogen_Tungsten_Angle',np.float64,('angle_dim'))

angle_variable[:] = Angle_space
Sputtering_H_W_variable[:] = Sputtering_Hydrogen_Tungsten_Angle

ncFile.close()
