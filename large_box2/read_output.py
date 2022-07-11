#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:53:18 2021

@author: jguterl
"""
import electronvolt_num as units
import numpy as np
from matplotlib import pyplot as plt
plt.ion()
FileNameSurface='/home/jguterl/Dropbox/python/pyGITR/examples/large_box2/output/surface.nc'
FileNameParticle='/home/jguterl/Dropbox/python/pyGITR/examples/large_box2/output/particleSource.nc'
import netCDF4
from netCDF4 import Dataset
SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")
D={}
vx = np.array(ParticleData.variables['vx'])
vy = np.array(ParticleData.variables['vy'])
vz = np.array(ParticleData.variables['vz'])
amu = np.array(ParticleData.variables['amu'])
E0=1/2*amu*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2

Distrib = np.array(SurfaceData.variables['surfEDist'])
GridE = np.array(SurfaceData.variables['gridE'])
GridA = np.array(SurfaceData.variables['gridA'])
E=np.sum(Distrib[:,:,:],(0,2))
print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))
plt.figure()
plt.plot(GridE,E)
plt.hist(E0,GridE)





