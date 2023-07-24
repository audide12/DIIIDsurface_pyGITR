#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 16:55:17 2022

@author: de
"""


import pyGITR
import numpy as np
from netCDF4 import Dataset

FileNameSurface='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/surface.nc'

SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")

E0=1/2*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2

Distrib = np.array(SurfaceData.variables['surfEDist'])
GridE = np.array(SurfaceData.variables['gridE'])
GridA = np.array(SurfaceData.variables['gridA'])
E=np.sum(Distrib[:,:,:],(0,2))
A=np.sum(Distrib[:,:,:],(0,1))
S=np.sum(Distrib[:,:,:],(1,2))
print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))
fig = plt.figure()
ax = fig.subplots(1,4)
ax[0].plot(GridE,E)
ax[0].hist(E0,GridE)
ax[1].plot(GridA,A)
ax[2].hist(vx)
ax[2].hist(vy,alpha=0.2)
ax[2].hist(vz,alpha=0.2)
ax[3].plot(np.arange(0,S.shape[0]),S)


#SurfaceData.variables
#SurfaceData.variables.keys()
#SurfaceData.variables['surfEDist']
#Distrib


