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
FileNameSurface='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/surface.nc'
FileNameParticle='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/particleSource.nc'
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
#%%
FileNameHistory='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/history.nc'
HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
x = np.array(HistoryData.variables['x'])
z = np.array(HistoryData.variables['z'])
y = np.array(HistoryData.variables['y'])

nT = HistoryData.dimensions['nT'].size
nP = HistoryData.dimensions['nP'].size
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(0,nP):
    ax.plot(x[i,:],y[i,:],z[i,:])





