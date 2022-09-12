#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 17:28:07 2022

@author: de
"""

from netCDF4 import Dataset
import numpy as np
import electronvolt_num as units
from matplotlib import pyplot as plt

#%%

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/history.nc'
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
    

#%%    
FileNameSurface='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/surface.nc'
FileNameParticle='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/particleSource.nc'    


SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")
D={}
vx = np.array(ParticleData.variables['vx'])
vy = np.array(ParticleData.variables['vy'])
vz = np.array(ParticleData.variables['vz'])
#amu = np.array(ParticleData.variables['amu'])
#E0=1/2*amu*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2
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

#%%
# History plotting

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/history.nc'
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
g.Plot(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)    

#%% 
# Setting up the initial C_C in surface_evolution_C.nc

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
import matplotlib.pyplot as plt    

surface_index_C = 0
surface_index_W = 0
plt.figure()
plt.plot(Surface_time,C_C[surface_index_C,:],'k',label='C_C')
plt.plot(Surface_time,C_W[surface_index_W,:],'r',label='C_W')
plt.legend()
plt.show()


#%%

FileName='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/particleConf_C.nc'
FileName2='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/particleSource_test.nc'

ParticleData = Dataset(FileName, "r", format="NETCDF4")
ParticleData2 = Dataset(FileName2, "r", format="NETCDF4")




    