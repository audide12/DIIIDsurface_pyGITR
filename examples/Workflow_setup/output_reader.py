#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:31:27 2022

@author: de
"""

# FileNameSurface='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/surface.nc'
# FileNameBoundary='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/boundary_values.nc'
# FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/history.nc'
# FileNamePositions='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/positions.nc'
# FileNameSpec='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/spec.nc'




# SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
# BoundaryData = Dataset(FileNameBoundary, "r", format="NETCDF4")
# HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
# PositionData = Dataset(FileNamePositions, "r", format="NETCDF4")
# SpecData = Dataset(FileNameSpec, "r", format="NETCDF4")


# ParticleInputData='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/particleConf.nc'

# ParticleInputData = Dataset(ParticleInputData, "r", format="NETCDF4")



Surface_time = [2,3]
C_C = np.ones((len(x1),len(Surface_time)))

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



# FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/input/surface_evolution_C.nc'

# SurfaceData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")


# C_C = SurfaceData['surface_concentration'][:,:]
# surface_time = SurfaceData['time'][:]
# surface_number = SurfaceData['surface_number'][:]
# counter = len(surface_time)

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit = np.array[PositionData['surfaceHit']]
count = 0
for i in surfacehit:
    if i == 2:
        count+=1
print(count)



# Edist = np.array(SurfaceData.variables['surfEDist'])
# gridE = np.array(SurfaceData.variables['gridE'])
# gridA = np.array(SurfaceData.variables['gridA'])


# E=np.sum(Edist[:,:,:],(0,2))
# A=np.sum(Edist[:,:,:],(0,1))
# S=np.sum(Edist[:,:,:],(1,2))  # 