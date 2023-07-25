# import electronvolt_num as units
# import numpy as np
# from matplotlib import pyplot as plt
# from netCDF4 import Dataset
# plt.ion()
# FileNameSurface='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/surface.nc'
# FileNameParticle='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/particleSource.nc'
# import netCDF4
# import os
# from SimManager import rget


# FileNameSurface='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/surface.nc'
# FileNameBoundary='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/boundary_values.nc'
# FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/history.nc'
# FileNamePositions='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/positions.nc'
# FileNameSpec='/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/output/spec.nc'


# ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")


# SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
# BoundaryData = Dataset(FileNameBoundary, "r", format="NETCDF4")
# HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
# PositionData = Dataset(FileNamePositions, "r", format="NETCDF4")
# SpecData = Dataset(FileNameSpec, "r", format="NETCDF4")


# def GetDistribSurface(SurfaceData):
#     Data = {}
#     Data['Distrib'] = np.array(SurfaceData.variables['surfEDist'])
#     Data['GridE'] = np.array(SurfaceData.variables['gridE'])
#     Data['GridA'] = np.array(SurfaceData.variables['gridA'])
#     Data['E'] = np.sum(Data['Distrib'][:,:,:],(0,2))
#     Data['A'] = np.sum(Data['Distrib'][:,:,:],(0,1))
#     Data['S'] = np.sum(Data['Distrib'][:,:,:],(1,2))
#     Data['Tot'] = np.sum(Data['Distrib'][:,:,:],(0,1,2))
#     return Data

# def GetDistribParticle(ParticleData):
#     Data={}
#     for k in ParticleData.variables.keys():
#         Data[k] = np.array(ParticleData.variables.get(k))
#     v = np.sqrt(Data['vz']**2+Data['vy']**2+Data['vx']**2)
#     vp = np.sqrt(Data['vy']**2+Data['vx']**2)
#     #Data['E'] = 1/2*Data['amu']*(v**2)*units.mp/units.eV*(units.V)**2
#     Data['E'] = 1/2*(v**2)*units.mp/units.eV*(units.V)**2
#     Data['theta'] = np.arccos(np.abs(Data['vz'])/v)
#     Data['phi'] = np.arccos(np.abs(Data['vx'])/vp)
#     return Data

# def GetHistoryParticle(ParticleData):
#     Data={}
#     for k in ParticleData.variables.keys():
#         Data[k] = np.array(ParticleData.variables.get(k))
#     return Data




# class DistributionAnalysis():
    
#     SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
#     ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")
#     D={}
#     vx = np.array(ParticleData.variables['vx'])
#     vy = np.array(ParticleData.variables['vy'])
#     vz = np.array(ParticleData.variables['vz'])
#     #amu = np.array(ParticleData.variables['amu'])
#     #E0=1/2*amu*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2
#     E0=1/2*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2
    
#     Distrib = np.array(SurfaceData.variables['surfEDist'])
#     GridE = np.array(SurfaceData.variables['gridE'])
#     GridA = np.array(SurfaceData.variables['gridA'])
#     E=np.sum(Distrib[:,:,:],(0,2))
#     A=np.sum(Distrib[:,:,:],(0,1))
#     S=np.sum(Distrib[:,:,:],(1,2))
#     print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))
#     fig = plt.figure()
#     ax = fig.subplots(1,4)
#     ax[0].plot(GridE,E)
#     ax[0].hist(E0,GridE)
#     ax[1].plot(GridA,A)
#     ax[2].hist(vx)
#     ax[2].hist(vy,alpha=0.2)
#     ax[2].hist(vz,alpha=0.2)
#     ax[3].plot(np.arange(0,S.shape[0]),S)
    