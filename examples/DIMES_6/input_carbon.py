#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 12:38:45 2023

@author: de
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 17:38:55 2022

@author: de
"""

import os
import pyGITR
from pyGITR import Input
import numpy as np

#os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_4/input/particleConf.nc")

#os.system("mv /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_C.nc /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf.nc")

ParticleFile='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/particleConf_C.nc'
GeometryFile='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/gitrGeom.cfg'
B0 = 2.25
nP=1000000
dt=1e-8
nT=1e4



def make_input(nP,dt,nT,ParticleFile='particleConf_C.nc',GeometryFile='gitrGeom.cfg',folder='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/'):

    B0 = 2.25
    thetaB = -2
    phiB = 0
    # B0 = -2.25
    # thetaB = 2
    # phiB = 0

    # Initiallize input object
    i = Input()

    # Add structures to configuration file
    i.SetBField(B0=B0, theta = thetaB, phi = phiB)
    i.SetTimeStep(dt=dt,nT=nT)
    i.SetGeometryFile(GeometryFile)
    i.SetParticleSource(ParticleFile, nP=nP, Zmax=6, M=12, Z=0)
    i.SetSurfaces()
    i.SetDiagnostics()
    i.SetBackgroundPlasmaProfiles()
    i.SetSurfaceModel()
    i.SetGeomHash()
    i.SetGeomSheath()

    i.Input['flags']['USE_CUDA'] = 1
    
    # Set the standard flags
    i.Input['flags']['BIASED_SURFACE'] = 0
    i.Input['flags']['USE_SURFACE_POTENTIAL'] = 0
    i.Input['flags']['USETHERMALFORCE'] = 1
    i.Input['flags']['USEPERPDIFFUSION'] = 1
    i.Input['flags']['USESURFACEMODEL'] = 1 ##  changed from 1
    i.Input['flags']['USESHEATHEFIELD'] = 1 # 1 or 0, on or off
    i.Input['flags']['USEPRESHEATHEFIELD'] = 1  # 1 or 0, on or off
    i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
    i.Input['flags']['USEFRICTION'] = 1
    i.Input['flags']['USEANGLESCATTERING'] = 1
    i.Input['flags']['USEHEATING'] = 1
    i.Input['flags']['FORCE_EVAL'] = 0
    i.Input['flags']['USECYLSYMM'] = 1 # rotates the plasma with cylindrical geometry
    i.Input['flags']['USE3DTETGEOM'] = 1  # causes errors  for 3D simulations
    i.Input['flags']['SPECTROSCOPY'] = 2
    
    i.Input['flags']['PARTICLE_TRACKS'] = 1   #PARTICLE_TRACKS turns on/off even producing a history.nc file
    i.Input['diagnostics']['trackSubSampleFactor'] = 5e2
    
    # Set the INTERP flags
    i.Input['flags']['BFIELD_INTERP'] = 2
    i.Input['flags']['EFIELD_INTERP'] = 0
    i.Input['flags']['PRESHEATH_INTERP'] = 0 # 2 # 3
    i.Input['flags']['DENSITY_INTERP'] = 2
    i.Input['flags']['TEMP_INTERP'] = 2
    i.Input['flags']['FLOWV_INTERP'] = 2 #2 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
    i.Input['flags']['GRADT_INTERP'] = 2 #2 # 3, 1 R, 2 R+Z, 3 R+Z+Y

    # Set other options
    i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # in fields.cpp, only needs to be >0
    i.Input['impurityParticleSource']['Z'] = 6
    i.Input['impurityParticleSource']['source_material_Z'] = 6
    i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
    i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
    i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'
    i.Input['backgroundPlasmaProfiles']['Diffusion']['Dperp'] = 0.1
    i.Input['surfaceModel']['fileString'] = 'surface_model_C-C.nc'
    i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVr'] = 20000   # redundant
    
    i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVrString'] = 'vr'   # added
    i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVzString'] = 'vz'   # added
    i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVtString'] = 'vt'   # added
    i.Input['backgroundPlasmaProfiles']['gradT']['gradTeYString'] = 'gradTeY'   # added
    i.Input['backgroundPlasmaProfiles']['gradT']['gradTiYString'] = 'gradTiY'   # added
    

    i.Input['impurityParticleSource']['ionization']['fileString'] = 'ADAS_Rates_C.nc'
    i.Input['impurityParticleSource']['recombination']['fileString'] = 'ADAS_Rates_C.nc'



    i.Input['surfaces']['flux']['nE'] = 200
    i.Input['surfaces']['flux']['E0'] = 0
    i.Input['surfaces']['flux']['E'] = 200
    i.Input['surfaces']['flux']['nA'] = 90
    i.Input['surfaces']['flux']['A0'] = 0
    i.Input['surfaces']['flux']['A'] = 90

    # i.Input['diagnostics']['trackSubSampleFactor'] = 10
    # i.Input['diagnostics']['netx0'] = 1.38
    # i.Input['diagnostics']['netx1'] = 1.58
    # i.Input['diagnostics']['nX'] = 250
    # i.Input['diagnostics']['nety0'] = -0.1
    # i.Input['diagnostics']['nety1'] = 0.1
    # i.Input['diagnostics']['nY'] = 250
    # i.Input['diagnostics']['netz0'] = -1.25
    # i.Input['diagnostics']['netz1'] = -1.05
    # i.Input['diagnostics']['nZ'] = 250
    # i.Input['diagnostics']['densityChargeBins'] = 6

    # Write input file
    i.WriteInputFile(Folder=folder,OverWrite=True)

if __name__ == '__main__':
    make_input(nP,dt,nT)
#%%

import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/output_C_2/history.nc'
HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
x = np.array(HistoryData.variables['x'])
z = np.array(HistoryData.variables['z'])
y = np.array(HistoryData.variables['y'])

nT = HistoryData.dimensions['nT'].size
nP = HistoryData.dimensions['nP'].size
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(100,1000):
    ax.plot(x[i,:],y[i,:],z[i,:])
ax.set_zlabel('Z-Axis')
ax.set_xlabel('X-Axis')
ax.set_ylabel('Y-Axis')
    
#g.Plot_Geom(["DiMES"],fig=fig, ax=ax)
#g.Plot_Geom(["DiMES","BoundBox"],fig=fig, ax=ax)
#g.Plot(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)    

#np.argwhere(np.isnan(Gamma_C_ero_global))

