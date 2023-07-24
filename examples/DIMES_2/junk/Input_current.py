#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 14:51:41 2022

@author: de
"""
# -*- coding: utf-8 -*-
#import electronvolt_num as units
import pyGITR
import numpy as np
#from single_microtrench  import Lxbc,Lybc
ParticleFile='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/particleConf.nc'
GeometryFile='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/gitrGeom.cfg'
B0 = 2.25
#B0 = -0.0002
nP=10000#10000


thetaB = 0#2
phiB = 0
AxisrotB = [0,1,0]

mi=2.0
Elem='C'

MassElem ={'C' : 12, 'D' : 2, 'W' : 184, 'Si' : 28}
Mimp = MassElem.get(Elem)
charge = 20 #20
Zmax = 6
#V= -60
V= -1

#Plasma parameter
Te = 25
Ti = Te
Timp = Te
vthi = 1e4#np.sqrt(eV*Ti/mp/mi)
vimp = 1e4#np.sqrt(eV*Ti/mp/Mimp)
cs = 1e4#np.sqrt(eV*2*Ti/mp/mi)



eV = 1.6e-19
mp = 1.67e-27

##

p = pyGITR.ParticleDistribution()
p.SetAttr('Np', nP+1)

# Set positions of particles
#p.SetAttr('x','Uniform') #set values of y and z with uniformly distributed values between -0.05 and 0.05
#p.SetAttr('y','Uniform')

p.SetAttr('x',1.5) 
p.SetAttr('y',0.0)
p.SetAttr('z',0.1)

# Set velocities of particles
p.SetAttr(['vz'],'Gaussian',sigma = 1.825e4,beta=3.16e14)
p.SetAttr(['vy','vx'],'Gaussian',sigma = 1.825e4,beta=3.16e14)

#p.Particles['vx'].mean()

# Rescale velocity by characteristic velocity
#vpara = vimp
#vperp = vimp
vpara = 10.0
vperp = 1.0

p.ScaleAttr(['vy','vx'],vperp)

p.ScaleAttr('vz',vpara)

#p.Rotate('v', AxisrotB, thetaB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB
vx = p.Particles['vx']
vy = p.Particles['vy']
vz = p.Particles['vz']
E= 1/2*Mimp*(vx**2+vy**2+vz**2)*mp/eV

# Write particle distribution in netcdf file
p.WriteParticleFile(ParticleFile)


#%%
import os
from pyGITR import Input

ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0
NP = 5000


# Initiallize input object
i = Input()

# Add structures to configuration file
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-11,nT=5e4)
i.SetGeometryFile(GeometryFile)
#i.SetParticleSource(ParticleFile, nP=NP, Zmax=74, M=183, Z=4)

i.SetParticleSource(ParticleFile, nP=NP, Zmax=6, M=12, Z=1)#Z 6 ask about this
i.SetSurfaces()
i.SetDiagnostics()
i.SetBackgroundPlasmaProfiles()
i.SetSurfaceModel()

# i.SetGeomHash()
# i.Input['flags']['GEOM_HASH'] = 1
# i.SetGeomSheath()
# i.Input['flags']['GEOM_HASH_SHEATH'] = 1

"""
# Set the standard flags
Hi there
i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # in fields.cpp, only needs to be >0

i.Input['flags']['BFIELD_INTERP'] = 2
i.Input['flags']['EFIELD_INTERP'] = 0
i.Input['flags']['PRESHEATH_INTERP'] = 2 # 3
i.Input['flags']['DENSITY_INTERP'] = 2
i.Input['flags']['TEMP_INTERP'] = 2
i.Input['flags']['FLOWV_INTERP'] = 2 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
i.Input['flags']['GRADT_INTERP'] = 2 # 3, 1 R, 2 R+Z, 3 R+Z+Y
"""


i.Input['flags']['BIASED_SURFACE'] = 0
i.Input['flags']['USETHERMALFORCE'] = 1
i.Input['flags']['USEPERPDIFFUSION'] = 1
i.Input['flags']['USE_SURFACE_POTENTIAL'] = 0
i.Input['flags']['FORCE_EVAL'] = 0

i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # 1 in fields.cpp, only needs to be >0
i.Input['flags']['BFIELD_INTERP'] = 2   #2

i.Input['flags']['GENERATE_LC'] = 0
i.Input['flags']['LC_INTERP'] = 0 # 0, 0 no Lc, 2 get from file, 3 consider flowV
i.Input['flags']['EFIELD_INTERP'] = 0
i.Input['flags']['PRESHEATH_INTERP'] = 0 # 3     changed
i.Input['flags']['DENSITY_INTERP'] = 2 # 0,3
i.Input['flags']['TEMP_INTERP'] = 2 # 0,3
i.Input['flags']['FLOWV_INTERP'] = 0 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D   changed
i.Input['flags']['GRADT_INTERP'] = 0 # 3, 1 R, 2 R+Z, 3 R+Z+Y         changed
i.Input['flags']['USECYLSYMM'] = 1 # 1 rotates the plasma with cylindrical geometry
#i.Input['flags']['USE3DTETGEOM'] = 0  # causing errors
i.Input['flags']['SPECTROSCOPY'] = 2

# i.SetConnectionLength()
# i.SetForceEvaluation()

# Set the standard options
i.Input['flags']['USESURFACEMODEL'] = 1
i.Input['surfaceModel']['fileString'] = 'surface_model_C-C.nc'

i.Input['impurityParticleSource']['ionization']['fileString'] = 'ADAS_Rates_C.nc'
i.Input['impurityParticleSource']['recombination']['fileString'] = 'ADAS_Rates_C.nc'


i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
i.Input['flags']['USEFRICTION'] = 1 # 1
i.Input['flags']['USEANGLESCATTERING'] = 1
i.Input['flags']['USESHEATHEFIELD'] = 1
i.Input['flags']['USEPRESHEATHEFIELD'] = 1  # 1 or 0, on or off
i.Input['flags']['USEHEATING'] = 1 # 1


i.Input['impurityParticleSource']['C'] = 6   # change this 
i.Input['impurityParticleSource']['source_material_C'] = 6  # change this 
i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000

i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'
i.Input['backgroundPlasmaProfiles']['Diffusion']['Dperp'] = 0.1

'''
i.Input['impurityParticleSource']['Z'] = 6
i.Input['impurityParticleSource']['source_material_Z'] = 6
i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'
i.Input['backgroundPlasmaProfiles']['Diffusion']['Dperp'] = 0.1
i.Input['surfaceModel']['fileString'] = 'surface_model_C-C.nc'
i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVr'] = 20000

i.Input['impurityParticleSource']['ionization']['fileString'] = 'ADAS_Rates_C.nc'
i.Input['impurityParticleSource']['recombination']['fileString'] = 'ADAS_Rates_C.nc'

'''


#why these

i.Input['surfaces']['flux']['nE'] = 200
i.Input['surfaces']['flux']['E0'] = 0
i.Input['surfaces']['flux']['E'] = 200
i.Input['surfaces']['flux']['nA'] = 90
i.Input['surfaces']['flux']['A0'] = 0
i.Input['surfaces']['flux']['A'] = 90


i.Input['diagnostics']['trackSubSampleFactor'] = 10
i.Input['diagnostics']['netx0'] = 1.38
i.Input['diagnostics']['netx1'] = 1.58
i.Input['diagnostics']['nX'] = 250
i.Input['diagnostics']['nety0'] = -0.1
i.Input['diagnostics']['nety1'] = 0.1
i.Input['diagnostics']['nY'] = 250
i.Input['diagnostics']['netz0'] = -1.25
i.Input['diagnostics']['netz1'] = -1.05
i.Input['diagnostics']['nZ'] = 250
i.Input['diagnostics']['densityChargeBins'] = 6


# Write input file

i.WriteInputFile(Folder='input',OverWrite=True)


#%%

import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/output/history.nc'
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
    
g.Plot_Geom(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)        
#g.Plot(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)    




#%% that works

import os
from pyGITR import Input

ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0
NP = 5000

# Initiallize input object
i = Input()

# Add structures to configuration file
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-11,nT=5e4)
i.SetGeometryFile(GeometryFile)
#i.SetParticleSource(ParticleFile, nP=NP, Zmax=74, M=183, Z=4)

i.SetParticleSource(ParticleFile, nP=NP, Zmax=6, M=12, Z=2)#Z 6
i.SetSurfaces()
i.SetDiagnostics()
i.SetBackgroundPlasmaProfiles()
i.SetSurfaceModel()

# i.SetGeomHash()
# i.Input['flags']['GEOM_HASH'] = 1
# i.SetGeomSheath()
# i.Input['flags']['GEOM_HASH_SHEATH'] = 1

i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 0 # 1 in fields.cpp, only needs to be >0
i.Input['flags']['BFIELD_INTERP'] = 0   #2

i.Input['flags']['GENERATE_LC'] = 0
i.Input['flags']['LC_INTERP'] = 0 # 0, 0 no Lc, 2 get from file, 3 consider flowV
i.Input['flags']['EFIELD_INTERP'] = 0
i.Input['flags']['PRESHEATH_INTERP'] = 0 # 3
i.Input['flags']['DENSITY_INTERP'] = 0 # 0,3
i.Input['flags']['TEMP_INTERP'] = 0 # 0,3
i.Input['flags']['FLOWV_INTERP'] = 0 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
i.Input['flags']['GRADT_INTERP'] = 0 # 3, 1 R, 2 R+Z, 3 R+Z+Y
i.Input['flags']['USECYLSYMM'] = 0 # 1 rotates the plasma with cylindrical geometry

# i.SetConnectionLength()
# i.SetForceEvaluation()

# Set the standard options
i.Input['flags']['USESURFACEMODEL'] = 0
i.Input['flags']['USECOULOMBCOLLISIONS'] = 0
i.Input['flags']['USEFRICTION'] = 0 # 1
i.Input['flags']['USEANGLESCATTERING'] = 0
i.Input['flags']['USESHEATHEFIELD'] = 0
i.Input['flags']['USEHEATING'] = 0 # 1
i.Input['impurityParticleSource']['C'] = 6   # change this 
i.Input['impurityParticleSource']['source_material_C'] = 6  # change this 
i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000

i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'

# Write input file
i.WriteInputFile(Folder='input',OverWrite=True)


#%%

import netCDF4
from netCDF4 import Dataset


FileName='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_model_C-C.nc'
SurfaceData = Dataset(FileName, "r", format="NETCDF4")


