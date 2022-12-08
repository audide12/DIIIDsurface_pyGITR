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
ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
#B0 = -0.0002
nP=1000#10000


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

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/output/history.nc'
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

import netCDF4
from netCDF4 import Dataset

ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/background/new_profile.nc', 'w', format='NETCDF4')
r_dim = ncFile.createDimension('r_dimension', len(new_profile['r'])) # radius dimension
z_dim = ncFile.createDimension('z_dimension', len(new_profile['z'])) # z dimension


r_nc = ncFile.createVariable('r', np.float32, ('r_dimension',))
z_nc = ncFile.createVariable('z', np.float32, ('z_dimension',))


te_nc = ncFile.createVariable('te',np.float64,('r_dimension','z_dimension'))
ti_nc = ncFile.createVariable('ti',np.float64,('r_dimension','z_dimension'))
ne_nc = ncFile.createVariable('ne',np.float64,('r_dimension','z_dimension'))
ni_nc = ncFile.createVariable('ni',np.float64,('r_dimension','z_dimension'))

r_nc[:] = new_profile['r']
z_nc[:] = new_profile['z']
counter = 0

for z in new_profile['z']:
    te_nc[:,counter] = new_profile['te'][:]
    ti_nc[:,counter] = new_profile['te'][:]
    ne_nc[:,counter] = new_profile['ne'][:]
    ni_nc[:,counter] = new_profile['ne'][:]
    counter+=1
        

ncFile.close()

#%%

FileName='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/background/new_profile.nc'

new_new = Dataset(FileName, "r", format="NETCDF4")
#%%

import json
  
# Opening JSON file
f = open('background/176492_B_field_dimes_json.dat')
  
# returns JSON object as 
# a dictionary
data = json.load(f)
  
# Iterating through the json
# list
# for i in data['emp_details']:
#     print(i)
  
# Closing file
f.close()