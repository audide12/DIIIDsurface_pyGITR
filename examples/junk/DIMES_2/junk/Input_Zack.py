#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 17:13:10 2022

@author: de
"""

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
i.SetTimeStep(dt=1e-8,nT=1e3)
i.SetGeometryFile(GeometryFile)
i.SetParticleSource(ParticleFile, nP=NP, Zmax=74, M=183, Z=4)
i.SetSurfaces()
i.SetDiagnostics()
i.SetBackgroundPlasmaProfiles()
i.SetSurfaceModel()

# i.SetGeomHash()
# i.Input['flags']['GEOM_HASH'] = 1
# i.SetGeomSheath()
# i.Input['flags']['GEOM_HASH_SHEATH'] = 1

i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # in fields.cpp, only needs to be >0
i.Input['flags']['BFIELD_INTERP'] = 2

i.Input['flags']['GENERATE_LC'] = 0
i.Input['flags']['LC_INTERP'] = 0 # 0, 0 no Lc, 2 get from file, 3 consider flowV
i.Input['flags']['EFIELD_INTERP'] = 0
i.Input['flags']['PRESHEATH_INTERP'] = 0 # 3
i.Input['flags']['DENSITY_INTERP'] = 0 # 3
i.Input['flags']['TEMP_INTERP'] = 0 # 3
i.Input['flags']['FLOWV_INTERP'] = 0 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
i.Input['flags']['GRADT_INTERP'] = 0 # 3, 1 R, 2 R+Z, 3 R+Z+Y
i.Input['flags']['USECYLSYMM'] = 1 # rotates the plasma with cylindrical geometry

# i.SetConnectionLength()
# i.SetForceEvaluation()

# Set the standard options
i.Input['flags']['USESURFACEMODEL'] = 0
i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
i.Input['flags']['USEFRICTION'] = 1
i.Input['flags']['USEANGLESCATTERING'] = 1
i.Input['flags']['USEHEATING'] = 1
i.Input['impurityParticleSource']['Z'] = 74
i.Input['impurityParticleSource']['source_material_Z'] = 74
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
#Reading geometry files

filename="input/gitrGeom.cfg"

with io.open(filename) as f:
    config = libconf.load(f)

x1 = np.array(config.geom.x1)
x2 = np.array(config.geom.x2)
x3 = np.array(config.geom.x3)
y1 = np.array(config.geom.y1)
y2 = np.array(config.geom.y2)
y3 = np.array(config.geom.y3)
z1 = np.array(config.geom.z1)
z2 = np.array(config.geom.z2)
z3 = np.array(config.geom.z3)
area = np.array(config.geom.area)
surf = np.array(config.geom.surface)
Z = np.array(config.geom.Z)
materialSurfaceInidces = np.nonzero(Z)
surfIndArray = np.asarray(materialSurfaceInidces)
a = np.array(config.geom.a)
b = np.array(config.geom.b)
c = np.array(config.geom.c)
d = np.array(config.geom.d)

in_direction = np.array(config.geom.inDir)
plane_norm = np.array(config.geom.plane_norm)
