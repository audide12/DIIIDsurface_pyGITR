# -*- coding: utf-8 -*-

from pyGITR import GITRInput
from pyGITR import ParticleDistribution

ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0


##
p = ParticleDistribution()
p.SetAttr('Np', 10000)
p.SetAttr(['z','y'],'Uniform',xmin=-0.05,xmax=0.05) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p.SetAttr('x',-0.01) # set all values of x to -0.01
p.SetAttr(['vx'],'Gaussian',beta=1.6)
p.SetAttr(['vy','vz'],'Gaussian')
vpara = 1e4
vperp = 1e5
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)
p.RotateAngle('v', theta = thetaB, phi = phiB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB

# Write particle distribution in netcdf file
p.WriteParticleFile(ParticleFile)

i = GITRInput()
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-8)
i.SetGeometryFile(GeometryFile)
i.SetParticleSource(ParticleFile, Zmax=74, M=183, Z=4)
i.WriteInputFile()





