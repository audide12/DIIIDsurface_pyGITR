#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:50:33 2022

@author: audide
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generation of particles distribution for GITR.
@author: guterl
"""
import pyGITR
from pyGITR.Particles import *
ParticleFile='particleConf2.nc'
p = ParticleDistribution()

# Attributes of particles x,y,z,v,x,vy,vz are stored in the dictionary p.Particles
# Display the list of attributes of particles.


# First, set numbers of particles. This is not affecting exisiting attributes
# but only the generation of those attributes
p.SetAttr('Np', 2)

p.ShowAttr()

# Set distribution for attributes of particles
# First, show list of available distribution pdf
p.ShowAvailablePdfs()

# Set positions of particles
p.SetAttr(['z','y'],'Uniform',xmin=-0.05,xmax=0.05) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p.SetAttr('x',-0.01) # set all values of x to -0.01

# Set velocities of particles
p.SetAttr(['vx'],'Gaussian',beta=1.6)
p.SetAttr(['vy','vz'],'Gaussian')

# Rescale velocity by characteristic velocity
vpara = 1e4
vperp = 1e5
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)

# Rotate velocity of particles to put in the reference frame of the magnetic field

# Define rotation
thetaB = 2 #degree by default, set Degree=False in RotateVector for radian
AxisrotB = [0,1,0]

# Check results of the rotation
B0 = [1,0,0]
B = RotateVector(B0, AxisrotB, thetaB)
PlotVector(B)

# Rotate velocity vector
p.Rotate('v', AxisrotB, thetaB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB

# Write particle distribution in netcdf file
#p.WriteParticleFile('particleConf.nc')

# Distribution can be also defined with a user-defined pdf. Must be formatted as
# f(x, **kwargs)
def LevyDistrib(x, c=1, mu=0):
    return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)
p.SetAttr(['vx'],LevyDistrib, x=np.linspace(0.001,10,1000), c=2, mu=0)
p.PlotAttr('vx')
#p.WriteParticleFile('particleConf2.nc')