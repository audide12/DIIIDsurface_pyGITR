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
p.SetAttr('Np', 1000)
p.SetAttr(['z','y'],'Uniform',xmin=-0.005,xmax=0.005) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p.SetAttr('x',-0) # set all values of x to -0.01
p.SetAttr(['vx'],'Gaussian',beta=1.6)
p.SetAttr(['vy','vz'],'Gaussian')
vpara = 1e4
vperp = 1e5
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)
#p.RotateAngle('v', theta = thetaB, phi = phiB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB

# Write particle distribution in netcdf file
p.WriteParticleFile(ParticleFile)

i = GITRInput()
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-8)
i.SetGeometryFile(GeometryFile)
i.SetParticleSource(ParticleFile, Zmax=74, M=183, Z=4)

i.Input['flags'] = {
                    'USE_CUDA':0,
                    'USE_MPI':1,
                    'USE_OPENMP':1,
                    'USE_IONIZATION':1,
                    'USE_RECOMBINATION':1,
                    'USEPERPDIFFUSION':0,
                    'USEPARDIFFUSION':0,
                    'USECOULOMBCOLLISIONS':0,
                    'USEFRICTION':0,
                    'USEANGLESCATTERING':0,
                    'USEHEATING':0,
                    'USETHERMALFORCE':0,
                    'USESURFACEMODEL':0,
                    'USE_SURFACE_POTENTIAL':0,
                    'USE_SHEATHEFIELD':1,
                    'BIASED_SURFACE':0,
                    'USEPRESHEATHEFIELD':0,
                    'BFIELD_INTERP':0,
                    'LC_INTERP':0,
                    'GENERATE_LC':0,
                    'EFIELD_INTERP':0,
                    'PRESHEATH_INTERP':0,
                    'DENSITY_INTERP':0,
                    'TEMP_INTERP':0,
                    'FLOWV_INTERP':0,
                    'GRADT_INTERP':0,
                    'ODEINT':0,
                    'FIXED_SEEDS':1,
                    'PARTICLESEEDS':1,
                    'GEOM_TRACE':0,
                    'GEOM_HASH':0,
                    'GEOM_HASH_SHEATH':0,
                    'PARTICLE_TRACKS':1,
                    'PARTICLE_SOURCE_SPACE':0,
                    'PARTICLE_SOURCE_ENERGY':0,
                    'PARTICLE_SOURCE_ANGLE':0,
                    'PARTICLE_SOURCE_FILE':1,
                    'SPECTROSCOPY':3,
                    'USE3DTETGEOM':1,
                    'USECYLSYMM':0,
                    'USEFIELDALIGNEDVALUES':0,
                    'FLUX_EA':1,
                    'FORCE_EVAL':0,
                    'USE_SORT':0,
                    'CHECK_COMPATIBILITY':1
                    }

i.WriteInputFile(Folder='input', OverWrite=True)





