
import os
from pyGITR import Input


os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf.nc")

os.system("mv /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_Si.nc /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf.nc")


ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0
NP = 10000

# Initiallize input object
i = Input()

# Add structures to configuration file
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-11,nT=5e4)
i.SetGeometryFile(GeometryFile)
#i.SetParticleSource(ParticleFile, nP=NP, Zmax=74, M=183, Z=4)

i.SetParticleSource(ParticleFile, nP=NP, Zmax=14, M=28, Z=2)#Z 6
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
i.Input['impurityParticleSource']['Si'] = 14   # change this 
i.Input['impurityParticleSource']['source_material_Si'] = 14  # change this 
i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000

i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'

# Write input file
i.WriteInputFile(Folder='input',OverWrite=True)





