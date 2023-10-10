# -*- coding: utf-8 -*-

import os
from pyGITR import Input

def make_input(nP,dt,nT,ParticleFile='particleConf.nc',GeometryFile='gitrGeom.cfg',folder='../input/'):

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
    # i.SetGeomHash()
    # i.SetGeomSheath()

    # Set the standard flags
    i.Input['flags']['BIASED_SURFACE'] = 0
    i.Input['flags']['USE_SURFACE_POTENTIAL'] = 0
    i.Input['flags']['USETHERMALFORCE'] = 1
    i.Input['flags']['USEPERPDIFFUSION'] = 1
    i.Input['flags']['USESURFACEMODEL'] = 1
    i.Input['flags']['USESHEATHEFIELD'] = 1 # 1 or 0, on or off
    i.Input['flags']['USEPRESHEATHEFIELD'] = 1  # 1 or 0, on or off
    i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
    i.Input['flags']['USEFRICTION'] = 1
    i.Input['flags']['USEANGLESCATTERING'] = 1
    i.Input['flags']['USEHEATING'] = 1
    i.Input['flags']['FORCE_EVAL'] = 0
    i.Input['flags']['USECYLSYMM'] = 1 # rotates the plasma with cylindrical geometry
    i.Input['flags']['USE3DTETGEOM'] = 0
    i.Input['flags']['SPECTROSCOPY'] = 2

    # Set the INTERP flags
    i.Input['flags']['BFIELD_INTERP'] = 2
    i.Input['flags']['EFIELD_INTERP'] = 0
    i.Input['flags']['PRESHEATH_INTERP'] = 2 # 3
    i.Input['flags']['DENSITY_INTERP'] = 2
    i.Input['flags']['TEMP_INTERP'] = 2
    i.Input['flags']['FLOWV_INTERP'] = 2 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
    i.Input['flags']['GRADT_INTERP'] = 2 # 3, 1 R, 2 R+Z, 3 R+Z+Y

    # Set other options
    i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # in fields.cpp, only needs to be >0
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
    i.WriteInputFile(Folder=folder,OverWrite=True)

if __name__ == '__main__':
    make_input(1000)
