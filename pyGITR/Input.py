#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 19:57:04 2021
@author: guterl
"""
import libconf, click
import os, io
from pyGITR.math_helper import *

class Input():
    def __init__(self):
        self.Input = {'backgroundPlasmaProfiles':{'BField':{}}}
        self.Input['flags'] = {
                    'USE_CUDA':0,
                    'USE_MPI':0,
                    'USE_OPENMP':0,
                    'USE_IONIZATION':1,
                    'USE_RECOMBINATION':1,
                    'USEPERPDIFFUSION':0,
                    # 'USEPARDIFFUSION':0,
                    'USECOULOMBCOLLISIONS':0,
                    'USEFRICTION':0,
                    'USEANGLESCATTERING':0,
                    'USEHEATING':0,
                    'USETHERMALFORCE':0,
                    'USESURFACEMODEL':0,
                    'USE_SURFACE_POTENTIAL':0,
                    'USESHEATHEFIELD':1,
                    'BIASED_SURFACE':0,
                    'USEPRESHEATHEFIELD':0,
                    'LC_INTERP':0,
                    'GENERATE_LC':0,
                    'BFIELD_INTERP':0,
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
                    'USE_ADAPTIVE_DT':0,
                    'CHECK_COMPATIBILITY':0
                    }

    def WriteInputFile(self, FileName='gitrInput.cfg', Folder='', OverWrite=False):
        FileName = os.path.abspath(os.path.join(Folder,FileName))
        print('Writing input config into {} ...'.format(FileName))

        if os.path.exists(FileName):
            if not OverWrite and not click.confirm('File {} exists.\n Do you want to overwrite it?'.format(FileName), default=True):
                return

        with io.open(FileName,'w') as f:
            libconf.dump(self.Input,f)

    def SetBField(self,Br=None, By=None, Bz=None, B0=None, theta=None, phi=None, Degree=True):
        if self.Input.get('BField') is None:
            self.Input['BField'] = {}

        if Br is not None and By is not None and Bz is not None:
            self.Input['BField']['r'] = Br
            self.Input['BField']['z'] = Bz
            self.Input['BField']['y'] = By
            self.Input['backgroundPlasmaProfiles']['BField']['r'] = Br
            self.Input['backgroundPlasmaProfiles']['BField']['z'] = Bz
            self.Input['backgroundPlasmaProfiles']['BField']['y'] = By

        elif B0 is not None and theta is not None and phi is not None:
            theta = 2 #degree by default, set Degree=False in RotateVector for radian
            Axisrotx = [1,0,0]
            Axisroty = [0,1,0]
            Axisrotz = [0,0,1]
            # Btot = [B0,0,0]
            Btot = [0,B0,0]
            # B1 = RotateVector(Btot, Axisroty, theta, Degree)
            # B  = RotateVector(B1, Axisrotz, phi, Degree)
            B  = RotateVector(Btot, Axisrotx, theta, Degree)
            self.Input['BField']['r'] = B[0]
            self.Input['BField']['z'] = B[2]
            self.Input['BField']['y'] = B[1]
            self.Input['backgroundPlasmaProfiles']['BField']['r'] = B[0]
            self.Input['backgroundPlasmaProfiles']['BField']['z'] = B[2]
            self.Input['backgroundPlasmaProfiles']['BField']['y'] = B[1]
        else:
            raise KeyError('Br, By and Bt not None or B, theta and phi not None')

        print('Bz={}; Br={}; Bt={}'.format(self.Input['BField']['z'],  self.Input['BField']['r'] ,self.Input['BField']['y'] ))

        self.Input['BField']['rString']= 'r'
        self.Input['BField']['zString']= 'z'
        self.Input['BField']['yString']= 'y'

    def SetSurfaceModel(self):
        self.Input['surfaceModel'] = {'fileString' : "ftridynSelf.nc",
            'nEsputtRefCoeffString' : "nE",
            'nAsputtRefCoeffString' : "nA",
            'nEsputtRefDistInString' : "nE",
            'nAsputtRefDistInString' : "nA",
            'nEsputtRefDistOutString' : "nEdistBins",
            'nEsputtRefDistOutStringRef' : "nEdistBinsRef",
            'nAsputtRefDistOutString' : "nAdistBins",
            'E_sputtRefCoeff' : "E",
            'A_sputtRefCoeff' : "A",
            'E_sputtRefDistIn' : "E",
            'A_sputtRefDistIn' : "A",
            'E_sputtRefDistOut' : "eDistEgrid",
            'E_sputtRefDistOutRef' : "eDistEgridRef",
            'Aphi_sputtRefDistOut' : "phiGrid",
            'Atheta_sputtRefDistOut' : "thetaGrid",
            'sputtYldString' : "spyld",
            'reflYldString' : "rfyld",
            'EDist_Y' : "energyDist",
            'AphiDist_Y' : "cosXDist",
            'AthetaDist_Y' : "cosYDist",
            'EDist_R' : "energyDistRef",
            'AphiDist_R' : "cosXDistRef",
            'AthetaDist_R' : "cosYDistRef"
            }

    def SetConnectionLength(self):
        self.Input['connectionLength'] = {'nTraceSteps' : 5000,
            'dr' : 0.0001,
            'netx0' : -0.076,
            'netx1' : 0.076,
            'nX' : 80,
            'nety0' : -0.076,
            'nety1' : 0.076,
            'nY' : 70,
            'netz0' : -0.05,
            'netz1' : 0.2,
            'nZ' : 100,
            'fileString' : "LcS.nc",
            'gridNrString' : "nR",
            'gridNyString' : "nY",
            'gridNzString' : "nZ",
            'gridRString' : "gridR",
            'gridYString' : "gridY",
            'gridZString' : "gridZ",
            'LcString' : "Lc",
            'SString' : "s",
            'noIntersectionString' : "noIntersection",
            }

    def SetGeomHash(self):
        self.Input['geometry_hash'] = {'nHashes' : 1,
            'hashX0' : -0.0625,
            'hashX1' : 0.0625,
            'hashY0' : -0.0625,
            'hashY1' : 0.0625,
            'hashZ0' : -0.0625,
            'hashZ1' : 0.0625,
            'nR_closeGeom' : 80,
            'nY_closeGeom' : 80,
            'nZ_closeGeom' : 80,
            'n_closeGeomElements' : 20,
            'fileString' : "geomHash0.nc",
            'gridNrString' : "nR",
            'gridNyString' : "nY",
            'gridNzString' : "nZ",
            'nearestNelementsString' : "n",
            'gridRString' : "gridR",
            'gridYString' : "gridY",
            'gridZString' : "gridZ",
            'closeGeomString' : "hash"
            }

    def SetGeomSheath(self):
        self.Input['geometry_sheath'] = {'nHashes' : 1,
            'hashX0' : -0.0625,
            'hashX1' : 0.0625,
            'hashY0' : -0.0625,
            'hashY1' : 0.0625,
            'hashZ0' : -0.0625,
            'hashZ1' : 0.0625,
            'nR_closeGeom' : 80,
            'nY_closeGeom' : 80,
            'nZ_closeGeom' : 80,
            'n_closeGeomElements' : 20,
            'fileString' : "geomHash_sheath.nc",
            'gridNrString' : "nR",
            'gridNyString' : "nY",
            'gridNzString' : "nZ",
            'nearestNelementsString' : "n",
            'gridRString' : "gridR",
            'gridYString' : "gridY",
            'gridZString' : "gridZ",
            'closeGeomString' : "hash"
            }

    def SetGeometryFile(self, FileName):
        if self.Input.get('geometry') is None:
            self.Input['geometry'] = {}
        self.Input['geometry']['fileString'] = FileName

    # def SetInput(self, Input, Value):
    #     Parameter = Input.split('.') [-1]

    #     if self.Input.get('geometry') is None:
    #         self.Input['geometry'] = {}
    #     self.Input['geometry']['fileString'] = FileName

    def SetBackgroundPlasmaProfiles(self, Voltage = 0):
        self.Input['backgroundPlasmaProfiles'] = {
    'Z': 1.0,
    'amu' : 2.0,
    'biasPotential' : Voltage,
    'Bfield':
        {
        'interpolation' : 0,
        # 'value' : 1.234,
        # 'filename' : "input/bField.nc",
        'r' : self.Input['BField']['r'],
        'z' : self.Input['BField']['z'],
        'y' : self.Input['BField']['y'],
        'rString' : "r",
        'zString' : "z",
        'yString' : "y",
        'fileString' : "bField_created.nc",
        'gridNrString' : "nR",
        'gridNyString' : "nY",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridYString' : "y",
        'gridZString' : "z",
        'radialComponentString' : "br",
        'axialComponentString' : "bz",
        'toroidalComponentString' : "bt"
        },
    'Efield' : 
        {
        'Er' : 0.0,
        'Ez' : 0.0,
        'Et' : 0.0,
        'fileString' : "profiles_created.nc", # LcS.nc
        'gridNrString' : "nR",
        'gridNyString' : "nY",
        'gridNzString' : "nZ",
        'gridRString' : "gridR",
        'gridYString' : "gridY",
        'gridZString' : "gridZ",
        'radialComponentString' : "Er",
        'axialComponentString' : "Ez",
        'toroidalComponentString' : "Et"
        },
    # 'dtsEfield' :
    #     {
    #     'dtsEr' : 0.0,
    #     'dtsEz' : 0.0,
    #     'dtsEt' : 0.0,
    #     'fileString' : "profiles.nc",
    #     'gridNrString' : "n_r_sheathDTS",
    #     'gridNzString' : "n_z_sheathDTS",
    #     'gridRString' : "gridRsheathDTS",
    #     'gridZString' : "gridZsheathDTS",
    #     'sheathDTS' : "sheathDTS",
    #     },
    'Temperature' :
        {
        'ti' : 20.0,
        'te' : 20.0,
        'fileString' : "profiles_created.nc",
        'gridNrString' : "nR",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridZString' : "z",
        'IonTempString' : "ti",
        'ElectronTempString' : "te",
        },
    'Density' :
        {
        'ni' : 1.0E+19,
        'ne' : 1.0E+19,
        'fileString' : "profiles_created.nc",
        'gridNrString' : "nR",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridZString' : "z",
        'IonDensityString' : "ni",
        'ElectronDensityString' : "ne",
        },
    'FlowVelocity' :
        {
        'interpolatorNumber' : 0,
        'flowVr' : 0.0,
        'flowVy' : 0.0,
        'flowVz' : 0.0,
        'fileString' : "profiles_created.nc", # LcS.nc
        'gridNrString' : "nR",
        'gridNyString' : "nY",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridYString' : "y",
        'gridZString' : "z",
        'flowVrString' : "flowVr",
        'flowVzString' : "flowVz",
        'flowVtString' : "flowVt",
        },
    # 'ConnectionLength' : 
    #     {    
    #     'interpolatorNumber' : 2,
    #     'Lc' : 10.0,
    #     's' : 1.0,
    #     'fileString' : "LcS.nc",
    #     'gridNrString' : "nR",
    #     'gridNyString' : "nY",
    #     'gridNzString' : "nZ",
    #     'gridRString' : "gridR",
    #     'gridYString' : "gridY",
    #     'gridZString' : "gridZ",
    #     'LcString' : "Lc",
    #     'SString' : "s",
    #     },
    'gradT' :
        {
        'gradTeR' : 0.0,
        'gradTeY' : 0.0,
        'gradTeZ' : 0.0,
        'gradTiR' : 0.0,
        'gradTiY' : 0.0,
        'gradTiZ' : 0.0,
        'fileString' : "profiles_created.nc",
        'gridNrString' : "nR",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridZString' : "z",
        'gradTiRString' : "gradTiR",
        'gradTiZString' : "gradTiZ",
        'gradTeRString' : "gradTeR",
        'gradTeZString' : "gradTeZ",
        },
    # 'Lc' : 
    #     {    
    #     'value' : 1.0,
    #     'fileString' : "profiles.nc",
    #     'gridNrString' : "nX_Lc",
    #     'gridNzString' : "nY_Lc",
    #     'gridRString' : "gridx_Lc",
    #     'gridZString' : "gridy_Lc",
    #     'variableString' : "Lc",
    #     },
    # 's' : 
    #     {    
    #     'value' : 1.0,
    #     'fileString' : "profiles.nc",
    #     'gridNrString' : "nX_s",
    #     'gridNzString' : "nY_s",
    #     'gridRString' : "gridx_s",
    #     'gridZString' : "gridy_s",
    #     'variableString' : "s",
    #     },
    'Diffusion' :
        {
        'Dperp' : 0.0,
        'fileString' : "profiles_created.nc",
        'gridNrString' : "nR",
        'gridNzString' : "nZ",
        'gridRString' : "r",
        'gridZString' : "z",
        'variableString' : "ni",
        }
}

    def SetParticleSource(self, FileName, nP, Zmax, M, Z):
        if self.Input.get('particleSource') is None:
            self.Input['particleSource'] = {}
        self.Input['particleSource']['ncFileString'] = FileName
        if self.Input.get('impurityParticleSource') is None:
            self.Input['impurityParticleSource'] = {}
        if self.Input['impurityParticleSource'].get('initialConditions') is None:
            self.Input['impurityParticleSource']['initialConditions'] = {}

        Dic = {
            'method':1,
            'nP':nP,
            'sourceStrength' : 1E+19,
            'Z' : 1.0,
            'source_material_Z' : 1,
            'source_material_SurfaceBindingEnergy' : 11.75,

            'initialConditions' :
            {
                # 'x_start' : 0.0,
                # 'y_start' : 0.0,
                # 'z_start' : 0.00001,
                # 'energy_eV_x_start' : 6.0,
                # 'energy_eV_y_start' : 0.0,
                # 'energy_eV_z_start' : 6.0,
                # 'impurity_amu' : 184.0,
                # 'impurity_Z' : 74.0,
                # 'charge' : 0.0,
                'energy_eV' :1,
                'theta':1,
                'phi':1
            },

            'ionization' :  {
            'fileString' : "ADAS_Rates_W.nc",
            'TempGridString' : "n_Temperatures_Ionize",
            'DensGridString' : "n_Densities_Ionize",
            'nChargeStateString' : "n_ChargeStates_Ionize",
            'TempGridVarName' : "gridTemperature_Ionization",
            'DensGridVarName' : "gridDensity_Ionization",
            'CoeffVarName' : "IonizationRateCoeff",
            },

            'recombination' : {
            'fileString' : "ADAS_Rates_W.nc",
            'TempGridString' : "n_Temperatures_Recombine",
            'DensGridString' : "n_Densities_Recombine",
            'nChargeStateString' : "n_ChargeStates_Recombine",
            'TempGridVarName' : "gridTemperature_Recombination",
            'DensGridVarName' : "gridDensity_Recombination",
            'CoeffVarName' : "RecombinationRateCoeff"
            }
        }
        self.Input['impurityParticleSource'].update(Dic)
        self.Input['impurityParticleSource']['initialConditions']['impurity_amu'] = M
        self.Input['impurityParticleSource']['initialConditions']['impurity_Z'] = Zmax
        self.Input['impurityParticleSource']['initialConditions']['charge'] = Z

    def SetTimeStep(self, dt = 1E-6, nPtsPerGyroOrbit = 10000.0, nT = 100000):
        if self.Input.get('timeStep') is None:
            self.Input['timeStep'] = {}
        self.Input['timeStep']['dt'] = dt
        self.Input['timeStep']['nPtsPerGyroOrbit'] = nPtsPerGyroOrbit
        self.Input['timeStep']['nT'] = nT
        self.Input['timeStep']['ionization_nDtPerApply'] = 1
        self.Input['timeStep']['collision_nDtPerApply'] = 5

    def SetSurfaces(self,Emin=0,Emax=1000,nE=1000,Amin=0,Amax=90,nA=90):
        self.Input['surfaces'] = {'useMaterialSurfaces':1,
            'flux' : {
                'nE':nE,
                'E0': Emin,
                'E':Emax,
                'nA' : nA,
                'A0' : Amin,
                'A' : Amax,
                }
            }

    def SetDiagnostics(self):
        self.Input['diagnostics'] ={
            'leakZ' : 1.0,
            'trackSubSampleFactor' : 1000,
            'netx0' : -0.03,
            'netx1' : 0.03,
            'nX' : 100,
            'nety0' : -0.03,
            'nety1' : 0.03,
            'nY' : 120,
            'netz0' : -0.015,
            'netz1' : 0.03,
            'nZ' : 150,
            'densityChargeBins' : 6,
            }

    def SetForceEvaluation(self):
        self.Input['forceEvaluation'] = {
            'X0' : -0.03,
            'X1' : 0.03,
            'Y0' : -0.03,
            'Y1' : 0.03,
            'Z0' : -0.015,
            'Z1' : 0.03,
            'nR' : 176,
            'nY' : 0,
            'nZ' : 372,
            'particleEnergy' : 10.0,
        }

    def SetFlags(self):
        code_flags = ''
        for k,v in self.Input['flags'].items():
            code_flags = code_flags + " -D{}={}".format(k,v)

        return code_flags

    def CompileGITR(self):
        print(self.SetFlags())
