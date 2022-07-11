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
                    'USEPARDIFFUSION':0,
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
                    'FIXEDSEEDS':1,
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

    def WriteInputFile(self, FileName='gitrInput.cfg', Folder='', OverWrite=False):
        FileName = os.path.abspath(os.path.join(Folder,FileName))
        print('Writing input config into {} ...'.format(FileName))

        if os.path.exists(FileName):
            if not OverWrite and not click.confirm('File {} exists.\n Do you want to overwrite it?'.format(FileName), default=True):
                return

        with io.open(FileName,'w') as f:
            libconf.dump(self.Input,f)

    def SetBField(self,Br=None, By=None, Bt=None, B0=None, theta=None, phi=None, Degree=True):
        if self.Input.get('BField') is None:
            self.Input['BField'] = {}

        if Br is not None and By is not None and Bt is not None:
            self.Input['BField']['r'] = Br
            self.Input['BField']['z'] = By
            self.Input['BField']['y'] = Bz
            self.Input['backgroundPlasmaProfiles']['BField']['r'] = Br
            self.Input['backgroundPlasmaProfiles']['BField']['z'] = By
            self.Input['backgroundPlasmaProfiles']['BField']['y'] = Bz

        elif B0 is not None and theta is not None and phi is not None:
            theta = 2 #degree by default, set Degree=False in RotateVector for radian
            Axisroty = [0,1,0]
            Axisrotz = [0,0,1]
            Btot = [B0,0,0]
            B1 = RotateVector(Btot, Axisroty, theta, Degree)
            B  = RotateVector(B1, Axisrotz, phi, Degree)
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
    'r' : self.Input['BField']['r'],
    'z' : self.Input['BField']['z'],
    'y' : self.Input['BField']['y'],
    'rString' : "r",
    'zString' : "z",
    'yString' : "y"
},
    'Temperature' :
        {
        'ti' : 20.0,
        'te' : 20.0,
        },
    'Density' :
        {
        'ni' : 1.0E+19,
        'ne' : 1.0E+19,
        },
        'FlowVelocity' :
        {
        'interpolatorNumber' : 0,
        'flowVr' : 0.0,
        'flowVy' : 0.0,
        'flowVz' : 0,
        },
        'gradT' :
        {
        'gradTeR' : -0.0,
        'gradTeY' : 0.0,
        'gradTeZ' : 0.0,
        'gradTiR' : -0.0,
        'gradTiY' : 0.0,
        'gradTiZ' : 0.0,
        },
    'Diffusion' :
        {
        'Dperp' : 0.0
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



        Dic ={   'method':1,
	'sourceStrength' : 1E+19,
	'Z' : 1.0,
    'source_material_Z' : 1,
    'nP':nP,
    'initialConditions' :
    {
        'energy_eV' :1,
        'theta':1,
        'phi':1

    },

        'ionization' :  {'fileString' : "ADAS_Rates_W.nc",
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
    'leakZ':1.0,
    'trackSubSampleFactor':1000
}

    def SetFlags(self):
        code_flags = ''
        for k,v in self.Input['flags'].items():
            code_flags = code_flags + " -D{}={}".format(k,v)

        return code_flags

    def CompileGITR(self):
        print(self.SetFlags())
















