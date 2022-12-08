import numpy as np
import io, libconf
import netCDF4
from netCDF4 import Dataset
import os

Flux_H = 1.0e20
alpha_c = 0.02       # Carbon concentration in the background plasma
Flux_C = alpha_c*Flux_H

Delta_t = 0.1 # in seconds
Delta_t_gitr = 1e-7
Delta_implant = 1e-5 # enter parameter value and units
amu_C = 12 #for carbon
amu_W = 184 #for tungsten
<<<<<<< HEAD
amu_Si = 28 #for silicon
=======
>>>>>>> 742051bb5c26808a9ff7ec7f08154e22d50ea9c2

n_atom = 6e22 # average number density
weight_gitr = Delta_t/Delta_t_gitr
Stopping_criteria = 0.1 # for C_C and C_W

Sputtering_yield_H_to_C = 0.005  # max 0.005
Sputtering_yield_C_to_C = 0.22   # 0.22 treated almost constant

Sputtering_yield_H_to_W = 0.002  # max 0.002 
Sputtering_yield_C_to_W = 0.5  # 0.5 treated almost constant

Reflection_yield_C_to_C = 0.005    # max 0.005 min 0.001  steady-state :  0.9
Reflection_yield_C_to_W = 0.75  # max 0.95 min 0.67    steady-state :  0.005

<<<<<<< HEAD
Sputtering_yield_H_to_SiC = 0.9459 # average from Eckstein data
Sputtering_yield_C_to_SiC = 0.5 * Sputtering_yield_C_to_C # approximation
Sputtering_yield_H_to_Si = 0.01 # average value from the sputtering plots

=======
>>>>>>> 742051bb5c26808a9ff7ec7f08154e22d50ea9c2
N_GITR = 10000 # number of GITR particles

def initDict(species,myDict):
    for specie in species:
        myDict[specie] = []
    return myDict


def getE(amu,vx,vy,vz):
    return np.array(0.5*amu*1.66e-27*(vx**2 + vy**2 + vz**2)/1.602e-19)

def getGeom(File):

    with io.open(File) as f:
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
    a = np.array(config.geom.a)
    b = np.array(config.geom.b)
    c = np.array(config.geom.c)
    d = np.array(config.geom.d)
    in_direction = np.array(config.geom.inDir)
    plane_norm = np.array(config.geom.plane_norm)

    return x1,x2,x3,y1,y2,y3,z1,z2,z3,area,\
        surf,Z,a,b,c,d,in_direction,plane_norm 


def makeInitNC(dim,area,Conc):

    N_GITR = 10000 # number of GITR particles
    initial_token_flux = 1.0e19  # tunable parameter

    Flux_proportionality = {}
<<<<<<< HEAD
    Total_Area = 0
    for k in range(dim):
        Total_Area += area[k]
    for Z in Conc.keys(): 
        Flux_proportionality[Z] = initial_token_flux*Total_Area*Delta_t_gitr/N_GITR
=======
    for Z in Conc.keys():
        Flux_proportionality[Z] = 0
        for k in range(dim):
            Flux_proportionality[Z] += initial_token_flux*area[k]*Delta_t_gitr/N_GITR
>>>>>>> 742051bb5c26808a9ff7ec7f08154e22d50ea9c2

    Surface_time = np.full((1,1),0.0)
    Surface_number = np.array(range(dim))

<<<<<<< HEAD
    ncFile = netCDF4.Dataset('surface_evolution_C_Si.nc', 'w', format='NETCDF4')
=======
    ncFile = netCDF4.Dataset('surface_evolution_C_W.nc', 'w', format='NETCDF4')
>>>>>>> 742051bb5c26808a9ff7ec7f08154e22d50ea9c2
    s_number_dim = ncFile.createDimension('surface_dim', dim) # surface number dimension
    s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

    s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
    s_time = ncFile.createVariable('time', np.float32, ('time_dim',))

    s_concentration = {}
    flux_proportionality = {}
    for Z in Conc.keys():
        s_concentration[Z] = ncFile.createVariable('surface_concentration_{}'.format(Z), np.float64, ('surface_dim','time_dim'))
        flux_proportionality[Z] = ncFile.createVariable('Flux_Conversion_{}'.format(Z),np.float64,('time_dim'))

    s_number[:] = np.linspace(1,dim,dim)
    s_time[:] = Surface_time
    for Z in Conc.keys():
        s_concentration[Z][:,:] = Conc[Z]
        flux_proportionality[Z][:] = Flux_proportionality[Z]

    ncFile.close()
<<<<<<< HEAD
    os.system("mv surface_evolution_C_Si.nc /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/")
=======
    os.system("mv surface_evolution_C_W.nc /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/")
>>>>>>> 742051bb5c26808a9ff7ec7f08154e22d50ea9c2
