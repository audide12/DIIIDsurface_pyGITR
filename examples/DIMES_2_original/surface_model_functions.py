import numpy as np
import io, libconf
import netCDF4
from netCDF4 import Dataset
import os



def sign (p1_x, p1_y, p2_x, p2_y, p3_x,p3_y):
    return (p1_x - p3_x) * (p2_y - p3_y) - (p2_x - p3_x) * (p1_y - p3_y)


def PointInTriangle (pt_x, pt_y, v1_x, v1_y, v2_x, v2_y, v3_x, v3_y):

    d1 = sign(pt_x,pt_y,v1_x,v1_y,v2_x,v2_y)
    d2 = sign(pt_x,pt_y,v2_x,v2_y,v3_x,v3_y)
    d3 = sign(pt_x,pt_y,v3_x,v3_y,v1_x,v1_y)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not(has_neg and has_pos)


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
    initial_token_flux = Flux_C_Background[0]#1.0e19  # tunable parameter

    Flux_proportionality = {}

    Total_Area = 0
    for k in range(dim):
        Total_Area += area[k]
    for Z in Conc.keys(): 
        Flux_proportionality[Z] = initial_token_flux*Total_Area*Delta_t_gitr/N_GITR

    for Z in Conc.keys():
        Flux_proportionality[Z] = 0
        for k in range(dim):
            Flux_proportionality[Z] += initial_token_flux*area[k]*Delta_t_gitr/N_GITR

    Surface_time = np.full((1,1),0.0)
    Surface_number = np.array(range(dim))


    ncFile = netCDF4.Dataset('surface_evolution_C_Si.nc', 'w', format='NETCDF4')

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

    os.system("mv surface_evolution_C_Si.nc /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/")



Flux_H = 1.0e20
alpha_c = 0.02      # Carbon concentration in the background plasma
Flux_C = alpha_c*Flux_H

Delta_t = 0.1 # in seconds 
Delta_t_gitr = 1e-7
Delta_implant = 250e-6 # enter parameter value and units
amu_C = 12 #for carbon
amu_W = 184 #for tungsten
k_B = 1.38 * 10**(-23)   # in m^2 kg s^-2 K^-1  Boltzmann constant
M_H = 1.67*10**(-27)  # in kg
M_C = 12*M_H  
amu_Si = 28 #for silicon
n_atom = 6e22 # average number density in cm^-3
n_atom_C = 1.12e29 # in m^-3
n_atom_Si = 5e28 # in m^-3
n_atom_SiC_crystal = 4.8e28 # in m^-3

Delta_implant_amorphous = 250e-6 # in metres

weight_gitr = Delta_t/Delta_t_gitr
Stopping_criteria = 0.1 # for C_C and C_W

# Processing of background plasma 

BackgroundPlasmaFile='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/profiles_created.nc'
BackgroundPlasma = netCDF4.Dataset(BackgroundPlasmaFile)

T_e = np.array(BackgroundPlasma['te'])
T_i = np.array(BackgroundPlasma['ti'])
R = np.array(BackgroundPlasma['r'])
Z = np.array(BackgroundPlasma['z'])
ni = np.array(BackgroundPlasma['ni'])

GeomFile = "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/gitrGeom.cfg"

#x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surf,indir,Atomic_no = loadCFG(geomFile=GeomFile)
x1,x2,x3,y1,y2,y3,z1,z2,z3,area,surf,Atomic_no,a,b,c,d,in_direction,plane_norm = getGeom(GeomFile)
Zs = []
Surfaces = []
idx = np.arange(0,len(surf))
for surface,z,i in zip(surf,Atomic_no,idx):
    if surface!=0:
        
        Zs.append(z)
        Surfaces.append(i)
Zs = np.unique(Zs)


Flux_H_Background = np.zeros((len(Surfaces),1))
Flux_C_Background = np.zeros((len(Surfaces),1))
YHtoSiC_Flux_H_in = np.zeros((len(Surfaces),1))
YCtoSiC_Flux_C_in = np.zeros((len(Surfaces),1))
YHtoSi_Flux_H_in  = np.zeros((len(Surfaces),1))
YHtoC_Flux_H_in   = np.zeros((len(Surfaces),1))
YCtoSi_Flux_C_in  = np.zeros((len(Surfaces),1))
YCtoC_Flux_C_in   = np.zeros((len(Surfaces),1))

beta_depC1 = np.zeros((len(Surfaces),1))
beta_depC2 = np.zeros((len(Surfaces),1))

Numbers = np.zeros((len(Surfaces),1))
sr_object = Sputtering_and_reflection()


for j in Surfaces:
    for idr,r in enumerate(R): 
        for idz,z in enumerate(Z):
            
            # print(z2[j])
        
            if(PointInTriangle (R[idr], Z[idz], x1[j], z1[j], x2[j], z2[j], x3[j], z3[j])):
            
                #print("done")
                Energy_ion_H = 2*T_i[idz,idr]+3*1*T_e[idz,idr]
                Energy_ion_C = 2*T_i[idz,idr]+3*6*T_e[idz,idr]  # check this
            
                
                Gamma_H = 0.61 * ni[idz,idr] * np.sqrt(k_B*T_e[idz,idr]*11600/M_H) * np.sin(1.5*np.pi/180) # ni(m^-3),Te(eV),1.5 degrees DiMES angle, 1eV =  11600 K, in m^-2 s^-1
                Gamma_C = alpha_c * Gamma_H  # in m^-2 s^-1
            
                #print(Gamma_H)
            
                Flux_H_Background[j,0] = Flux_H_Background[j,0] + Gamma_H
                Flux_C_Background[j,0] = Flux_C_Background[j,0] + Gamma_C
                
                sr_object = Sputtering_and_reflection()
                YHtoSiC_Flux_H_in[j,0] = YHtoSiC_Flux_H_in[j,0] + sr_object.Calculate_PhysicalSputteringParameters('H','SiC',Energy_ion_H,60)*Gamma_H
                YCtoSiC_Flux_C_in[j,0] = YCtoSiC_Flux_C_in[j,0] + sr_object.Calculate_PhysicalSputteringParameters('C','SiC',Energy_ion_C,60)*Gamma_C
                YHtoSi_Flux_H_in[j,0]  = YHtoSi_Flux_H_in[j,0]  + sr_object.Calculate_PhysicalSputteringParameters('H','Si',Energy_ion_H,60)*Gamma_H
                YHtoC_Flux_H_in[j,0]   = YHtoC_Flux_H_in[j,0]   + sr_object.Calculate_PhysicalSputteringParameters('H','C',Energy_ion_H,60)*Gamma_H
                YCtoSi_Flux_C_in[j,0]  = YCtoSi_Flux_C_in[j,0]  + sr_object.Calculate_PhysicalSputteringParameters('C','Si',Energy_ion_C,60)*Gamma_C
                YCtoC_Flux_C_in[j,0]   = YCtoC_Flux_C_in[j,0]   + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_ion_C,60)*Gamma_C
                
                sr_object = Sputtering_and_reflection()
                beta_depC1[j,0] = beta_depC1[j,0] + (1-sr_object.Calculate_ReflectionCoefficients('C','C',Energy_ion_C))*Gamma_C
                beta_depC2[j,0] = beta_depC2[j,0] + (1-sr_object.Calculate_ReflectionCoefficients('C','C',Energy_ion_C))*Gamma_C
                
                Numbers[j,0] = Numbers[j,0] + 1
                
                
Flux_H_Background = Flux_H_Background/Numbers
Flux_C_Background = Flux_C_Background/Numbers
YHtoSiC_Flux_H_in = YHtoSiC_Flux_H_in/Numbers
YCtoSiC_Flux_C_in = YCtoSiC_Flux_C_in/Numbers
YHtoSi_Flux_H_in  = YHtoSi_Flux_H_in/Numbers
YHtoC_Flux_H_in   = YHtoC_Flux_H_in/Numbers
YCtoSi_Flux_C_in  = YCtoSi_Flux_C_in/Numbers
YCtoC_Flux_C_in   = YCtoC_Flux_C_in/Numbers
beta_depC1 = beta_depC1/Numbers
beta_depC2 = beta_depC2/Numbers

# Sputtering_yield_H_to_C = 0.005  # max 0.005
# Sputtering_yield_C_to_C = 0.22   # 0.22 treated almost constant

# Sputtering_yield_H_to_W = 0.002  # max 0.002 
# Sputtering_yield_C_to_W = 0.5  # 0.5 treated almost constant

# Reflection_yield_C_to_C = 0.005    # max 0.005 min 0.001  steady-state :  0.9
# Reflection_yield_C_to_W = 0.75  # max 0.95 min 0.67    steady-state :  0.005


# Sputtering_yield_H_to_SiC = 0.9459 # average from Eckstein data
# Sputtering_yield_C_to_SiC = 0.5 * Sputtering_yield_C_to_C # approximation
# Sputtering_yield_H_to_Si = 0.01 # average value from the sputtering plots


N_GITR = 10000 # number of GITR particles
