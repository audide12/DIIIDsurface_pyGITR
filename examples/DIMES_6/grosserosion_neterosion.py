#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:58:04 2023

@author: de
"""

# run surface_model_functions.py
#Reading geometry files

GeomFile = "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/gitrGeom.cfg"
x1,x2,x3,y1,y2,y3,z1,z2,z3,area,surf,Atomic_no,a,b,c,d,in_direction,plane_norm = getGeom(GeomFile)
#x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surf,indir,Atomic_no = loadCFG(geomFile=GeomFile)

# Initialize the surface_evolution netcdf file
# Only care about surfaces


Zs = []

Surfaces = []
idx = np.arange(0,len(surf))
for surface,z,i in zip(surf,Atomic_no,idx):
    if surface!=0:
        Zs.append(z)
        Surfaces.append(i)
Zs = np.unique(Zs)
Zs = np.append(Zs,14)  # Adding Silicon

print(Zs,"make up the", len(Surfaces),"surface mesh elements")
# print(Surfaces)

#Reading the surface features from the surface evolution netcdf file

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/surface_evolution_C_Si.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")

# Record concentrations of all surface elements and their initial Z
Concentration = {}
Flux_proportionality = {}
for z in Zs:
    Concentration[z] = SurfaceConcentrationData['surface_concentration_{}'.format(z)][:,:]
    Flux_proportionality[z] = SurfaceConcentrationData['Flux_Conversion_{}'.format(z)][:]

Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)




#%% 

total_runs = 19 ## changed this

time_index = 0
Gamma_C_grosserosion  = np.zeros((len(Surfaces),total_runs))
Gamma_C_neterosion    = np.zeros((len(Surfaces),total_runs))
Gamma_Si_grosserosion = np.zeros((len(Surfaces),total_runs))
Gamma_Si_neterosion   = np.zeros((len(Surfaces),total_runs))
Time                  = np.zeros((1,total_runs))

for l in range(total_runs):
    number = l+2# changed this
    
    time_index = time_index + 1
    
    
    FileNameHistoryCarbon='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/output_C_'+str(number)+'/positions.nc'
    
    PositionData = netCDF4.Dataset(FileNameHistoryCarbon)

    surfacehit_C = np.array(PositionData['surfaceHit'])
    surface_vx_C = np.array(PositionData['vx'])
    surface_vy_C = np.array(PositionData['vy'])
    surface_vz_C = np.array(PositionData['vz'])
    Energy_particles_C = np.array(0.5*amu_C*1.66e-27*(surface_vx_C**2 + surface_vy_C**2 + surface_vz_C**2)/1.602e-19) # make sure that this energy is in eV
    Angles_particles_C = np.arctan(surface_vx_C/surface_vz_C)*(180/np.pi)   # in degrees

    # Reading position files of Silicon

    FileNameHistorySilicon='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/output_Si_'+str(number)+'/positions.nc'
    
    PositionData = Dataset(FileNameHistorySilicon, "r", format="NETCDF4")

    surfacehit_Si = np.array(PositionData['surfaceHit'])
    surface_vx_Si = np.array(PositionData['vx'])
    surface_vy_Si = np.array(PositionData['vy'])
    surface_vz_Si = np.array(PositionData['vz'])
    Energy_particles_Si = np.array(0.5*amu_Si*1.66e-27*(surface_vx_Si**2 + surface_vy_Si**2 + surface_vz_Si**2)/1.602e-19) # make sure that this energy is in eV
    Angles_particles_Si = np.arctan(surface_vx_Si/surface_vz_Si)*(180/np.pi)    #  in degrees
    
    Gamma_C_redep = np.zeros((len(Surfaces),1))
    Y_CSiC_Gamma_C_redep = np.zeros((len(Surfaces),1))
    Y_CSi_Gamma_C_redep = np.zeros((len(Surfaces),1))
    Y_CC_Gamma_C_redep = np.zeros((len(Surfaces),1))


    for i in range(len(Energy_particles_C)):
        if surfacehit_C[i] != -1:

            surface_index = int(surfacehit_C[i])
            sr_object = Sputtering_and_reflection()

            for j in Surfaces:
                if j == surface_index:

                    #print("yes")
                    Flux_C_local = Flux_proportionality[6][time_index]/(Delta_t_gitr*area[surface_index])  # we start with initial weight 1 (uniform)
                    #print(Flux_C_local)
                    
                    Gamma_C_redep[surface_index] = Gamma_C_redep[surface_index] + Flux_C_local  # check this
                    Y_CSiC_Gamma_C_redep[surface_index] = Y_CSiC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','SiC',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
                    Y_CC_Gamma_C_redep[surface_index] = Y_CC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
                    Y_CSi_Gamma_C_redep[surface_index] = Y_CSi_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','Si',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local

    Gamma_Si_redep = np.zeros((len(Surfaces),1))
    Y_SiSi_Gamma_Si_redep = np.zeros((len(Surfaces),1))
    Y_SiSiC_Gamma_Si_redep = np.zeros((len(Surfaces),1))
    Y_Si_C_Gamma_Si_redep = np.zeros((len(Surfaces),1))    # Y_Si_to_C
            
            
    for i in range(len(Energy_particles_Si)):
        if surfacehit_Si[i] != -1:
            surface_index = int(surfacehit_Si[i])

            sr_object = Sputtering_and_reflection()
            
            for j in Surfaces:
                if j == surface_index:

                    Flux_Si_local = Flux_proportionality[14][time_index]/(Delta_t_gitr*area[surface_index]) # we start with initial weight 1 (uniform)
                    
                    Gamma_Si_redep[surface_index] = Gamma_Si_redep[surface_index] + Flux_Si_local  # check this
                    Y_SiSi_Gamma_Si_redep[surface_index] = Y_SiSi_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','Si',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
                    Y_SiSiC_Gamma_Si_redep[surface_index] = Y_SiSiC_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','SiC',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
                    Y_Si_C_Gamma_Si_redep[surface_index] = Y_Si_C_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','C',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local



    beta_eroSi1 = YHtoSi_Flux_H_in
    beta_eroSi2 = YCtoSi_Flux_C_in 
    beta_eroSi3 = Y_CSi_Gamma_C_redep
    beta_eroSi4 = Y_SiSi_Gamma_Si_redep

    beta_eroC1 = YHtoC_Flux_H_in
    beta_eroC2 = YCtoC_Flux_C_in
    beta_eroC3 = Y_CC_Gamma_C_redep
    beta_eroC4 = Y_Si_C_Gamma_Si_redep

    beta_SiC1 = YHtoSiC_Flux_H_in + YCtoSiC_Flux_C_in
    beta_SiC = beta_SiC1 + Y_CSiC_Gamma_C_redep + Y_SiSiC_Gamma_Si_redep
       
    Gamma_SiC_ero_global = np.reshape(Concentration[20][:,time_index],(len(Surfaces),1))*beta_SiC

    Gamma_Si_ero_exclusive = np.reshape(Concentration[14][:,time_index],(len(Surfaces),1))*(beta_eroSi1 + beta_eroSi2 + beta_eroSi3 + beta_eroSi4)
    Gamma_Si_ero_global = Gamma_SiC_ero_global + Gamma_Si_ero_exclusive


    Gamma_C_ero_exclusive = np.reshape(Concentration[6][:,time_index],(len(Surfaces),1))*(beta_eroC1 + beta_eroC2 + beta_eroC3 + beta_eroC3)
    Gamma_C_ero_global = Gamma_SiC_ero_global + Gamma_C_ero_exclusive
    
    Gamma_C_ero_global[np.argwhere(np.isnan(Gamma_C_ero_global))] = Gamma_C_ero_global[np.argwhere(np.isnan(Gamma_C_ero_global))+1] #new
    
    Gamma_C_dep_global = Gamma_C_redep + np.reshape(Concentration[14][:,time_index],(len(Surfaces),1))*Flux_C_Background + np.reshape(Concentration[6][:,time_index],(len(Surfaces),1))*beta_depC1 +  np.reshape(Concentration[20][:,time_index],(len(Surfaces),1))*beta_depC2     
    Gamma_Si_dep_global = Gamma_Si_redep
    
    
    # Gamma_C_grosserosion  = Gamma_C_grosserosion +  Gamma_C_ero_global
    # Gamma_Si_grosserosion = Gamma_Si_grosserosion + Gamma_Si_ero_global
    
    # Gamma_C_neterosion  = Gamma_C_neterosion +  (Gamma_C_ero_global - Gamma_C_redep)
    # Gamma_Si_neterosion = Gamma_Si_neterosion + (Gamma_Si_ero_global - Gamma_C_redep)
    
    Gamma_C_grosserosion[:,time_index-1]  = Gamma_C_ero_global[:,0]
    Gamma_Si_grosserosion[:,time_index-1] = Gamma_Si_ero_global[:,0]
    
    Gamma_C_neterosion[:,time_index-1]    = (Gamma_C_ero_global - Gamma_C_redep)[:,0]
    Gamma_Si_neterosion[:,time_index-1]   = (Gamma_Si_ero_global - Gamma_Si_redep)[:,0]
    
    Time[0,time_index-1] = Surface_time[time_index]
    print(Surface_time[time_index],"  seconds")
    
#%%
import matplotlib.pyplot as plt    
surface_in_question = 40

time_fit = Time[0,:]
gross_fit = 8.15e19+1e19*np.exp(-1.04*time_fit)
net_fit = 3.8e19+1e19*(1-0.03*time_fit)
plt.figure()
plt.plot(Time[0,:],Gamma_Si_grosserosion[surface_in_question,:],marker='^',label='GITR + Surface model')
plt.plot(time_fit,net_fit,marker='+',label='FIT')
plt.ylabel("Erosion (m^-2 s^-1)")
plt.xlabel("Time (s)")
#plt.xlim(0.0,0.12)

plt.legend()
plt.title("Si Gross Erosion in Surface Element %d" % surface_in_question)
plt.show()  

#%%
import matplotlib.pyplot as plt    

surface_in_question = 47
plt.figure()

plt.plot(Time[0,:],Gamma_Si_neterosion[surface_in_question,:],marker='^',label='GITR + Surface model')
#plt.plot(Time_experiment_SiC_crystal-1.9, Si_grosserosion_experiment_SiC_crystal, label='DiMES Experiment - SiC crystalline')
#plt.plot(Time_experiment1, Si_grosserosion_experiment1, label='DiMES Experiment - - SiC amorphous')

plt.ylabel("Erosion (m^-2 s^-1)")
plt.xlabel("Time (s)")

plt.legend()
plt.title("Si Gross Erosion in Surface Element %d" % surface_in_question)
plt.show()


#%%
import pandas as pd

df1 = pd.read_excel(r'/Users/de/Downloads/176508_Sigrosserosion_DiMESTV_forStefanDmitry_v2+DR.xlsx', sheet_name='176488')

Time_experiment_SiC_crystal = np.array(df1.loc[:,"Time (ms)"])
Time_experiment_SiC_crystal = Time_experiment_SiC_crystal - 1.51622
Si_grosserosion_experiment_SiC_crystal = np.array(df1.loc[:,"Erosion (cm-2 s-1)"])
Si_grosserosion_experiment_SiC_crystal = Si_grosserosion_experiment_SiC_crystal*1e4
df = pd.read_excel('/Users/de/Downloads/176508_Sigrosserosion_DiMESTV_forStefanDmitry_v2+DR.xlsx',sheet_name='176508')
Time_experiment = np.array(df.loc[:,"Time (ms)"])
Time_experiment1 = Time_experiment - 1.51622
Si_grosserosion_experiment = np.array(df.loc[:,"Erosion (cm-2 s-1)"])
Si_grosserosion_experiment1 = Si_grosserosion_experiment*1e4

#%%
import matplotlib.pyplot as plt    
surface_in_question = 190
plt.figure()

plt.plot(Time[0,:],Gamma_Si_neterosion[surface_in_question,:],marker='^',label='Net')
plt.plot(Time[0,:],Gamma_Si_grosserosion[surface_in_question,:],marker='^',label='Gross')

plt.ylabel("Erosion (m^-2 s^-1)")
plt.xlabel("Time (s)")

plt.legend()
plt.title("Si Net Erosion in Surface Element %d" % surface_in_question)
plt.show()






