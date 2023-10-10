#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:36:51 2023

@author: de
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 10:45:19 2023

@author: de
"""

g = GeomSetup('DiMES_SiC_amorphous.msh', Verbose=True)
# Plot mesh for some groups
g.Plot_Geom(["DiMES"])

g.ShowCentroids_Annotated("DiMES")

# Measurement after 4.5 seconds of shots
GeomFile = "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_6/input/gitrGeom.cfg"
x1,x2,x3,y1,y2,y3,z1,z2,z3,area,surf,Atomic_no,a,b,c,d,in_direction,plane_norm = getGeom(GeomFile)

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

mesh_tracking = [43,41,42,30,26,33,29,47,46,45,44]

mesh_tracking_paper = [1,2,3,4,5,6,7,8,9,10,11,12,13]

total_runs = 9

time_index = 0
total_time = 0.0

Gamma_C_grosserosion  = np.zeros((len(Surfaces),total_runs))
Gamma_C_neterosion    = np.zeros((len(Surfaces),total_runs))
Gamma_Si_grosserosion = np.zeros((len(Surfaces),total_runs))
Gamma_Si_neterosion   = np.zeros((len(Surfaces),total_runs))
Time                  = np.zeros((1,total_runs))

Gamma_C_grosserosion_absolute  = np.zeros((len(Surfaces),1))
Gamma_C_neterosion_absolute    = np.zeros((len(Surfaces),1))
Gamma_Si_grosserosion_absolute = np.zeros((len(Surfaces),1))
Gamma_Si_neterosion_absolute   = np.zeros((len(Surfaces),1))

Gamma_C_ero_global_last =  np.zeros((len(Surfaces),1))
Gamma_Si_ero_global_last =  np.zeros((len(Surfaces),1))

for l in range(total_runs):
    number = l+1
        
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
    
    Gamma_C_grosserosion[:,time_index]  = Gamma_C_ero_global[:,0]
    Gamma_Si_grosserosion[:,time_index] = Gamma_Si_ero_global[:,0]
    
    Gamma_C_neterosion[:,time_index]    = (Gamma_C_ero_global_last - Gamma_C_redep)[:,0]
    Gamma_Si_neterosion[:,time_index]   = (Gamma_Si_ero_global_last - Gamma_Si_redep)[:,0]
    
    Gamma_C_ero_global_last[:,0] = Gamma_C_ero_global[:,0]
    Gamma_Si_ero_global_last[:,0] = Gamma_Si_ero_global[:,0]
    
    Time[0,time_index] = Surface_time[time_index]
    
    
    Delta_tt = Surface_time[time_index+1] - Surface_time[time_index]
    print(Surface_time[time_index],"  in seconds")
    
    Gamma_C_grosserosion_absolute[:,0]  = Gamma_C_grosserosion_absolute[:,0] + Delta_tt* Gamma_C_ero_global[:,0]    # NOT RELEVANT
    Gamma_Si_grosserosion_absolute[:,0] = Gamma_Si_grosserosion_absolute[:,0] + Delta_tt* Gamma_Si_ero_global[:,0]  # NOT RELEVANT
    
    Gamma_C_neterosion_absolute[:,0]    = Gamma_C_neterosion_absolute[:,0] + Delta_tt* (Gamma_C_neterosion[:,time_index])
    Gamma_Si_neterosion_absolute[:,0]   = Gamma_Si_neterosion_absolute[:,0] + Delta_tt* (Gamma_Si_neterosion[:,time_index])
    
    
    Time[0,time_index] = Surface_time[time_index]
    
    time_index = time_index + 1
    total_time = total_time + Delta_tt


# Gamma_C_neterosion_absolute = Gamma_C_neterosion_absolute/total_time
# Gamma_Si_neterosion_absolute = Gamma_Si_neterosion_absolute/total_time  

# Gamma_C_grosserosion_absolute = Gamma_C_grosserosion_absolute/total_time
# Gamma_Si_grosserosion_absolute = Gamma_Si_grosserosion_absolute/total_time    

#%%    
Gamma_Si_neterosion_comparison = np.zeros((len(mesh_tracking),1))
Gamma_Si_grosserosion_comparison = np.zeros((len(mesh_tracking),1))
index = 0    
for k in mesh_tracking:
    Gamma_Si_grosserosion_comparison[index,0] = Gamma_Si_grosserosion_absolute[k,0]
    Gamma_Si_neterosion_comparison[index,0] = Gamma_Si_neterosion_absolute[k,0]
    index = index+1
    
#mesh_tracking = [43,41,42,30,17,33,29,47,46,45,44]

mesh_tracking = [43,41,42,30,35,34,29,47,46,45,44]    
surface_index = np.array([1,2,3,4,5,6,7,8,9,10,11])

#mesh_tracking = [43,41,42,30,31,30,29]    
#surface_index = np.array([1,2,3,4,5,6,7])

surface_index_ticks = ['1','2','3','4','5','6','7','8','9','11','12'] 

#surface_index_ticks = ['1','2','3','4','5','6','7'] 
import matplotlib.pyplot as plt


fig, ax = plt.subplots(1,1) 
ax.bar(surface_index,(Gamma_Si_neterosion_comparison/1e21).reshape(11), width = 0.2, color ='maroon')

# Set number of ticks for x-axis
ax.set_xticks(surface_index)
# Set ticks labels for x-axis
ax.set_xticklabels(surface_index_ticks, fontsize=18)
plt.yticks(fontsize=10)
plt.xticks(fontsize=10)
plt.ylabel(r"Erosion $\times(10^{17} cm^{-2}) $",fontsize=18,color='blue')
plt.xlabel("Surface Index",fontsize=18,color='blue')
plt.title("Si Erosion comparison with RBS measurements",fontsize=20,color='blue')
plt.show() 

#%%
total_runs = 19

time_index = 0
Gamma_C_grosserosion  = np.zeros((len(Surfaces),total_runs))
Gamma_C_neterosion    = np.zeros((len(Surfaces),total_runs))
Gamma_Si_grosserosion = np.zeros((len(Surfaces),total_runs))
Gamma_Si_neterosion   = np.zeros((len(Surfaces),total_runs))
Time                  = np.zeros((1,total_runs))

for l in range(total_runs):
    number = l+1    
    
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
    
    Gamma_C_grosserosion[:,time_index]  = Gamma_C_ero_global[:,0]
    Gamma_Si_grosserosion[:,time_index] = Gamma_Si_ero_global[:,0]
    
    #Gamma_C_neterosion[:,time_index]    = (Gamma_C_ero_global - Gamma_C_redep)[:,0]
    #Gamma_Si_neterosion[:,time_index]   = (Gamma_Si_ero_global - Gamma_Si_redep)[:,0]
    
    Time[0,time_index] = Surface_time[time_index]
    print(Surface_time[time_index],"  seconds")
    time_index = time_index + 1
    
#%% Averaging over the DiMES cap    

mesh_tracking = [47,46,45,44]

Gamma_C_grosserosion_DiMES  = np.zeros((1,total_runs))
Gamma_Si_grosserosion_DiMES = np.zeros((1,total_runs))
total_area = 0.0

for surface_index in Surfaces:#mesh_tracking:Surfaces:
    Gamma_C_grosserosion_DiMES[0,:]  = Gamma_C_grosserosion_DiMES[0,:] + area[surface_index]*Gamma_C_grosserosion[surface_index,:]
    Gamma_Si_grosserosion_DiMES[0,:]  = Gamma_Si_grosserosion_DiMES[0,:] + area[surface_index]*Gamma_Si_grosserosion[surface_index,:]
    total_area = total_area + area[surface_index]
    
Gamma_C_grosserosion_DiMES = Gamma_C_grosserosion_DiMES/total_area
Gamma_Si_grosserosion_DiMES = Gamma_Si_grosserosion_DiMES/total_area

import matplotlib.pyplot as plt    

plt.figure()
plt.plot(Time[0,:],Gamma_Si_grosserosion_DiMES[0,:],marker='^',label='GITR + Surface model')
plt.ylabel("Erosion (m^-2 s^-1)")
plt.xlabel("Time (s)")

plt.legend()
plt.title("Si Gross Erosion")
plt.show()