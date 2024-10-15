#%% Startup
import pyGITR
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/make_particleSource.py')
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/Physical_Sputtering.py')
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/Geom.py')
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/Particles.py')
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/particleSource_functions.py')
runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/pyGITR/Eckstein_angle_sputtering_parameters.py')

#%% Subsequent surface runs
import time

start = time.time()
alpha_c = 0.02
# Stopping_criteria = 0.075
# Delta_implant_amorphous = 30e-9
#alpha_amor = 0.18

runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')

runcell('Reading geometry files', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_initial_run.py')

runcell('Initiallize all surfaces with concentrations equal to the Z at that surface', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_initial_run.py')

for run_no in [2,3,4,5,6]:#,6,7,8]:
    if run_no == 2:
        Stopping_criteria = 0.06
        alpha_amor = 0.18
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')    
        
    elif run_no == 3:
        Stopping_criteria = 0.08
        alpha_amor = 0.2
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')    
    elif run_no == 4:
        Stopping_criteria = 0.08#0.08
        alpha_amor = 0.2
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')                                              
    elif run_no == 5:
        Stopping_criteria = 0.12#0.1
        alpha_amor = 0.08
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')        
    elif run_no == 6:
        Stopping_criteria = 0.2#0.13
        alpha_amor = 0.08
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py') 
        
    elif run_no == 7:
        Stopping_criteria = 0.4#0.13
        alpha_amor = 0.1
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')     
    
    elif run_no == 8:
        Stopping_criteria = 0.4#0.13
        alpha_amor = 0.1
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')    
        
    elif run_no == 8:
        Stopping_criteria = 0.13
        alpha_amor = 0.14
        Delta_implant_amorphous = 30e-9
        runcell(0, '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_functions.py')        
        
    else:
        Stopping_criteria = 0.5
        alpha_amor = 0.6

    
    runcell('Reading position files of Carbon and Silicon', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')

    runcell('Reading geometry files', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')

    runcell('Calculation of erosion and deposition fluxes for Carbon and Silicon for each GITR particle', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')

    runcell('The following arrays will keep track of entries to be made in the next GITR run.', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')
    
    #if run_no == 8:
    #    runcell('making particle sources', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')

    runcell('Estimating the total time evolution for the surface model', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')

    runcell('Appending time to all the surface characteristics', '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/surface_model_subsequent_runs.py')
    
print("--------------------")
end = time.time()
print(end - start)          

#print("Characteristic time in seconds")    
#surface_in_question = 80
#1/np.polyfit(np.log(Surface_time), Concentration[20][surface_in_question,:], 1)    

#%%
surface_in_question = 34

import matplotlib.pyplot as plt    

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_evolution_C_Si.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")
# Record concentrations of all surface elements and their initial Z
Flux_proportionality = {}
for z in Zs:
    Concentration[z] = SurfaceConcentrationData['surface_concentration_{}'.format(z)][:,:]
    Flux_proportionality[z] = SurfaceConcentrationData['Flux_Conversion_{}'.format(z)][:]

Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)


plt.figure(figsize=(20, 20))
plt.plot(Surface_time,Concentration[6][surface_in_question,:],marker='^',label=r"$C_C$",linewidth=3)
plt.plot(Surface_time,Concentration[14][surface_in_question,:],marker='*',label=r"$C_{Si}$",linewidth=3)
plt.plot(Surface_time,Concentration[20][surface_in_question,:],marker='+',label=r"$C_{SiC}$",linewidth=3)
plt.plot(Surface_time,Concentration[26][surface_in_question,:],marker='+',label="C Enrichment",linewidth=3)
plt.plot(Surface_time,Concentration[34][surface_in_question,:],marker='+',label="Si Enrichment",linewidth=3)
plt.legend(fontsize=25)
plt.ylabel("Concentrations",fontsize=30,color='blue')
plt.xlabel("Time (s)",fontsize=30,color='blue')
#plt.xlim(0,7)
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2
plt.title("Representative surface mesh element: "+str(surface_in_question),fontsize=30)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.grid(color='gray', linestyle='-.',alpha=0.4)
#ax1.tick_params('both', length=20)

plt.show()
#%% Reading position files of Carbon and Silicon
#run_no = 2
#Stopping_criteria = 0.8

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/output_C_'+str(run_no)+'/positions.nc'
#dict_keys(['x', 'y', 'z', 'vx', 'vy', 'vz', 'transitTime', 'hitWall', 'surfaceHit', 'weight', 'charge', 'hasLeaked', 'distTraveled', 'time', 'dt'])
PositionData = netCDF4.Dataset(FileNameHistory)

surfacehit_C = np.array(PositionData['surfaceHit'])
surface_vx_C = np.array(PositionData['vx'])
surface_vy_C = np.array(PositionData['vy'])
surface_vz_C = np.array(PositionData['vz'])

Energy_particles_C = np.array(0.5*amu_C*1.66e-27*(surface_vx_C**2 + surface_vy_C**2 + surface_vz_C**2)/1.602e-19) # make sure that this energy is in eV
Angles_particles_C = np.arctan(surface_vx_C/surface_vz_C)*(180/np.pi)   # in degrees

count_C = 0

for i in surfacehit_C:
    if i != -1:
        count_C+=1
    
print(count_C," C have hit a mesh element (not necessarily a surface)")
# Reading position files of Silicon

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/output_Si_'+str(run_no)+'/positions.nc'
PositionData = Dataset(FileNameHistory, "r", format="NETCDF4")

surfacehit_Si = np.array(PositionData['surfaceHit'])
surface_vx_Si = np.array(PositionData['vx'])
surface_vy_Si = np.array(PositionData['vy'])
surface_vz_Si = np.array(PositionData['vz'])


Energy_particles_Si = np.array(0.5*amu_Si*1.66e-27*(surface_vx_Si**2 + surface_vy_Si**2 + surface_vz_Si**2)/1.602e-19) # make sure that this energy is in eV
Angles_particles_Si = np.arctan(surface_vx_Si/surface_vz_Si)*(180/np.pi)    #  in degrees

count_Si = 0
for i in surfacehit_Si:
    if i != -1:
        count_Si+=1
print(count_Si,"Si have hit a mesh element (not necessarily a surface)")



#%% Reading geometry files

GeomFile = "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/gitrGeom.cfg"
x1,x2,x3,y1,y2,y3,z1,z2,z3,area,surf,Atomic_no,a,b,c,d,in_direction,plane_norm = getGeom(GeomFile)
#x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surf,indir,Atomic_no = loadCFG(geomFile=GeomFile)

# Initialize the surface_evolution netcdf file
# Only care about surfaces


Zs = []

Surfaces = []

Surfaces_C = []
Surfaces_SiC = []

idx = np.arange(0,len(surf))
for surface,z,i in zip(surf,Atomic_no,idx):
    if surface!=0:
        Zs.append(z)
        Surfaces.append(i)
        if z==6:
            Surfaces_C.append(i)
        elif z==20:
            Surfaces_SiC.append(i)    

Zs = np.unique(Zs)
            
Zs = np.append(Zs,6)  # Adding Carbon
Zs = np.append(Zs,14)  # Adding Silicon
Zs = np.append(Zs,26)  # Adding enriched Carbon (C+SiC)
Zs = np.append(Zs,34)  # Adding enriched Silicon (Si+SiC)

print(Zs,"make up the", len(Surfaces),"surface mesh elements")
# print(Surfaces)


#Reading the surface features from the surface evolution netcdf file

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_evolution_C_Si.nc'

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

#%%  Calculation of erosion and deposition fluxes for Carbon and Silicon for each GITR particle

#import time

#start = time.time()

Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CSiC_Si_Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CSiC_C_Gamma_C_redep = np.zeros((len(Surfaces),1))
 
Y_CSi_Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CC_Gamma_C_redep = np.zeros((len(Surfaces),1))

Y_CSien_Gamma_C_redep = np.zeros((len(Surfaces),1))
Y_CCen_Gamma_C_redep = np.zeros((len(Surfaces),1))

Target_1 = 'SiC,C-crystalline' #averaged over three crystal orientations
Target_2 = 'SiC,Si-crystalline'
#Target_1 = 'SiC,C,Crystalline'
#Target_2 = 'SiC,Si,Crystalline'

Target_1_amor = 'SiC,C,Amorphous'             
Target_2_amor = 'SiC,Si,Amorphous' 

#Target_Si = 'Si-ENRICHEDSiC'
Target_Si = 'Si'

#Target_1 = 'SiC'
#Target_2 = 'SiC'

for i in range(len(Energy_particles_C)):
    
    if surfacehit_C[i] != -1 and int(surfacehit_C[i]) in Surfaces:
        
        
        surface_index = int(surfacehit_C[i])
        sr_object = Sputtering_and_reflection()
        
        Flux_C_local = Flux_proportionality[6][-1]/(Delta_t_gitr*area[surface_index])  # we start with initial weight 1 (uniform)
                
        Gamma_C_redep[surface_index] = Gamma_C_redep[surface_index] + Flux_C_local  # check this
        
        Y_CSiC_Si_Gamma_C_redep[surface_index] = Y_CSiC_Si_Gamma_C_redep[surface_index] + alpha_amor*sr_object.Calculate_PhysicalSputteringParameters('C',Target_2_amor,Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local + (1-alpha_amor)*sr_object.Calculate_PhysicalSputteringParameters('C',Target_2,Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
        Y_CSiC_C_Gamma_C_redep[surface_index] = Y_CSiC_C_Gamma_C_redep[surface_index] + alpha_amor*sr_object.Calculate_PhysicalSputteringParameters('C',Target_1_amor,Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local + (1-alpha_amor)*sr_object.Calculate_PhysicalSputteringParameters('C',Target_1,Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
                
        Y_CC_Gamma_C_redep[surface_index] = Y_CC_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
        Y_CSi_Gamma_C_redep[surface_index] = Y_CSi_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','Si',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
                
        Y_CSien_Gamma_C_redep[surface_index] = Y_CSien_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C',Target_Si,Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local
        Y_CCen_Gamma_C_redep[surface_index] = Y_CCen_Gamma_C_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('C','C',Energy_particles_C[i],Angles_particles_C[i])*Flux_C_local

print("done")

Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_SiSiC_Si_Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_SiSiC_C_Gamma_Si_redep = np.zeros((len(Surfaces),1))

Y_SiSi_Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_Si_C_Gamma_Si_redep = np.zeros((len(Surfaces),1))    # Y_Si_to_C

Y_SiSien_Gamma_Si_redep = np.zeros((len(Surfaces),1))
Y_SiCen_Gamma_Si_redep = np.zeros((len(Surfaces),1))    # Y_Si_to_C
        
        
for i in range(len(Energy_particles_Si)):
    if surfacehit_Si[i] != -1 and int(surfacehit_Si[i]) in Surfaces:
        surface_index = int(surfacehit_Si[i])

        sr_object = Sputtering_and_reflection()                        
                
        Flux_Si_local = Flux_proportionality[14][-1]/(Delta_t_gitr*area[surface_index]) # we start with initial weight 1 (uniform)
                
        Gamma_Si_redep[surface_index] = Gamma_Si_redep[surface_index] + Flux_Si_local  # check this
        Y_SiSiC_Si_Gamma_Si_redep[surface_index] = Y_SiSiC_Si_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','SiC',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
        Y_SiSiC_C_Gamma_Si_redep[surface_index] = Y_SiSiC_C_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','SiC',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
                
        Y_SiSi_Gamma_Si_redep[surface_index] = Y_SiSi_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','Si',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local                
        Y_Si_C_Gamma_Si_redep[surface_index] = Y_Si_C_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','C',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
                
        # Y_SiSien_Gamma_Si_redep[surface_index] = Y_SiSien_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','Si-ENRICHEDSiC',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local                
        # Y_SiCen_Gamma_Si_redep[surface_index] = Y_SiCen_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','C-ENRICHEDSiC',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local
        Y_SiSien_Gamma_Si_redep[surface_index] = Y_SiSien_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','Si',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local                
        Y_SiCen_Gamma_Si_redep[surface_index] = Y_SiCen_Gamma_Si_redep[surface_index] + sr_object.Calculate_PhysicalSputteringParameters('Si','C',Energy_particles_Si[i],Angles_particles_Si[i])*Flux_Si_local


print("done")

beta_eroSi1 = YHtoSi_Flux_H_in
beta_eroSi2 = YCtoSi_Flux_C_in 
beta_eroSi3 = Y_CSi_Gamma_C_redep
beta_eroSi4 = Y_SiSi_Gamma_Si_redep

beta_eroC1 = YHtoC_Flux_H_in
beta_eroC2 = YCtoC_Flux_C_in
beta_eroC3 = Y_CC_Gamma_C_redep
beta_eroC4 = Y_Si_C_Gamma_Si_redep

beta_SiC1_Si = YHtoSiC_Si_Flux_H_in + YCtoSiC_Si_Flux_C_in
beta_SiC1_C = YHtoSiC_C_Flux_H_in + YCtoSiC_C_Flux_C_in

beta_SiC_Si = YHtoSiC_Si_Flux_H_in + YCtoSiC_Si_Flux_C_in + Y_CSiC_Si_Gamma_C_redep + Y_SiSiC_Si_Gamma_Si_redep
beta_SiC_C = YHtoSiC_C_Flux_H_in + YCtoSiC_C_Flux_C_in + Y_CSiC_C_Gamma_C_redep + Y_SiSiC_C_Gamma_Si_redep + 0.1*XDtoCchem_Flux_D_in
   
Gamma_SiC_ero_global_Si = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))*beta_SiC_Si
Gamma_SiC_ero_global_C = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))*beta_SiC_C

#print(sum(Gamma_SiC_ero_global))

Gamma_enrichedSi_ero = np.reshape(Concentration[34][:,-1],(len(Surfaces),1))*(YDtoSien_Flux_D_in + YCtoSien_Flux_C_in + Y_CSien_Gamma_C_redep + Y_SiSien_Gamma_Si_redep)

Gamma_Si_ero_exclusive = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*(beta_eroSi1 + beta_eroSi2 + beta_eroSi3 + beta_eroSi4)
Gamma_Si_ero_global = Gamma_SiC_ero_global_Si + Gamma_Si_ero_exclusive + Gamma_enrichedSi_ero


Gamma_enrichedC_ero = np.reshape(Concentration[26][:,-1],(len(Surfaces),1))*(YDtoCen_Flux_D_in + YCtoCen_Flux_C_in + Y_CCen_Gamma_C_redep + Y_SiCen_Gamma_Si_redep +XDtoCchem_Flux_D_in)

Gamma_C_ero_exclusive = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))*(beta_eroC1 + beta_eroC2 + beta_eroC3 + beta_eroC3 + XDtoCchem_Flux_D_in)
Gamma_C_ero_global = Gamma_SiC_ero_global_C + Gamma_C_ero_exclusive + Gamma_enrichedC_ero

Gamma_C_dep_global = Gamma_C_redep + np.reshape(Concentration[14][:,-1],(len(Surfaces),1))*Flux_C_Background + np.reshape(Concentration[6][:,-1],(len(Surfaces),1))*beta_depC1 +  np.reshape(Concentration[20][:,-1],(len(Surfaces),1))*beta_depC2 + np.reshape(Concentration[26][:,-1],(len(Surfaces),1))*beta_depC1+np.reshape(Concentration[34][:,-1],(len(Surfaces),1))*Flux_C_Background 
Gamma_Si_dep_global = Gamma_Si_redep

Gamma_C_ero_global[np.argwhere(np.isnan(Gamma_C_ero_global))] = Gamma_C_ero_global[np.argwhere(np.isnan(Gamma_C_ero_global))+1]

print(np.argwhere(np.isnan(Gamma_C_ero_global)))

#end = time.time()
#print(end - start)          


#%% The following arrays will keep track of entries to be made in the next GITR run. 

nP_C_global = 0 #tracks total number of eroded particles
nP_Si_global = 0 #tracks total number of eroded particles

prop_Si = 0
prop_C = 0
prop_SiC = 1
prop_Sien = 1
prop_Cen = 1

for surface_index in range(len(Surfaces)):
    prop_Si = prop_Si + Gamma_Si_ero_global[surface_index]*area[surface_index]*Delta_t_gitr
    prop_C = prop_C + Gamma_C_ero_global[surface_index]*area[surface_index]*Delta_t_gitr

prop_Si = prop_Si/N_GITR
prop_C = prop_C/N_GITR

particleSourceDict_C = {}
for i,surface in enumerate(Surfaces):

    num_particles = np.array(Gamma_C_ero_global[i]*Delta_t_gitr*area[surface]/prop_C).item()
    if num_particles!=0: 
        print("Surface:",surface,"C particles:",int(num_particles))
        particleSourceDict_C[surface] = round(num_particles)
        nP_C_global+=round(num_particles)
        
print("C particles generated")        
particleSourceDict_Si = {}
for i,surface in enumerate(Surfaces):
    num_particles = np.array(Gamma_Si_ero_global[i]*Delta_t_gitr*area[surface]/prop_Si).item()
    
    if num_particles!=0: 
        print("Surface:",surface,"Si particles:",int(num_particles))
        particleSourceDict_Si[surface] = round(num_particles)
        nP_Si_global+=round(num_particles)            
print("Si particles generated")  
#%%  making particle sources

import time

start = time.time()

makeParticleSource('C',particleSourceDict_C, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/particleConf_C.nc")
if nP_Si_global>0:     
    makeParticleSource('Si',particleSourceDict_Si, "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/gitrGeom.cfg", "/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/particleConf_Si.nc")
            
end = time.time()
print(end - start)          

#%%  Estimating the total time evolution for the surface model

last_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))
last_entry_Si = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))
last_entry_SiC = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))
last_entry_Cenriched = np.reshape(Concentration[26][:,-1],(len(Surfaces),1))
last_entry_Sienriched = np.reshape(Concentration[34][:,-1],(len(Surfaces),1))


Gamma_C_net = Gamma_C_dep_global - Gamma_C_ero_exclusive
Gamma_Si_net = Gamma_Si_dep_global - Gamma_Si_ero_exclusive
Gamma_SiC_net = -np.maximum(Gamma_SiC_ero_global_Si,Gamma_SiC_ero_global_C)
Gamma_Cenr_net = np.zeros((len(Surfaces),1))
Gamma_Sienr_net = np.zeros((len(Surfaces),1))

Gamma_C_bulk = np.zeros((len(Surfaces),1))
Gamma_Si_bulk = np.zeros((len(Surfaces),1))
Gamma_SiC_bulk = np.zeros((len(Surfaces),1))
Gamma_Sienr_bulk = np.zeros((len(Surfaces),1))
Gamma_Cenr_bulk = np.zeros((len(Surfaces),1))



#print(Gamma_C_net)

for surface_index in range(len(Surfaces)):
    
    enrichment_flux = abs(Gamma_SiC_ero_global_C[surface_index]-Gamma_SiC_ero_global_Si[surface_index])
    if (Gamma_SiC_ero_global_Si[surface_index]>Gamma_SiC_ero_global_C[surface_index]):
        print("ERROR")
        Gamma_Cenr_net[surface_index] = enrichment_flux#-Gamma_enrichedC_ero[surface_index]
        Gamma_Sienr_net[surface_index] = 0.0#-Gamma_enrichedSi_ero[surface_index]
    else:  #Gamma_SiC_ero_global_Si<Gamma_SiC_ero_global_C
        #Gamma_Sienr_net[surface_index] = enrichment_flux-Gamma_enrichedSi_ero[surface_index]
        Gamma_Sienr_net[surface_index] = enrichment_flux#-Gamma_enrichedSi_ero[surface_index]
        Gamma_Cenr_net[surface_index] = 0.0#-Gamma_enrichedC_ero[surface_index]
            
        
    Total_net_flux = Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index]+Gamma_Cenr_net[surface_index] + Gamma_Sienr_net[surface_index]
    
    if Total_net_flux > 0: # deposition regime
        #print("deposition")
        
        
        Gamma_C_bulk[surface_index] = last_entry_C[surface_index,0]*Total_net_flux
        Gamma_Si_bulk[surface_index] = last_entry_Si[surface_index,0]*Total_net_flux
        Gamma_SiC_bulk[surface_index] = last_entry_SiC[surface_index,0]*Total_net_flux
        
        Gamma_Cenr_bulk[surface_index] = last_entry_Cenriched[surface_index,0]*Total_net_flux
        Gamma_Sienr_bulk[surface_index] = last_entry_Sienriched[surface_index,0]*Total_net_flux
    
    else:  #  erosion regime
        #print("erosion")
        Gamma_C_bulk[surface_index] = 0.0#(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
        Gamma_Si_bulk[surface_index] = 0.0
        Gamma_SiC_bulk[surface_index] = Total_net_flux
        
        Gamma_Cenr_bulk[surface_index] = 0.0
        Gamma_Sienr_bulk[surface_index] = 0.0
    

RHS_C = np.abs(Gamma_C_net - Gamma_C_bulk)
RHS_Si = np.abs(Gamma_Si_net - Gamma_Si_bulk)
RHS_SiC = np.abs(Gamma_SiC_net - Gamma_SiC_bulk)
RHS_Cenr = np.abs(Gamma_Cenr_net - Gamma_Cenr_bulk)
RHS_Sienr = np.abs(Gamma_Sienr_net - Gamma_Sienr_bulk)

#RHS_C + RHS_Si + RHS_SiC + RHS_Cenr + RHS_Sienr 
# Gamma_C_bulk + Gamma_Si_bulk + Gamma_SiC_bulk + Gamma_Cenr_bulk + Gamma_Sienr_bulk

# Stopping_criteria = 0.2 #changed from 0.1
       
Delta_t_surface_estimate_C = (Delta_implant_amorphous*n_atom*Stopping_criteria)/RHS_C
Delta_t_surface_estimate_Si = (Delta_implant_amorphous*n_atom*Stopping_criteria)/RHS_Si
Delta_t_surface_estimate_SiC = (Delta_implant_amorphous*n_atom*Stopping_criteria)/RHS_SiC

# Delta_t_surface = min(np.amin(Delta_t_surface_estimate_C),np.amin(Delta_t_surface_estimate_Si),np.amin(Delta_t_surface_estimate_SiC))        
Delta_t_surface = min(Delta_t_surface_estimate_SiC)      

Time = Delta_t_surface
Time_steps = 1e2
Delta_Time = Delta_t_surface/Time_steps #Delta_t/Time_steps   This is the time step variable
Delta_t_Stopping = 0
#print("Till here")

new_entry_C = np.reshape(Concentration[6][:,-1],(len(Surfaces),1))   # populating it with last entry
new_entry_Si = np.reshape(Concentration[14][:,-1],(len(Surfaces),1))
new_entry_SiC = np.reshape(Concentration[20][:,-1],(len(Surfaces),1))
new_entry_Cenr = np.reshape(Concentration[26][:,-1],(len(Surfaces),1))
new_entry_Sienr = np.reshape(Concentration[34][:,-1],(len(Surfaces),1))
        
for t in range(1,int(Time_steps)):
    
    print(t," out of ",Time_steps)
    
    Gamma_SiC_ero_global_Si = new_entry_SiC*beta_SiC_Si
    Gamma_SiC_ero_global_C = new_entry_SiC*beta_SiC_C

    #print(sum(Gamma_SiC_ero_global))

    Gamma_enrichedSi_ero = new_entry_Sienr*(YDtoSien_Flux_D_in + YCtoSien_Flux_C_in + Y_CSien_Gamma_C_redep + Y_SiSien_Gamma_Si_redep)

    Gamma_Si_ero_exclusive = new_entry_Si*(beta_eroSi1 + beta_eroSi2 + beta_eroSi3 + beta_eroSi4)
    Gamma_Si_ero_global = Gamma_SiC_ero_global_Si + Gamma_Si_ero_exclusive + Gamma_enrichedSi_ero


    Gamma_enrichedC_ero = new_entry_Cenr*(YDtoCen_Flux_D_in + YCtoCen_Flux_C_in + Y_CCen_Gamma_C_redep + Y_SiCen_Gamma_Si_redep + XDtoCchem_Flux_D_in)

    Gamma_C_ero_exclusive = new_entry_C*(beta_eroC1 + beta_eroC2 + beta_eroC3 + beta_eroC3 + XDtoCchem_Flux_D_in)
    Gamma_C_ero_global = Gamma_SiC_ero_global_C + Gamma_C_ero_exclusive + Gamma_enrichedC_ero

    Gamma_C_dep_global = Gamma_C_redep + new_entry_Si*Flux_C_Background + new_entry_C*beta_depC1 +  new_entry_SiC*beta_depC2 + new_entry_Cenr*beta_depC1 +new_entry_Sienr*Flux_C_Background   
    Gamma_Si_dep_global = Gamma_Si_redep
    
    
    # determining erosion or deposition
    Gamma_C_net = Gamma_C_dep_global - Gamma_C_ero_exclusive
    Gamma_Si_net = Gamma_Si_dep_global - Gamma_Si_ero_exclusive
    Gamma_SiC_net = -np.maximum(Gamma_SiC_ero_global_Si,Gamma_SiC_ero_global_C)
    Gamma_Cenr_net = np.zeros((len(Surfaces),1))
    Gamma_Sienr_net = np.zeros((len(Surfaces),1))

    Gamma_C_bulk = np.zeros((len(Surfaces),1))
    Gamma_Si_bulk = np.zeros((len(Surfaces),1))
    Gamma_SiC_bulk = np.zeros((len(Surfaces),1))
    Gamma_Sienr_bulk = np.zeros((len(Surfaces),1))
    Gamma_Cenr_bulk = np.zeros((len(Surfaces),1))

    
    for surface_index in range(len(Surfaces)):

        enrichment_flux = abs(Gamma_SiC_ero_global_Si[surface_index]-Gamma_SiC_ero_global_C[surface_index])
        if (Gamma_SiC_ero_global_Si[surface_index]>Gamma_SiC_ero_global_C[surface_index]):
            # Gamma_Cenr_net[surface_index] = enrichment_flux-Gamma_enrichedC_ero[surface_index]
            # Gamma_Sienr_net[surface_index] = -Gamma_enrichedSi_ero[surface_index]
            Gamma_Cenr_net[surface_index] = enrichment_flux#-Gamma_enrichedC_ero[surface_index]
            Gamma_Sienr_net[surface_index] = 0.0#-Gamma_enrichedSi_ero[surface_index]
        else:
            Gamma_Sienr_net[surface_index] = enrichment_flux#-Gamma_enrichedSi_ero[surface_index]
            Gamma_Cenr_net[surface_index] = 0.0#-Gamma_enrichedC_ero[surface_index]
        
        Total_net_flux = Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index]+Gamma_Cenr_net[surface_index] + Gamma_Sienr_net[surface_index]    

        if Total_net_flux > 0: # deposition regime
            #print("deposition")
            Gamma_C_bulk[surface_index] = last_entry_C[surface_index,0]*Total_net_flux
            Gamma_Si_bulk[surface_index] = last_entry_Si[surface_index,0]*Total_net_flux
            Gamma_SiC_bulk[surface_index] = last_entry_SiC[surface_index,0]*Total_net_flux
            
            Gamma_Cenr_bulk[surface_index] = last_entry_Cenriched[surface_index,0]*Total_net_flux
            Gamma_Sienr_bulk[surface_index] = last_entry_Sienriched[surface_index,0]*Total_net_flux
        
        else:  #  erosion regime
            #print("erosion")
            Gamma_C_bulk[surface_index] = 0.0#(Gamma_C_net[surface_index]+Gamma_Si_net[surface_index]+Gamma_SiC_net[surface_index])
            Gamma_Si_bulk[surface_index] = 0.0
            Gamma_SiC_bulk[surface_index] = Total_net_flux
            
            Gamma_Cenr_bulk[surface_index] = 0.0
            Gamma_Sienr_bulk[surface_index] = 0.0
               
    
    #print(t)
    new_entry_C = new_entry_C + Delta_Time*(Gamma_C_net - Gamma_C_bulk)/(Delta_implant_amorphous*n_atom)    
    new_entry_Si = new_entry_Si + Delta_Time*(Gamma_Si_net - Gamma_Si_bulk)/(Delta_implant_amorphous*n_atom)    
    new_entry_SiC = new_entry_SiC + Delta_Time*(Gamma_SiC_net - Gamma_SiC_bulk)/(Delta_implant_amorphous*n_atom)
    new_entry_Cenr = new_entry_Cenr + Delta_Time*(Gamma_Cenr_net - Gamma_Cenr_bulk)/(Delta_implant_amorphous*n_atom)
    new_entry_Sienr = new_entry_Sienr + Delta_Time*(Gamma_Sienr_net - Gamma_Sienr_bulk)/(Delta_implant_amorphous*n_atom)
    
    Delta_t_Stopping += Delta_Time
        
    # if (np.abs(new_entry_C-last_entry_C)>Stopping_criteria).any() or (np.abs(new_entry_Si-last_entry_Si)>Stopping_criteria).any() or (np.abs(new_entry_SiC-last_entry_SiC)>Stopping_criteria).any():
    #     print(Delta_t_Stopping," Delta_t_Stopping ", t)
    #     break
    if (np.abs(new_entry_SiC-last_entry_SiC)>Stopping_criteria).any():
        print(Delta_t_Stopping," Delta_t_Stopping ", t)
        break    
    
# last_entry_C + last_entry_Si + last_entry_SiC+last_entry_Cenriched + last_entry_Sienriched
# new_entry_C + new_entry_Si + new_entry_SiC + new_entry_Cenr + new_entry_Sienr
#%% Appending time to all the surface characteristics

Concentration[6] = np.concatenate((Concentration[6],new_entry_C),axis=1)
Concentration[14] = np.concatenate((Concentration[14],new_entry_Si),axis=1)
Concentration[20] = np.concatenate((Concentration[20],new_entry_SiC),axis=1)
Concentration[26] = np.concatenate((Concentration[26],new_entry_Cenr),axis=1)
Concentration[34] = np.concatenate((Concentration[34],new_entry_Sienr),axis=1)


Flux_proportionality[6] = np.append(Flux_proportionality[6],prop_C)
Flux_proportionality[14] = np.append(Flux_proportionality[14],prop_Si)
Flux_proportionality[20] = np.append(Flux_proportionality[20],prop_SiC)
Flux_proportionality[26] = np.append(Flux_proportionality[26],prop_Cen)
Flux_proportionality[34] = np.append(Flux_proportionality[34],prop_Sien)

Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t_Stopping)

#Writing the surface features with time
os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_evolution_C_Si.nc")
ncFile = netCDF4.Dataset('/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_evolution_C_Si.nc', 'w', format='NETCDF4')

s_number_dim = ncFile.createDimension('surface_dim', len(Surfaces)) # surface number dimension
s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
s_time = ncFile.createVariable('time', np.float32, ('time_dim',))


s_concentration = {}
flux_proportionality = {}
for z in Concentration.keys():
    s_concentration[z] = ncFile.createVariable('surface_concentration_{}'.format(z), np.float64, ('surface_dim','time_dim'))
    flux_proportionality[z] = ncFile.createVariable('Flux_Conversion_{}'.format(z),np.float64,('time_dim'))
    print(z)


s_number[:] = np.linspace(1,len(Surfaces),len(Surfaces))
s_time[:] = Surface_time


for z in Zs:
   s_concentration[z][:,:] = Concentration[z]
   flux_proportionality[z][:] = Flux_proportionality[z]

ncFile.close()


#%% PLOTTING FROM HERE ON

import matplotlib.pyplot as plt    

FileNameSurfaceConcentration='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/surface_evolution_C_Si.nc'

SurfaceConcentrationData = Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")

# Record concentrations of all surface elements and their initial Z
Flux_proportionality = {}
for z in Zs:
    Concentration[z] = SurfaceConcentrationData['surface_concentration_{}'.format(z)][:,:]
    Flux_proportionality[z] = SurfaceConcentrationData['Flux_Conversion_{}'.format(z)][:]

Surface_time = SurfaceConcentrationData['time'][:]
Surface_number = SurfaceConcentrationData['surface_number'][:]
counter = len(Surface_time)

surface_in_question = 83

plt.figure()
plt.plot(Surface_time,Concentration[6][surface_in_question,:],marker='^',label='C_C')
plt.plot(Surface_time,Concentration[14][surface_in_question,:],marker='*',label='C_Si')
plt.plot(Surface_time,Concentration[20][surface_in_question,:],marker='+',label='C_SiC')

plt.legend()
plt.title("Surface Element %d" % surface_in_question)
plt.xlabel("Time (s)")
plt.ylabel("Concentrations")
plt.show()

#%% plotting

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

 
cubic_model1 = interp1d(Surface_time, Concentration[6][surface_index_C,:], kind = "cubic")
Y_C = cubic_model1(Surface_time)

cubic_model2 = interp1d(Surface_time, Concentration[14][surface_index_C,:], kind = "cubic")
Y_Si = cubic_model2(Surface_time)

cubic_model1 = interp1d(Surface_time, Concentration[6][surface_index_C,:], kind = "cubic")
Y_SiC = cubic_model1(Surface_time)
 
# Plotting the Graph
X_=np.linspace(x.min(), x.max(), 500)
Y_=cubic_interploation_model(X_)



#%%
import matplotlib.pyplot as plt    

fig, ax = plt.subplots(figsize=(20, 20))

ax.plot(Surface_time,Concentration[6][surface_in_question,:],marker='^',label=r"$C_C$",linewidth=3)
ax.plot(Surface_time,Concentration[14][surface_in_question,:],marker='*',label=r"$C_{Si}$",linewidth=3)
ax.plot(Surface_time,Concentration[20][surface_in_question,:],marker='+',label=r"$C_{SiC}$",linewidth=3)
ax.plot(Surface_time,Concentration[26][surface_in_question,:],marker='+',label="C Enrichment",linewidth=3)
ax.plot(Surface_time,Concentration[34][surface_in_question,:],marker='+',label="Si Enrichment",linewidth=3)
ax.legend(fontsize=25)
ax.set_ylabel("Concentrations",fontsize=30,color='black')
ax.set_xlabel("Time (s)",fontsize=30,color='black')
plt.xlim(0,5)
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2
#ax.set_title("Representative surface mesh element: "+str(surface_in_question),fontsize=30)
ax.tick_params(axis='both', labelsize=25,length=10)
ax.grid(color='gray', linestyle='-.',alpha=0.4)
#ax1.tick_params('both', length=20)

plt.show()
