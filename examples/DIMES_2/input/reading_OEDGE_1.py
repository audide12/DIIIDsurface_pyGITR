#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:17:00 2023

@author: de
"""

import pandas as pd
import netCDF4
import numpy as np

df = pd.read_table("/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/d-176492-bg-a5.plasma.dat", sep="\s+")


R  = np.array(df['R(m)'])
Z  = np.array(df['Z(m)'])
Ne = np.array(df['ne(m-3)'])
Te = np.array(df['Te(eV)'])
Ti = np.array(df['Ti(eV)'])
V_para = np.array(df['Vpara(m/s)'])
E_para = np.array(df['Epara(V/m)'])
B_tot = np.array(df['Btot(T)'])
B_R = np.array(df['B_R'])
B_Z = np.array(df['B_Z'])
B_T = np.array(df['B_T'])


B_R = B_R*B_tot
B_Z = B_Z*B_tot
B_T = B_T*B_tot
# Converting 1D arrays to 2D arrays
nR = 251
nZ = 251

Z_1D = np.zeros(nZ)
R_1D = np.zeros(nR)
Ne_2D = np.zeros((nZ,nR))
Te_2D = np.zeros((nZ,nR))
Ti_2D = np.zeros((nZ,nR))
V_para_2D = np.zeros((nZ,nR))
E_para_2D = np.zeros((nZ,nR))
B_tot_2D = np.zeros((nZ,nR))
B_R_2D = np.zeros((nZ,nR))
B_Z_2D = np.zeros((nZ,nR))
B_T_2D = np.zeros((nZ,nR))


count = 0
for r in range(nR):
    R_1D[r] = R[count]
    count = count+nZ
    
count = 0
for z in range(nZ):
    Z_1D[z] = Z[count]
    count = count+1    

count = 0
for r in range(nR):
    for z in range(nZ):
        Ne_2D[z,r] = Ne[count]
        Te_2D[z,r] = Te[count]
        Ti_2D[z,r] = Ti[count]
        V_para_2D[z,r] = V_para[count]
        E_para_2D[z,r] = E_para[count]
        B_tot_2D[z,r] = B_tot[count]
        B_R_2D[z,r] = B_R[count]
        B_Z_2D[z,r] = B_Z[count]
        B_T_2D[z,r] = B_T[count]
        count = count+1

# Adjusting R and Z

r_shift = 1.485
z_shift = -1.250
R_1D = R_1D + r_shift
Z_1D = Z_1D + z_shift
    
        
#   Gradient Calculation

gradTe_R_2D = np.zeros((nZ,nR))
gradTi_R_2D = np.zeros((nZ,nR))
gradTe_Z_2D = np.zeros((nZ,nR))
gradTi_Z_2D = np.zeros((nZ,nR))

for r in range(nR-1):
    gradTe_R_2D[:,r] = (Te_2D[:,r+1] - Te_2D[:,r])/(R_1D[r+1]-R_1D[r])
    gradTi_R_2D[:,r] = (Ti_2D[:,r+1] - Ti_2D[:,r])/(R_1D[r+1]-R_1D[r])

for z in range(nZ-1):
    gradTe_Z_2D[z,:] = (Te_2D[z+1,:] - Te_2D[z,:])/(Z_1D[z+1]-Z_1D[z])
    gradTi_Z_2D[z,:] = (Ti_2D[z+1,:] - Ti_2D[z,:])/(Z_1D[z+1]-Z_1D[z])        
        


rootgrp = netCDF4.Dataset("/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/profiles_created.nc", "w", format="NETCDF4")

# Create dimensions
nR_1 = rootgrp.createDimension("nR", nR)
nZ_1 = rootgrp.createDimension("nZ", nZ)


r    = rootgrp.createVariable("r","f8",("nR"))
z    = rootgrp.createVariable("z","f8",("nR"))
te   = rootgrp.createVariable("te","f8",("nZ","nR"))
ti   = rootgrp.createVariable("ti","f8",("nZ","nR"))
ne   = rootgrp.createVariable("ne","f8",("nZ","nR"))
ni   = rootgrp.createVariable("ni","f8",("nZ","nR"))
Btot = rootgrp.createVariable("Btot","f8",("nZ","nR"))
br   = rootgrp.createVariable("br","f8",("nZ","nR"))
bt   = rootgrp.createVariable("bt","f8",("nZ","nR"))
bz   = rootgrp.createVariable("bz","f8",("nZ","nR"))
Velocity = rootgrp.createVariable("Velocity","f8",("nZ","nR"))
vr   = rootgrp.createVariable("vr","f8",("nZ","nR"))
vt   = rootgrp.createVariable("vt","f8",("nZ","nR"))
vz   = rootgrp.createVariable("vz","f8",("nZ","nR"))
efield = rootgrp.createVariable("efield","f8",("nZ","nR"))
Er   = rootgrp.createVariable("Er","f8",("nZ","nR"))
Et   = rootgrp.createVariable("Et","f8",("nZ","nR"))
Ez   = rootgrp.createVariable("Ez","f8",("nZ","nR"))
gradTeR = rootgrp.createVariable("gradTeR","f8",("nZ","nR"))
gradTeZ = rootgrp.createVariable("gradTeZ","f8",("nZ","nR"))

gradTeY = rootgrp.createVariable("gradTeY","f8",("nZ","nR"))   # added : dummy zero values

gradTiR = rootgrp.createVariable("gradTiR","f8",("nZ","nR"))
gradTiZ = rootgrp.createVariable("gradTiZ","f8",("nZ","nR"))

gradTiY = rootgrp.createVariable("gradTiY","f8",("nZ","nR"))   # added : dummy zero values

r[:] = R_1D
z[:] = Z_1D
te[:] = Te_2D
ti[:] = Ti_2D
ne[:] = Ne_2D
ni[:] = Ne_2D # Set equal to ne

Btot[:] = B_tot_2D
br[:] = B_R_2D
bt[:] = B_T_2D
bz[:] = B_Z_2D

gradTeR[:] = gradTe_R_2D
gradTeZ[:] = gradTe_Z_2D

gradTiR[:] = gradTi_R_2D
gradTiZ[:] = gradTi_Z_2D

np.zeros((nZ,nR))

gradTeY[:] = np.zeros((nZ,nR))
gradTiY[:] = np.zeros((nZ,nR))


V_r_2D = (V_para_2D*B_R_2D)/B_tot_2D
V_t_2D = (V_para_2D*B_T_2D)/B_tot_2D
V_z_2D = (V_para_2D*B_Z_2D)/B_tot_2D

Velocity[:] = V_para_2D
vr[:] = V_r_2D
vt[:] = V_t_2D
vz[:] = V_z_2D

efield_r_2D = (E_para_2D*B_R_2D)/B_tot_2D
efield_t_2D = (E_para_2D*B_T_2D)/B_tot_2D
efield_z_2D = (E_para_2D*B_Z_2D)/B_tot_2D
    
efield[:] = E_para_2D
Er[:] = efield_r_2D
Et[:] = efield_t_2D
Ez[:] = efield_z_2D

rootgrp.close()

#   Writing the magnetic field

datafile = netCDF4.Dataset("/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES_2/input/bField_created.nc", "w", format="NETCDF4")

# Create dimensions
nR_1 = datafile.createDimension("nR", nR)
nZ_1 = datafile.createDimension("nZ", nZ)

# Create variables
r = datafile.createVariable("r","f8",("nR"))
z = datafile.createVariable("z","f8",("nZ"))
br = datafile.createVariable("br","f8",("nZ","nR"))
bt = datafile.createVariable("bt","f8",("nZ","nR"))
bz = datafile.createVariable("bz","f8",("nZ","nR"))

# Initialize variables
r[:] = R_1D
z[:] = Z_1D
br[:] = B_R_2D
bt[:] = B_T_2D
bz[:] = B_Z_2D
    
datafile.close()

#%%  Plotting


import matplotlib.pyplot as plt
r_reference = 100
z_reference = 12    


plt.figure()       
plt.plot(R_1D,B_tot_2D[z_reference,:],label='Total B')
plt.plot(R_1D,B_R_2D[z_reference,:],label='Br')
plt.plot(R_1D,B_Z_2D[z_reference,:],label='Bz')
plt.plot(R_1D,B_T_2D[z_reference,:],label='Bt')
plt.legend()
plt.xlabel('R (m)')
plt.ylabel('B (T)')
plt.title('B vs R')

plt.figure()       
plt.plot(Z_1D,B_tot_2D[:,r_reference],label='Total B')
plt.plot(Z_1D,B_R_2D[:,r_reference],label='Br')
plt.plot(Z_1D,B_Z_2D[:,r_reference],label='Bz')
plt.plot(Z_1D,B_T_2D[:,r_reference],label='Bt')
plt.legend()
plt.xlabel('Z (m)')
plt.ylabel('B (T)')
plt.title('B vs Z')

plt.figure()        
plt.plot(Z_1D,Te_2D[:,r_reference],label='Te')
plt.plot(Z_1D,Ti_2D[:,r_reference],label='Ti')
plt.legend()
plt.xlabel('Z (m)')
plt.ylabel('T(eV)')
plt.title('T vs Z')

plt.figure()
plt.plot(R_1D,Te_2D[z_reference,:],label='Te')
plt.plot(R_1D,Ti_2D[z_reference,:],label='Ti')
plt.legend()
plt.xlabel('R (m)')
plt.ylabel('T(eV')
plt.title('T vs R')

plt.figure()
plt.plot(R_1D,Ne_2D[z_reference,:],label='Ne')
#plt.plot(R_1D,Ni_2D[z_reference,:],label='Ni')
plt.legend()
plt.xlabel('R (m)')
plt.ylabel(r'$Ne(m^{-3})$')
plt.title('Ne vs R')


#%% Plotting Background
r_reference = 100
z_reference = 0 

plt.figure()
plt.plot(R_1D,Te_2D[z_reference,:],label=r'$T_e$ (eV)',linewidth=5,color='black')
plt.plot(R_1D,Ne_2D[z_reference,:]/1e18,label=r'$n_{e} \times 10^{18} (m^{-3})$',linewidth=5,color='blue')
#plt.plot(R_1D,Ti_2D[z_reference,:],label='Ti',linewidth=5)
plt.legend()
plt.xlabel('R (m)',fontsize=20)
#plt.ylabel('T(eV)')
plt.yticks(fontsize=20)

plt.xticks([1.42,1.45,1.46,1.474,1.5,1.51,1.6])
plt.xticks(fontsize=20)
plt.axvline(x=1.474,linewidth=5)
plt.legend(fontsize=20,loc='upper right')

#plt.title('T vs R')
plt.xlim(1.42, 1.6) 
plt.show()

plt.figure()       
plt.plot(R_1D,B_tot_2D[z_reference,:],label='Total B')
plt.plot(R_1D,B_R_2D[z_reference,:],label='Br')
plt.plot(R_1D,B_Z_2D[z_reference,:],label='Bz')
plt.plot(R_1D,B_T_2D[z_reference,:],label='Bt')
plt.legend()
plt.xlabel('R (m)')
plt.ylabel('B (T)')
plt.title('B vs R')
#%% 
r_reference = 100
z_reference = 0 

plt.figure()       
plt.plot(R_1D,(Te_2D[z_reference,:])**2/(Ne_2D[z_reference,:]/1e18),linewidth=5,color='black')
#plt.plot(R_1D,Ne_2D[z_reference,:]/1e18,label=r'$n_{e} \times 10^{18} (m^{-3})$',linewidth=5,color='blue')
plt.axvline(x=1.474,linewidth=5)
plt.legend()
plt.xlabel('R (m)')
plt.ylabel('T^2/Ne')
plt.title('Hypothesis')