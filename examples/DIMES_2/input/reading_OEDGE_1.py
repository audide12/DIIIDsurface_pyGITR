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

gradTiR = rootgrp.createVariable("gradTiR","f8",("nZ","nR"))
gradTiZ = rootgrp.createVariable("gradTiZ","f8",("nZ","nR"))


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

gradTeR[:] = gradTe_R
gradTeZ[:] = gradTe_Z

gradTiR[:] = gradTi_R
gradTiZ[:] = gradTi_Z



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



