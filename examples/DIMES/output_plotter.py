#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:30:34 2022

@author: de
"""

# History plotting

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np


FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/output/history.nc'
HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
x = np.array(HistoryData.variables['x'])
z = np.array(HistoryData.variables['z'])
y = np.array(HistoryData.variables['y'])

nT = HistoryData.dimensions['nT'].size
nP = HistoryData.dimensions['nP'].size

fig = plt.figure() 
ax = fig.add_subplot(111,projection='3d') 
# ax.axes.set_xlim3d(left=1.45, right=1.55) 
# ax.axes.set_ylim3d(bottom=-0.04, top=0.1) 
# ax.axes.set_zlim3d(bottom=0.0, top=0.001) 

# for i in range(0,100):
#     ax.plot(x[i,:],y[i,:],z[i,:])

k=1012
ax.plot(x[k,:],y[k,:],z[k,:])    
    
#g.Plot_output(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)    
plt.show()
#%%

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf.nc'
Particle = Dataset(FileNameHistory, "r", format="NETCDF4")

vx_C = np.array(Particle['vx'])
vy_C = np.array(Particle['vy'])
vz_C = np.array(Particle['vz'])

Energy_particles = np.array(0.5*12*1.66e-27*(vx_C**2 + vz_C**2 + vy_C**2)/1.602e-19) # make sure that this energy is in eV
fig = plt.figure()
plt.hist(Energy_particles,bins = 100)
plt.show()


#%%
fig = plt.figure()
plt.hist(Energy_particles_C,bins = 20)
plt.show()


#%%


import matplotlib.pyplot as plt
import numpy as np

np.random.seed(42)
x = np.random.normal(size=1000)

plt.hist(x, density=True, bins=30)  # density=False would make counts
plt.ylabel('Probability')
plt.xlabel('Data')
plt.show()


#%%

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

FileNameHistory='/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/output_C_11/history.nc'
HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
x = np.array(HistoryData.variables['x'])
z = np.array(HistoryData.variables['z'])
y = np.array(HistoryData.variables['y'])

nT = HistoryData.dimensions['nT'].size
nP = HistoryData.dimensions['nP'].size

fig = plt.figure() 
ax = fig.add_subplot(111,projection='3d') 
# ax.axes.set_xlim3d(left=1.3, right=1.65) 
# ax.axes.set_ylim3d(bottom=-0.2, top=0.2) 
# ax.axes.set_zlim3d(bottom=0.0, top=0.3) 

# ax.axes.set_xlim3d(left=1.45, right=1.55) 
# ax.axes.set_ylim3d(bottom=-0.04, top=0.1) 
# ax.axes.set_zlim3d(bottom=0.0, top=0.001) 

for i in range(0,nP):
    ax.plot(x[i,:],y[i,:],z[i,:])

#g.Plot_output(ElemAttr='Z', Alpha=0.1, fig=fig, ax=ax)    
plt.show()


#%%

import matplotlib.pyplot as plt    

surface_in_question = 100
surface_index_C = surface_in_question
surface_index_Si = surface_in_question
surface_index_SiC = surface_in_question

plt.figure()
plt.plot(Surface_time,Concentration[6][surface_index_C,:],'k',label='C_C')
plt.plot(Surface_time,Concentration[14][surface_index_Si,:],'b',label='C_Si')
plt.plot(Surface_time,Concentration[20][surface_index_SiC,:],'g',label='C_SiC')

# plt.scatter(Surface_time,Concentration[6][surface_index_C,:],s=50,marker='^',label='C_C')
# plt.scatter(Surface_time,Concentration[14][surface_index_Si,:],s=50,marker='*',label='C_Si')
# plt.scatter(Surface_time,Concentration[20][surface_index_SiC,:],s=50,marker='+',label='C_SiC')

plt.legend()
plt.title("Surface Element %d" % surface_in_question)
plt.show()
