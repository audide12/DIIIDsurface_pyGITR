#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:33:44 2023

@author: de
"""

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np


Alpha_c_array = [0.01,0.0125,0.015,0.0175,0.02,0.025]

Delta_implant_array = [2,10,40,80,120]

tau_decay_array = [ [1.428,	0.8333,	0.8333,	0.33,	0.33,	0.33],
                    [7.14,	4,	3.333,	2,	2,	1],
                    [20,   20,	12.5,	10,	6.66666,	3.333],
                    [37.03,	33.333,	33.333,	20.0,	12.5,	8.333],
                    [58.823,	50,	50,	33.333,	25.0,	20]]


Gross_erosion_array = [ [4.10E+19,	4.81E+19,	5.52E+19,	6.28E+19,	6.93E+19,	8.35E+19],
                    [4.10E+19,	4.86E+19,	5.52E+19,	6.33E+19,	6.93E+19,	8.35E+19],
                    [4.10E+19,	4.81E+19,	5.46E+19,	6.28E+19,	6.93E+19,	8.35E+19],
                    [4.10E+19,	4.81E+19,	5.46E+19,	6.28E+19,	6.88E+19,	8.35E+19],
                    [4.15E+19,	4.86E+19,	5.57E+19,	6.28E+19,	6.93E+19,	8.35E+19]]

SiC_conc_array = [ [0.811,	0.810,	0.808,	0.805,	0.805,	0.799],
                    [0.811,	0.81,	0.806,	0.805,	0.802,	0.801],
                    [0.813,	0.81,	0.807,	0.802,	0.799,	0.802],
                    [0.813,	0.81,	0.807,	0.802,	0.799,	0.802],
                    [0.81,	0.81,	0.802,	0.802,	0.802,	0.803]]


C_conc_array = [[ 0.089,	0.091,	0.093,	0.09,	0.09,	0.09 ],
[0.091,	0.091,	0.091,	0.088,	0.09,	0.091],
[0.089,	0.089,	0.090,	0.09,	0.088,	0.091],
[0.087,	0.089,	0.089,	0.088,	0.085,	0.091],
[0.089,	0.089,	0.091,	0.088,	0.09,	0.092]]


Si_conc_array = [[0.077,	0.079,	0.082,	0.079,	0.082,	0.082],
[0.080,	0.079,	0.08,	0.079,	0.082,	0.082],
[0.079,	0.080,	0.083,	0.079,	0.079,	0.083],
[0.079,	0.079,	0.079,	0.079,	0.079,	0.084],
[0.079,	0.080,	0.079,	0.079,	0.079,	0.085]]

#%%  Plotting
fig, ax = plt.subplots()

# Function to show the heat map
im = ax.imshow( tau_decay_array , cmap = 'magma' )

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(Alpha_c_array)), labels=Alpha_c_array)
ax.set_yticks(np.arange(len(Delta_implant_array)), labels=Delta_implant_array)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
  
# Loop over data dimensions and create text annotations.
# for i in range(len(Alpha_c_array)):
#     for j in range(len(Delta_implant_array)):
#         text = ax.text(j, i, tau_decay_array[i, j],
#                        ha="center", va="center", color="w")

# Adding details to the plot
ax.set_title( "Characteristic Decay Time (s)",fontsize=20)
ax.set_xlabel(r"$\alpha_C$",fontsize=18,color='r')
ax.set_ylabel(r"$\Delta_{implant} (nm)$",fontsize=18,color='r')

fig.tight_layout()


# Adding a color bar to the plot
#ax.colorbar()
fig.colorbar(im)

plt.show()





fig, ax = plt.subplots()
im = ax.imshow( Gross_erosion_array , cmap = 'magma' )
ax.set_xticks(np.arange(len(Alpha_c_array)), labels=Alpha_c_array)
ax.set_yticks(np.arange(len(Delta_implant_array)), labels=Delta_implant_array)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_title( r"Si Gross Erosion ($m^{-2}s^{-1}$)",fontsize=20)
ax.set_xlabel(r"$\alpha_C$",fontsize=18,color='r')
ax.set_ylabel(r"$\Delta_{implant} (nm)$",fontsize=18,color='r')
fig.tight_layout()
fig.colorbar(im)
plt.show()




fig, ax = plt.subplots()
im = ax.imshow( SiC_conc_array , cmap = 'magma' )
ax.set_xticks(np.arange(len(Alpha_c_array)), labels=Alpha_c_array)
ax.set_yticks(np.arange(len(Delta_implant_array)), labels=Delta_implant_array)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_title( r"Silicon Carbide Fraction",fontsize=20)
ax.set_xlabel(r"$\alpha_C$",fontsize=18,color='r')
ax.set_ylabel(r"$\Delta_{implant} (nm)$",fontsize=18,color='r')
fig.tight_layout()
fig.colorbar(im)
plt.show()



fig, ax = plt.subplots()
im = ax.imshow( C_conc_array , cmap = 'magma' )
ax.set_xticks(np.arange(len(Alpha_c_array)), labels=Alpha_c_array)
ax.set_yticks(np.arange(len(Delta_implant_array)), labels=Delta_implant_array)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_title( r"Carbon Fraction",fontsize=20)
ax.set_xlabel(r"$\alpha_C$",fontsize=18,color='r')
ax.set_ylabel(r"$\Delta_{implant} (nm)$",fontsize=18,color='r')
fig.tight_layout()
fig.colorbar(im)
plt.show()




fig, ax = plt.subplots()
im = ax.imshow( Si_conc_array , cmap = 'magma' )
ax.set_xticks(np.arange(len(Alpha_c_array)), labels=Alpha_c_array)
ax.set_yticks(np.arange(len(Delta_implant_array)), labels=Delta_implant_array)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_title( r"Silicon Fraction",fontsize=20)
ax.set_xlabel(r"$\alpha_C$",fontsize=18,color='r')
ax.set_ylabel(r"$\Delta_{implant} (nm)$",fontsize=18,color='r')
fig.tight_layout()
fig.colorbar(im)
plt.show()




