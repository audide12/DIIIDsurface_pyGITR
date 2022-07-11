#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 20:28:29 2022

@author: audide
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:50:33 2022

@author: audide
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generation of particles distribution for GITR.
@author: guterl
"""
import pyGITR
from pyGITR.Particles import *


class Physical_Sputtering_Reflection_Plots():
    

    @classmethod
    def Sputtering_yields(cls, Projectile, Target):
        
        s_plot = Sputtering_and_reflection()

        Energies = np.linspace(460,1e5,10000)
        Sputtering = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Sputtering[counter] = s_plot.Calculate_PhysicalSputteringParameters(Projectile,Target,energy).real
            counter = counter + 1

        return Sputtering
    
    @classmethod
    def Reflection_yields(cls, Projectile, Target):
        
        s_plot = Sputtering_and_reflection()

        Energies = np.linspace(20,1e5,10000)
        Reflection = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Reflection[counter] = s_plot.Calculate_ReflectionCoefficients(Projectile,Target,energy).real
            counter = counter + 1

        return Reflection
       #plt.plot(Energies,Sputtering)
       
       
       
       
       
       

Energies = np.linspace(0.1,1e5,10000)
Sputtering = Physical_Sputtering_Reflection_Plots.Sputtering_yields('C', 'C')


plt.loglog(Energies,Sputtering.real,'ro')
plt.xlim(0.1,1e5)
plt.ylim(1e-5,1)

plt.xlabel('E_incident(eV)')
plt.ylabel('Sputering Yield (Y)')
plt.title(r'$C \rightarrow C$')


# Reflection = Physical_Sputtering_Reflection_Plots.Reflection_yields('C', 'C')


# plt.loglog(Energies,Reflection.real,'ro')
# plt.xlim(20,1e5)
# #plt.ylim(1e-2,15)

# plt.xlabel('E_incident(eV)')
# plt.ylabel('Reflection Coefficients (Y)')
# plt.title(r'$C \rightarrow C$')
