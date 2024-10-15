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
#import pyGITR
from pyGITR.Particles import *

    
    
class Physical_Sputtering_Reflection_Plots():
    

    @classmethod
    def Sputtering_yields(cls, Projectile, Target,Energies: np.ndarray = np.linspace(460,1e5,10000)):
        
        s_plot = Sputtering_and_reflection()

        Sputtering = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Sputtering[counter] = s_plot.Calculate_PhysicalSputteringParameters(Projectile,Target,energy).real
            counter = counter + 1

        return Sputtering
    
    @classmethod
    def Sputtering_yields_angular(cls, Projectile, Target,energy:float = 100.0, Angles:np.ndarray = np.linspace(0,90,num=100)):
        
        s_plot = Sputtering_and_reflection()

        Sputtering = np.zeros(Angles.size)
        counter = 0
        
        
        for angle in Angles:
            
            
            Sputtering[counter] = s_plot.Calculate_PhysicalSputteringParameters(Projectile,Target,energy,angle).real
            counter = counter + 1

        return Sputtering
    
    @classmethod
    def Reflection_yields(cls, Projectile, Target,  Energies:np.ndarray = np.linspace(20,1e5,10000)):
        
        s_plot = Sputtering_and_reflection()

        Reflection = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Reflection[counter] = s_plot.Calculate_ReflectionCoefficients(Projectile,Target,energy).real
            counter = counter + 1

        return Reflection
       #plt.plot(Energies,Sputtering)
       
       
#%%       

import matplotlib.pyplot as plt              
plt.figure()

Projectile = 'H'
Target = 'Si'
       
Energies = np.linspace(0.1,1e3,10000)
Sputtering = Physical_Sputtering_Reflection_Plots.Sputtering_yields(Projectile,Target,Energies)

plt.loglog(Energies,Sputtering.real,'ro')
# plt.xlim(0.1,1e5)
# plt.ylim(1e-5,2)

plt.xlabel('E_incident(eV)')
plt.ylabel('Sputering Yield (Y)')
plt.title(r'$%s \rightarrow %s \quad $'%(Projectile,Target)+'(normal incidence)')


# Reflection = Physical_Sputtering_Reflection_Plots.Reflection_yields('C', 'W')


# plt.loglog(Energies,Reflection.real,'ro')
# plt.xlim(20,1e5)
# #plt.ylim(1e-2,15)

# plt.xlabel('E_incident(eV)')
# plt.ylabel('Reflection Coefficients (Y)')
# plt.title(r'$C \rightarrow W$')


#%%
import matplotlib.pyplot as plt       
plt.figure()       

Angles:np.ndarray = np.linspace(0,90,num=100)
Energy = 1000

Projectile = 'H'
Target = 'C'

Sputtering = Physical_Sputtering_Reflection_Plots.Sputtering_yields_angular(Projectile,Target,Energy,Angles)

plt.plot(Angles,Sputtering.real,'ro')


plt.xlabel('Angle of Incidence(degrees)')
plt.ylabel('Sputering Yield (Y)')
plt.title(r'$%s \rightarrow %s \quad (Energy = %.0f eV)$'%(Projectile,Target,Energy))

    






