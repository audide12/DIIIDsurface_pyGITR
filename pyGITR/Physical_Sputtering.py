import numpy as np
from pyGITR.math_helper import *
from typing import Callable
import matplotlib.pyplot as plt
import pydoc
import netCDF4
import os


class Sputtering_and_reflection():


    def ShowAvailableProjectiles(self):
        for D in ['H', 'D', 'T', 'He4', 'Si', 'C', 'W']:
            print(D)
            
            
            
    def ShowAvailableTargets(self):
        for D in ['C', 'Si', 'W']:
            print(D)



    # these parameter values are hard-coded. Want to put it in a text file separately 
    @classmethod
    def Set_PhysicalSputteringParameters(cls, Projectile, Target):
        if Projectile == 'H' and Target =='C':
            cls.lambda_parameter = 1.3533
            cls.q_parameter = 0.0241
            cls.mu_parameter = 1.4103
            cls.Eth_parameter = 38.630 # in eV
        elif Projectile == 'D' and Target =='C':
            cls.lambda_parameter = 1.2848
            cls.q_parameter = 0.0539
            cls.mu_parameter = 1.1977
            cls.Eth_parameter = 27.770 # in eV
        elif Projectile == 'T' and Target =='C':
            cls.lambda_parameter = 1.9050
            cls.q_parameter = 0.0718
            cls.mu_parameter = 1.1512
            cls.Eth_parameter = 23.617 # in eV
        elif Projectile == 'He4' and Target =='C':
            cls.lambda_parameter = 4.5910
            cls.q_parameter = 0.1951
            cls.mu_parameter = 1.7852
            cls.Eth_parameter = 19.124 # in eV
        elif Projectile == 'C' and Target =='C':
            cls.lambda_parameter = 13.9666
            cls.q_parameter = 0.7015
            cls.mu_parameter = 2.0947
            cls.Eth_parameter = 21.4457 # in eV    
        elif Projectile == 'W' and Target =='C':  # Carolina Bjorkas
            cls.lambda_parameter = 1.2625
            cls.q_parameter = 1.3902
            cls.mu_parameter = 3.5395
            cls.Eth_parameter = 114.9398 # in eV 
           
        elif Projectile == 'Si' and Target =='C': # Find data
            print('Choose a different set of projectile and target')
        
        
        elif Projectile == 'H' and Target =='Si':
            cls.lambda_parameter = 0.4819
            cls.q_parameter = 0.0276
            cls.mu_parameter = 0.9951
            cls.Eth_parameter = 49.792 # in eV  
        elif Projectile == 'D' and Target =='Si':
            cls.lambda_parameter = 0.5326
            cls.q_parameter = 0.0569
            cls.mu_parameter = 1.6537
            cls.Eth_parameter = 24.543 # in eV  
        elif Projectile == 'T' and Target =='Si':
            cls.lambda_parameter = 0.4112
            cls.q_parameter = 0.0816
            cls.mu_parameter = 0.9325
            cls.Eth_parameter = 21.298 # in eV 
        elif Projectile == 'He4' and Target =='Si':
            cls.lambda_parameter = 0.2524
            cls.q_parameter = 0.2319
            cls.mu_parameter = 1.4732
            cls.Eth_parameter = 18.899 # in eV 
        elif Projectile == 'Si' and Target =='Si':
            cls.lambda_parameter = 0.6726
            cls.q_parameter = 2.6951
            cls.mu_parameter = 1.7584
            cls.Eth_parameter = 20.035 # in eV     
        elif Projectile == 'C' and Target =='Si':
            print('Choose a different set of projectile and target')
        elif Projectile == 'W' and Target =='Si':
            print('Choose a different set of projectile and target')    
        elif Projectile == 'H' and Target =='W':
            cls.lambda_parameter = 1.0087
            cls.q_parameter = 0.0075
            cls.mu_parameter = 1.2046
            cls.Eth_parameter = 457.42 # in eV   
        elif Projectile == 'D' and Target =='W':
            cls.lambda_parameter = 0.3583
            cls.q_parameter = 0.0183
            cls.mu_parameter = 1.4410
            cls.Eth_parameter = 228.84 # in eV   
        elif Projectile == 'T' and Target =='W':
            cls.lambda_parameter = 0.2870
            cls.q_parameter = 0.0419
            cls.mu_parameter = 1.5802
            cls.Eth_parameter = 153.8842 # in eV   
        elif Projectile == 'He4' and Target =='W':
            cls.lambda_parameter = 0.1692
            cls.q_parameter = 0.1151
            cls.mu_parameter = 1.7121
            cls.Eth_parameter = 120.56 # in eV   
        elif Projectile == 'C' and Target =='W':  # Carolina Bjorkas
            cls.lambda_parameter = 0.0447
            cls.q_parameter = 1.5622
            cls.mu_parameter = 1.0200 
            cls.Eth_parameter = 59.1980  # in eV   
        elif Projectile == 'W' and Target =='W': 
            cls.lambda_parameter = 2.2697
            cls.q_parameter = 18.6006
            cls.mu_parameter = 3.1273
            cls.Eth_parameter = 24.9885 # in eV   
        elif Projectile == 'Si' and Target =='W':
            print('Choose a different set of projectile and target')    
        
        
       
        
        
    @classmethod
    def Set_Mass_AtomicN(cls, Projectile, Target):
        if Projectile == 'H':    
            cls.Projectile_Mass = 1
            cls.Projectile_AtomicN = 1
        elif Projectile == 'D':
            cls.Projectile_Mass = 2.014
            cls.Projectile_AtomicN = 1
        elif Projectile == 'T':
            cls.Projectile_Mass = 3.016
            cls.Projectile_AtomicN= 1
        elif Projectile == 'He4':
            cls.Projectile_Mass = 4
            cls.Projectile_AtomicN = 2
        elif Projectile == 'Si':
            cls.Projectile_Mass = 28.0855
            cls.Projectile_AtomicN = 14
        elif Projectile == 'C':
            cls.Projectile_Mass = 12.011
            cls.Projectile_AtomicN = 6
        elif Projectile == 'W':
            cls.Projectile_Mass = 183.84
            cls.Projectile_AtomicN = 74
           
            
            
            
            
            
        if Target == 'C':    
            cls.Target_Mass = 12.011
            cls.Target_AtomicN = 6
        elif Target == 'Si':
            cls.Target_Mass = 28.0855
            cls.Target_AtomicN = 14
        elif Target == 'W':
            cls.Target_Mass = 183.84
            cls.Target_AtomicN = 74
           
            
            
            
    @classmethod        
    def Set_ReflectionParameters(cls, Projectile, Target):
        if Projectile == 'C' and Target =='C':
            cls.a_1 = -0.03022
            cls.a_2 = -1.107
            cls.a_3 = 6.544
            cls.a_4 = 0.1256
        elif Projectile == 'C' and Target =='W':
            cls.a_1 = 1.96
            cls.a_2 = 0.1
            cls.a_3 = 2.2
            cls.a_4 = 0.18   
        else:
            print('Choose a different set of projectile and target')    
            
    
            
    
    @classmethod
    def Calculate_PhysicalSputtering_RotationFactor(cls, Projectile, Target, Incident_Energy, Incident_Theta):
        
        
        flag = 1 # change it to zero if not found from the list 
        
        from scipy.interpolate import interp1d
        
        if Projectile == 'H' and Target =='C':
            Sputtering_Dictionary = Sputtering_Rotation_H_C
            
        elif Projectile == 'D' and Target =='C':
            Sputtering_Dictionary = Sputtering_Rotation_D_C
            
        elif Projectile == 'T' and Target =='C':
            Sputtering_Dictionary = Sputtering_Rotation_T_C    
        
        elif Projectile == 'C' and Target =='C':
            Sputtering_Dictionary = Sputtering_Rotation_C_C    
            
        elif Projectile == 'D' and Target =='Si':
            Sputtering_Dictionary = Sputtering_Rotation_D_Si    
                        
        elif Projectile == 'Si' and Target =='Si':
            Sputtering_Dictionary = Sputtering_Rotation_Si_Si
            
        elif Projectile == 'H' and Target =='W':
            Sputtering_Dictionary = Sputtering_Rotation_H_W
            
        elif Projectile == 'D' and Target =='W':
            Sputtering_Dictionary = Sputtering_Rotation_D_W
            
        elif Projectile == 'T' and Target =='W':
            Sputtering_Dictionary = Sputtering_Rotation_T_W
            
        elif Projectile == 'W' and Target =='W':
            Sputtering_Dictionary = Sputtering_Rotation_W_W  
            
        elif Projectile == 'W' and Target =='C':
            Sputtering_Dictionary = Sputtering_Rotation_W_C  
            
        elif Projectile == 'C' and Target =='W':
            Sputtering_Dictionary = Sputtering_Rotation_C_W      
            
        else:
            flag = 0
                
        if flag == 1:
            f_interpolate = interp1d(Sputtering_Dictionary['E0'], Sputtering_Dictionary['f'], kind='cubic')
            b_interpolate = interp1d(Sputtering_Dictionary['E0'], Sputtering_Dictionary['b'], kind='cubic')
            c_interpolate = interp1d(Sputtering_Dictionary['E0'], Sputtering_Dictionary['c'], kind='cubic')
            Theta0star_interpolate = interp1d(Sputtering_Dictionary['E0'], Sputtering_Dictionary['Theta0_star'], kind='cubic')
            
            cosine_factor = np.cos(np.pi*0.5*Incident_Theta/Theta0star_interpolate(Incident_Energy))**c_interpolate(Incident_Energy)
            
            factor = (cosine_factor)**(-f_interpolate(Incident_Energy)) * np.exp(b_interpolate(Incident_Energy)*(1-(1/cosine_factor) )) 
            
            return factor
        
        else:
            return 1
        
            
        
        
        
    
    
    def Calculate_PhysicalSputteringParameters(self, Projectile, Target, Incident_Energy, Incident_Theta:float=0):
        Sputtering_and_reflection.Set_PhysicalSputteringParameters(Projectile,Target)
        Sputtering_and_reflection.Set_Mass_AtomicN(Projectile,Target)
        
        electric_charge_square = 1.4399 # eV nm
        Bohr_radius = 0.0529177 # nm
        
        Lindhard_length = (9*np.pi**2/128)**(1/3) * Bohr_radius * (self.Target_AtomicN**(2/3)+ self.Projectile_AtomicN**(2/3))**(-0.5)  # nm 
        
        epsilon = Incident_Energy*(self.Target_Mass/(self.Projectile_Mass + self.Target_Mass))* (Lindhard_length/(self.Projectile_AtomicN*self.Target_AtomicN)*electric_charge_square) 
        Nuclear_stopping = (0.5*np.log(1+1.2288*epsilon))/(epsilon + 0.1728*np.sqrt(epsilon) + 0.008*epsilon**0.1504) 
        
        #Y = self.q_parameter*Nuclear_stopping*((Incident_Energy/self.Eth_parameter)-1)**self.mu_parameter/(self.lambda_parameter + ((Incident_Energy/self.Eth_parameter)-1)**self.mu_parameter)
        
        
        Numerator = ((Incident_Energy/self.Eth_parameter)-1)
        
        if Numerator < 0:
            Numerator = 0.0
        
        Y = self.q_parameter*Nuclear_stopping*Numerator**self.mu_parameter/(self.lambda_parameter + Numerator**self.mu_parameter)
        
<<<<<<< HEAD
        #Y = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, Target, Incident_Energy, Incident_Theta)*Y
=======
        
        if Incident_Theta != 0:    
            Y = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, Target, Incident_Energy, Incident_Theta)*Y
>>>>>>> refs/remotes/origin/main
        
        return Y
    
    def Calculate_ReflectionCoefficients(self, Projectile, Target, Incident_Energy):
        Sputtering_and_reflection.Set_ReflectionParameters(Projectile,Target)
        Sputtering_and_reflection.Set_Mass_AtomicN(Projectile,Target)
        
        electric_charge_square = 1.4399 # eV nm
        Bohr_radius = 0.0529177 # nm
        
        Lindhard_length = (9*np.pi**2/128)**(1/3) * Bohr_radius * (self.Target_AtomicN**(2/3)+ self.Projectile_AtomicN**(2/3))**(-0.5)  # nm 
        
        epsilon_L = ((self.Target_Mass + self.Projectile_Mass)/self.Target_Mass)*(self.Target_AtomicN*self.Projectile_AtomicN*electric_charge_square/Lindhard_length)
        
        epsilon = Incident_Energy/epsilon_L
        
        R_N = np.exp(self.a_1*epsilon**self.a_2)/(1+np.exp(self.a_3*epsilon**self.a_4))
    
        return R_N.real
    
    
    
    
#s = Sputtering_and_reflection()
#s.ShowAvailableTargets()
    
    
    
    
    
    
    