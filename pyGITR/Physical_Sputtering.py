import numpy as np
#from pyGITR.math_helper import *
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
        for D in ['C', 'Si', 'W', 'SiC','SiC,Si,Crystalline','SiC,C,Crystalline','SiC,Si,Amorphous','SiC,C,Amorphous','SiC,Si-crystalline','SiC,C-crystalline']:
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
           
        elif Projectile == 'Si' and Target =='C': # Fit with RUSTBCA simulations
            cls.lambda_parameter = 1.0
            cls.q_parameter = 1.8
            cls.mu_parameter = 0.3
            cls.Eth_parameter = 10.457 # in eV  
        
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
        elif Projectile == 'C' and Target =='Si':  # Fit with RUSTBCA simulations
            cls.lambda_parameter = 1.0
            cls.q_parameter = 1.4
            cls.mu_parameter = 12.0
            cls.Eth_parameter = 14 # in eV     
        elif Projectile == 'W' and Target =='Si':
            print('Choose a different set of projectile and target, data not available for W->Si')    
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
            
        elif Projectile == 'C' and Target =='SiC':  # from RustBCA
            cls.lambda_parameter = 1.0
            cls.q_parameter = 2.0
            cls.mu_parameter = 1.3
            cls.Eth_parameter = 25 # in eV   
        
        elif Projectile == 'Si' and Target =='SiC': # from RustBCA
            cls.lambda_parameter = 2.0
            cls.q_parameter = 4.0
            cls.mu_parameter = 1.3
            cls.Eth_parameter = 70 # in eV   
            
        elif Projectile == 'H' and Target =='SiC': # from RustBCA
            cls.lambda_parameter = 2.2
            cls.q_parameter = 0.031
            cls.mu_parameter = 1.5
            cls.Eth_parameter = 40 # in eV      
            
        elif Projectile == 'C' and Target =='SiC,Si,Crystalline':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 100.2
            cls.q_parameter = 0.41
            cls.mu_parameter = 2.1
            cls.Eth_parameter = 40 # in eV
            

        elif Projectile == 'C' and Target =='SiC,C,Crystalline':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 100.2
            cls.q_parameter = 0.85
            cls.mu_parameter = 1.6
            cls.Eth_parameter = 10 # in eV    
             
        elif Projectile == 'C' and Target =='SiC,Si,Amorphous':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 100.2
            cls.q_parameter = 1.65
            cls.mu_parameter = 2.6
            cls.Eth_parameter = 16 # in eV  
             
        elif Projectile == 'C' and Target =='SiC,C,Amorphous':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 100.2
            cls.q_parameter = 1.85
            cls.mu_parameter = 2.2
            cls.Eth_parameter = 10 # in eV
             
             
        elif Projectile == 'D' and Target =='SiC,Si,Crystalline':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 100.2
            cls.q_parameter = 0.013
            cls.mu_parameter = 2
            cls.Eth_parameter = 47 # in eV    
             
        elif Projectile == 'D' and Target =='SiC,C,Crystalline':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 160.2
            cls.q_parameter = 0.0365
            cls.mu_parameter = 2.4
            cls.Eth_parameter = 20 # in eV    
        
        elif Projectile == 'D' and Target =='SiC,Si,Amorphous':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 220.2
            cls.q_parameter = 0.09
            cls.mu_parameter = 2.4
            cls.Eth_parameter = 10 # in eV    
            
        elif Projectile == 'D' and Target =='SiC,C,Amorphous':  # for 45 incidence from SDTrimSP
            cls.lambda_parameter = 200.2
            cls.q_parameter = 0.11
            cls.mu_parameter = 2
            cls.Eth_parameter = 5 # in eV    
            
        elif Projectile == 'D' and Target =='SiC,Si-crystalline':  # for 45 incidence from SDTrimSP, averaged over three crystal orientations
            cls.lambda_parameter = 30.2
            cls.q_parameter = 0.0135
            cls.mu_parameter = 2
            cls.Eth_parameter = 80 # in eV      

        elif Projectile == 'D' and Target =='SiC,C-crystalline':  # for 45 incidence from SDTrimSP, averaged over three crystal orientations
            cls.lambda_parameter = 550
            cls.q_parameter = 0.04
            cls.mu_parameter = 2.8
            cls.Eth_parameter = 14 # in eV
            
        elif Projectile == 'C' and Target =='SiC,Si-crystalline':  # for 45 incidence from SDTrimSP, averaged over three crystal orientations
            cls.lambda_parameter = 80.2
            cls.q_parameter = 0.39
            cls.mu_parameter = 2.1
            cls.Eth_parameter = 40 # in eV            
            
        elif Projectile == 'C' and Target =='SiC,C-crystalline':  # for 45 incidence from SDTrimSP, averaged over three crystal orientations
            cls.lambda_parameter = 100.2
            cls.q_parameter =  0.88
            cls.mu_parameter = 1.8
            cls.Eth_parameter = 10 # in eV        
            
        elif Projectile == 'Si' and Target =='W':
            print('Choose a different set of projectile and target, data not available for Si-> W')    
        else:
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
        elif Target == 'SiC':
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20    
        elif Target == 'SiC,Si,Crystalline':
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20
        elif Target == 'SiC,C,Crystalline':
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20    
        elif Target == 'SiC,Si,Amorphous':
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20    
        elif Target == 'SiC,C,Amorphous':
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20    
        elif Target == 'SiC,Si-crystalline':# averaged over three crystal orientations
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20   
        elif Target == 'SiC,C-crystalline':# averaged over three crystal orientations
            cls.Target_Mass = 40
            cls.Target_AtomicN = 20      
            
           
            
         
            
            
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
        elif Projectile == 'W' and Target =='W':
            cls.a_1 = -3.685
            cls.a_2 = 0.02781
            cls.a_3 = 0.7825e-4
            cls.a_4 = -1.109
        elif Projectile == 'Si' and Target =='Si':
            cls.a_1 = -3.908
            cls.a_2 = -0.1721
            cls.a_3 = 1.974
            cls.a_4 = -0.5695 
        elif Projectile == 'D' and Target =='C':
            cls.a_1 = 0.1526
            cls.a_2 = -0.2304
            cls.a_3 = 0.2113
            cls.a_4 = 1.287
        elif Projectile == 'D' and Target =='Si':
            cls.a_1 = 0.2381
            cls.a_2 = -0.1662
            cls.a_3 = 0.1552
            cls.a_4 = 1.535   
            
        else:
            print('Choose a different set of projectile and target reflection')    
            
    
            
    
    @classmethod
    def Calculate_PhysicalSputtering_RotationFactor(cls, Projectile, Target, Incident_Energy, Incident_Theta):
        
        
        flag = 0 # change it to zero if not found from the list 
        
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
        # CHANGE THESE    
        elif Projectile == 'D' and Target =='SiC':    # verify these
            Sputtering_Dictionary = Sputtering_Rotation_H_SiC      
        
        elif Projectile == 'C' and Target =='SiC':
            Sputtering_Dictionary = Sputtering_Rotation_C_SiC     
            
            
        
        
        else:
            flag = 1
                
        if flag == 0:
            E0_array =  np.array(Sputtering_Dictionary['E0'])
            f_array = np.array(Sputtering_Dictionary['f'])
            b_array = np.array(Sputtering_Dictionary['b'])
            c_array = np.array(Sputtering_Dictionary['c'])
            Theta0_star_array = np.array(Sputtering_Dictionary['Theta0_star'])
            
            
            f_interpolate = interp1d(E0_array,f_array,bounds_error=False,fill_value="extrapolate")
            b_interpolate = interp1d(E0_array, b_array,bounds_error=False,fill_value="extrapolate")
            c_interpolate = interp1d(E0_array, c_array,bounds_error=False,fill_value="extrapolate")
            Theta0star_interpolate = interp1d(E0_array, Theta0_star_array,bounds_error=False,fill_value="extrapolate")
            
            cosine_factor = np.cos((np.pi*0.5*Incident_Theta/Theta0star_interpolate(Incident_Energy))**c_interpolate(Incident_Energy))
            
            factor = (cosine_factor)**(-f_interpolate(Incident_Energy)) * np.exp(b_interpolate(Incident_Energy)*(1-(1/cosine_factor) )) 
            
            if factor>10:
                factor = 1
            elif np.isnan(factor):
                factor = 1            
            #print("factor is   ",factor)
            return factor
        
        else:
            return 1
        
            
        
        
        
    
    
    def Calculate_PhysicalSputteringParameters(self, Projectile, Target, Incident_Energy, Incident_Theta:float=0.0):
        Sputtering_and_reflection.Set_PhysicalSputteringParameters(Projectile,Target)
        Sputtering_and_reflection.Set_Mass_AtomicN(Projectile,Target)
        
        electric_charge_square = 1.4399 # eV nm
        Bohr_radius = 0.0529177 # nm
        
        Lindhard_length = (9*np.pi**2/128)**(1/3) * Bohr_radius * (self.Target_AtomicN**(2/3)+ self.Projectile_AtomicN**(2/3))**(-0.5)  # nm 
        
        epsilon = Incident_Energy*(self.Target_Mass/(self.Projectile_Mass + self.Target_Mass))* (Lindhard_length/(self.Projectile_AtomicN*self.Target_AtomicN)*electric_charge_square) 
        Nuclear_stopping = (0.5*np.log(1+1.2288*epsilon))/(epsilon + 0.1728*np.sqrt(epsilon) + 0.008*epsilon**0.1504) 
        
        Y = self.q_parameter*Nuclear_stopping*((Incident_Energy/self.Eth_parameter)-1)**self.mu_parameter/(self.lambda_parameter + ((Incident_Energy/self.Eth_parameter)-1)**self.mu_parameter)
        
        
        Numerator = ((Incident_Energy/self.Eth_parameter)-1)
        
        if Numerator < 0:
            Numerator = 0.0

        Y = self.q_parameter*Nuclear_stopping*Numerator**self.mu_parameter/(self.lambda_parameter + Numerator**self.mu_parameter)
        
        if (Incident_Theta>0.0):
            Y = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, Target, Incident_Energy, Incident_Theta)*Y
            
        Rotation_database_factor = 1
        Rotation_database = 45
        
        # if Projectile == 'C' and Target =='SiC,Si,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'C' and Target =='SiC,Si,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor 
        # elif Projectile == 'C' and Target =='SiC,C,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'C' and Target =='SiC,C,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
            
            
        # elif Projectile == 'D' and Target =='SiC,Si,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'D' and Target =='SiC,Si,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor 
        # elif Projectile == 'D' and Target =='SiC,C,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'D' and Target =='SiC,C,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor    
        
        # elif Projectile == 'D' and Target =='SiC,Si,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'D' and Target =='SiC,Si,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor 
        # elif Projectile == 'D' and Target =='SiC,C,Crystalline':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor
        # elif Projectile == 'D' and Target =='SiC,C,Amorphous':  # for 45 incidence from SDTrimSP 
        #     Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'SiC', Incident_Energy, Rotation_database)*Rotation_database_factor    
            
        if Projectile == 'D' and Target =='Si-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
            Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'Si', Incident_Energy, Rotation_database)*Rotation_database_factor           
        elif Projectile == 'D' and Target =='C-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
            Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'C', Incident_Energy, Rotation_database)*Rotation_database_factor            
        elif Projectile == 'C' and Target =='Si-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
            Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'Si', Incident_Energy, Rotation_database)*Rotation_database_factor                        
        elif Projectile == 'C' and Target =='C-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
            Rotation_database_factor = Sputtering_and_reflection.Calculate_PhysicalSputtering_RotationFactor(Projectile, 'C', Incident_Energy, Rotation_database)*Rotation_database_factor
            
        
        Y = Y/Rotation_database_factor
        
        return Y
    
    
    def Calculate_ChemicalSputtering(self, Projectile, Target, Incident_Energy, Incident_Flux, Surface_Temperature): # ROTH Formula
        
        #print("WARNING: CHEMICAL SPUTTERING CURRENTLY ONLY FROM D->C")
        Sputtering_and_reflection.Set_Mass_AtomicN(Projectile,Target)
        
        k_B = 1.38 * 10**(-23)   # in m^2 kg s^-2 K^-1  Boltzmann constant
        
        ETHC = 27.0   #Deuterium
        ETFC = 447.0  #Deuterium    
        QC   = 0.1    #Deuterium
        ETHC_star = 1.0   #Deuterium
        D = 125   #Deuterium
        
        
        SNC = 0.5*np.log(1.+1.2288*Incident_Energy/ETFC)/(Incident_Energy/ETFC+ 0.1728*np.sqrt(Incident_Energy/ETFC)+ 0.008*pow(Incident_Energy/ETFC,0.1504))
                                                          
        if Incident_Energy>ETHC:
            YPHYS = QC*SNC*(1-(ETHC/Incident_Energy)**(2./3.))*(1-ETHC/Incident_Energy)**2
        else:
            YPHYS = 0.0
        
        CSURF  = 1/(1+1E13*np.exp(-2.45*11604/Surface_Temperature)) #Surface_Temperature in K
        CSP3   = CSURF*(2E-32*Incident_Flux+np.exp(-1.7*11604/Surface_Temperature))/(2E-32*Incident_Flux+(1+2E29/Incident_Flux*np.exp(-1.8*11604/Surface_Temperature))*np.exp(-1.7*11604/Surface_Temperature))
        
        YTHERM = CSP3*0.033*np.exp(-1.7*11604/Surface_Temperature)/(2E-32*Incident_Flux+np.exp(-1.7*11604/Surface_Temperature)) # Thermal Activated Process
        
        if Incident_Energy>ETHC_star:
            YSURF = CSP3*QC*SNC*(1-(ETHC_star/Incident_Energy)**(2./3.))*(1-ETHC_star/Incident_Energy)**2/(1.0+np.exp((min(90.0,Incident_Energy)-90)/50))
        else:
            YSURF = 0
            
        return YSURF + YTHERM * (1 + D * YPHYS)    
         
    
        
        
         
        
        
    
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
    
    
class Physical_Sputtering_Reflection_Plots():
    

    @classmethod
    def Sputtering_yields(cls, Projectile, Target,Energies):
        
        s_plot = Sputtering_and_reflection()

        #Energies = np.linspace(460,1e5,10000)
        Sputtering = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Sputtering[counter] = s_plot.Calculate_PhysicalSputteringParameters(Projectile,Target,energy).real
            counter = counter + 1

        return Sputtering
    
    @classmethod
    def Reflection_yields(cls, Projectile, Target, Energies):
        
        s_plot = Sputtering_and_reflection()

        #Energies = np.linspace(20,1e5,10000)
        Reflection = np.zeros(Energies.size)
        counter = 0
        
        
        for energy in Energies:
            
            
            Reflection[counter] = s_plot.Calculate_ReflectionCoefficients(Projectile,Target,energy).real
            counter = counter + 1

        return Reflection
       #plt.plot(Energies,Sputtering)
          
    
#s = Sputtering_and_reflection()
#s.ShowAvailableTargets()
    
    
# elif Projectile == 'D' and Target =='Si-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
#     cls.lambda_parameter = 400.2
#     cls.q_parameter = 0.036
#     cls.mu_parameter = 3
#     cls.Eth_parameter = 24 # in eV    
    
# elif Projectile == 'D' and Target =='C-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
#     cls.lambda_parameter = 10.2
#     cls.q_parameter = 0.1
#     cls.mu_parameter = 1
#     cls.Eth_parameter = 30 # in eV    
    
# elif Projectile == 'C' and Target =='Si-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
#     cls.lambda_parameter = 10.2
#     cls.q_parameter = 0.7
#     cls.mu_parameter = 2
#     cls.Eth_parameter = 40 # in eV    
    
# elif Projectile == 'C' and Target =='C-ENRICHEDSiC':  # for 45 incidence from SDTrimSP
#     cls.lambda_parameter = 100.2
#     cls.q_parameter = 1
#     cls.mu_parameter = 1.4
#     cls.Eth_parameter = 7 # in eV       
    
    
    
    