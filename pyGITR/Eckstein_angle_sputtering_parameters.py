#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 19:14:41 2022

@author: de
"""

import numpy as np

# Creating a dictionary for the rotation parameters

List_Paramters = ['E0','f','b','c','Theta0_star']




#  Silicon to Carbon -- using RUSTBCA (done)
Sputtering_Rotation_Si_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_Si_C['E0'] = [10, 27.82, 46.4159, 129.1550, 215.4435, 359.3814, 599.4853]

Sputtering_Rotation_Si_C['f'] =  [68, 63, 63, 64, 64, 50, 29 ]

Sputtering_Rotation_Si_C['b'] =  [51, 49, 49, 49, 49, 37, 25  ]

Sputtering_Rotation_Si_C['c'] =  [0.8, 0.8, 0.8, 0.8, 0.8, 1.4,  0.1]

Sputtering_Rotation_Si_C['Theta0_star'] = [220, 220, 220, 220, 220, 220 , 220]





#  Carbon to Silicon -- using RUSTBCA (done)
Sputtering_Rotation_C_Si = dict.fromkeys(List_Paramters)

Sputtering_Rotation_C_Si['E0'] = [10, 27.82, 46.4159, 129.1550, 215.4435, 359.3814, 599.4853]

Sputtering_Rotation_C_Si['f'] =  [40, 52, 43, 44, 44, 37, 12]

Sputtering_Rotation_C_Si['b'] =  [10, 38, 33, 33, 33, 28, 10]

Sputtering_Rotation_C_Si['c'] =  [0.05, 1.0, 1.2, 1.2, 1.2, 1.4, 0.1]

Sputtering_Rotation_C_Si['Theta0_star'] = [90.56, 200.56, 200.56, 200.56, 200.56, 200.56, 200.56 ]





#  Hydrogen to Silicon   -- using RUSTBCA (very bad fit -- done)
Sputtering_Rotation_H_Si = dict.fromkeys(List_Paramters)

Sputtering_Rotation_H_Si['E0'] = [10, 27.82, 46.4159, 129.1550, 215.4435, 359.3814, 599.4853]

Sputtering_Rotation_H_Si['f'] =  [26, 17, 40, 38, 38 , 39, 37  ]

Sputtering_Rotation_H_Si['b'] =  [20, 18, 25.8, 26, 26, 26, 26  ]

Sputtering_Rotation_H_Si['c'] =  [0.8,  0.5, 0.08, 0.08, 0.08, 0.08 , 0.08]

Sputtering_Rotation_H_Si['Theta0_star'] = [190, 260, 260, 260, 260, 260, 260  ]

#  Hydrogen to Carbon

Sputtering_Rotation_H_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_H_C['E0'] = [40, 50, 70, 100, 140, 200, 300, 500, 1000, 2000]

Sputtering_Rotation_H_C['f'] = [25.1377, 11.3083, 6.2750, 4.4547, 4.1224, 4.0616, 3.9732, 4.5761, 4.9595, 5.2103]

Sputtering_Rotation_H_C['b'] = [16.8180, 7.6579, 3.4977, 1.7320, 1.2451, 1.0525, 0.9309, 1.0786, 1.1242, 1.1392]

Sputtering_Rotation_H_C['c'] = [0.8164, 0.7594, 0.9006, 1.0321, 1.0706, 1.0682, 1.0465, 0.9697, 0.8424, 0.8804]

Sputtering_Rotation_H_C['Theta0_star'] = [98.98, 98.05, 96.82, 95.71, 94.83, 94.04, 93.30, 92.56, 91.81, 91.28]


#  Deuterium to Carbon

Sputtering_Rotation_D_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_D_C['E0'] = [30, 40, 50, 70, 100, 140, 200, 300, 350, 500, 1000, 2000]

Sputtering_Rotation_D_C['f'] = [18.7020, 13.9859, 9.1046, 6.8357, 6.1142, 5.5800, 5.4183, 5.4718, 5.3398, 5.5508, 5.2293, 5.1506]

Sputtering_Rotation_D_C['b'] = [12.7582, 8.3859, 4.7171, 2.9145, 2.1925, 1.7856, 1.6184, 1.5643, 1.4866, 1.5347, 1.2861, 1.1281]

Sputtering_Rotation_D_C['c'] = [0.5994, 0.7022, 0.8423, 0.9611, 1.0137, 1.0228, 1.0124, 0.9764, 0.9908, 0.9123, 0.8339, 0.8779]

Sputtering_Rotation_D_C['Theta0_star'] = [100.35, 98.98, 98.05, 96.82, 95.71, 94.83, 94.04, 93.30, 93.06, 92.56, 91.81, 91.28]



#  Tritium to Carbon

Sputtering_Rotation_T_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_T_C['E0'] = [25, 30, 40, 50, 70, 100, 140, 200, 300, 500, 1000]

Sputtering_Rotation_T_C['f'] = [22.0075, 16.7368, 12.8519, 11.1855, 8.9010, 8.2697, 7.2448, 7.2269, 6.7254, 6.1973, 6.3026]

Sputtering_Rotation_T_C['b'] = [14.0078, 9.7905, 6.5908, 5.2720, 3.6518, 3.1125, 2.4653, 2.3719, 2.0851, 1.7780, 1.7218]

Sputtering_Rotation_T_C['c'] = [0.5594, 0.6356, 0.7794, 0.8488, 0.9388, 0.9528, 0.9755, 0.9540, 0.9294, 0.8939, 0.7847]

Sputtering_Rotation_T_C['Theta0_star'] = [101.31, 100.35, 98.98, 98.05, 96.82, 95.71, 94.83, 94.04, 93.30, 92.56, 91.81]



#  Carbon to Carbon


Sputtering_Rotation_C_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_C_C['E0'] = [25, 30, 40, 50, 70, 100, 140, 200, 300, 500, 1000, 3000]

Sputtering_Rotation_C_C['f'] = [50.5913, 48.6907, 42.9993, 38.3944, 32.5113, 27.6668, 23.6811, 20.2369, 17.3540, 14.4607, 11.7081, 7.4737]

Sputtering_Rotation_C_C['b'] = [26.2291, 25.5043, 22.5078, 19.9900, 16.7415, 14.0728, 11.8531, 9.8676, 8.1805, 6.4252, 4.7573, 2.4240]

Sputtering_Rotation_C_C['c'] = [0.4926, 0.5177, 0.5206, 0.5194, 0.5356, 0.5588, 0.5893, 0.6211, 0.6569, 0.7055, 0.7358, 0.8103]

Sputtering_Rotation_C_C['Theta0_star'] = [118.57, 116.44, 113.29, 111.06, 108.02, 105.23, 102.96, 100.90, 98.93, 96.94, 94.92, 92.85]


#  Deuterium to Silicon

Sputtering_Rotation_D_Si = dict.fromkeys(List_Paramters)

Sputtering_Rotation_D_Si['E0'] = [30, 50, 100, 500, 1000]

Sputtering_Rotation_D_Si['f'] = [44.6047, 25.2643, 12.7099, 6.3060, 4.6790]

Sputtering_Rotation_D_Si['b'] = [26.3734, 14.1721, 6.0832, 2.0569, 1.2020]

Sputtering_Rotation_D_Si['c'] = [0.3164, 0.5693, 0.6848, 0.8725, 0.8521]

Sputtering_Rotation_D_Si['Theta0_star'] = [100.35, 98.05, 95.71, 92.56, 91.81]





# Silicon to Silicon


Sputtering_Rotation_Si_Si = dict.fromkeys(List_Paramters)

Sputtering_Rotation_Si_Si['E0'] = [200, 500, 2000] # Si(KrC) for 200 eV

Sputtering_Rotation_Si_Si['f'] = [15.7808, 11.1924, 8.2975]

Sputtering_Rotation_Si_Si['b'] = [7.9103, 4.9698, 3.1245]

Sputtering_Rotation_Si_Si['c'] = [0.6307, 0.7029, 0.7702]

Sputtering_Rotation_Si_Si['Theta0_star'] = [98.72, 95.54, 92.78]


#  Hydrogen to Tungsten

Sputtering_Rotation_H_W = dict.fromkeys(List_Paramters)

Sputtering_Rotation_H_W['E0'] = [500, 550, 600, 700, 800, 900, 1000, 2000, 4000]

Sputtering_Rotation_H_W['f'] = [2.9896, 2.9665, 2.1320, 1.8691, 2.2265, 1.6453, 1.4684, 1.2970, 1.7925]

Sputtering_Rotation_H_W['b'] = [1.5460, 1.5154, 1.0397, 0.8504, 0.9849, 0.6802, 0.5018, 0.1244, 0.1535]

Sputtering_Rotation_H_W['c'] = [0.9382, 0.9223, 0.9630, 0.9839, 0.9428, 0.9714, 1.0082, 1.0630, 1.0319]

Sputtering_Rotation_H_W['Theta0_star'] = [92.56, 92.44, 92.34, 92.16, 92.02, 91.91, 91.81, 91.28, 90.91]

#  Deuterium to Tungsten

Sputtering_Rotation_D_W = dict.fromkeys(List_Paramters)

Sputtering_Rotation_D_W['E0'] = [250, 270, 300, 350, 400, 500, 600, 700, 1000]

Sputtering_Rotation_D_W['f'] = [3.9269, 3.3257, 2.4382, 2.2277, 1.9710, 1.5306, 1.0222, 1.2293, 1.2531]

Sputtering_Rotation_D_W['b'] = [2.5577, 1.9540, 1.4079, 1.2410, 1.0176, 0.6843, 0.3585, 0.3613, 0.2141]

Sputtering_Rotation_D_W['c'] = [0.8242, 0.9277, 0.9691, 0.9674, 0.9898, 1.0177, 1.0601, 1.0416, 1.0543]

Sputtering_Rotation_D_W['Theta0_star'] = [93.62, 93.48, 93.30, 93.06, 92.86, 92.56, 92.34, 92.16, 91.81]


#  Tritium to Tungsten

Sputtering_Rotation_T_W = dict.fromkeys(List_Paramters)

Sputtering_Rotation_T_W['E0'] = [170, 180, 200, 250, 300, 400, 500, 700, 1000]

Sputtering_Rotation_T_W['f'] = [4.9883, 2.7166, 2.7572, 2.0253, 1.7568, 1.4228, 1.3307, 1.0453, 1.1420]

Sputtering_Rotation_T_W['b'] = [3.2765, 1.7643, 1.7937, 1.2204, 0.9622, 0.6303, 0.4476, 0.1849, 0.1278]

Sputtering_Rotation_T_W['c'] = [0.8428, 0.9784, 0.9434, 0.9902, 1.0136, 1.0284, 1.0410, 1.0952, 1.1154]

Sputtering_Rotation_T_W['Theta0_star'] = [94.39, 94.26, 94.04, 93.62, 93.30, 92.86, 92.56, 92.16, 91.81]


#  Tungsten to Tungsten

Sputtering_Rotation_W_W = dict.fromkeys(List_Paramters)

Sputtering_Rotation_W_W['E0'] = [35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 300, 350, 400, 500, 800, 1000, 2000, 2500, 5000]

Sputtering_Rotation_W_W['f'] = [38.9214, 38.7481, 36.5132, 34.9177, 31.7802, 29.3956, 25.6694, 23.3024, 21.4480, 18.3650, 15.5941, 14.6047, 14.0761, 12.7467, 10.8164, 9.8530, 7.9586, 7.7404, 6.9226]

Sputtering_Rotation_W_W['b'] = [19.3927, 19.6182, 18.4183, 18.0217, 16.4397, 15.2327, 13.3290, 12.1541, 11.2129, 9.6309, 8.1239, 7.5244, 7.2438, 6.4829, 5.2902, 4.6927, 3.5143, 3.3309, 2.7745]

Sputtering_Rotation_W_W['c'] = [0.4803, 0.5111, 0.5099, 0.5214, 0.5209, 0.5224, 0.5375, 0.5347, 0.5469, 0.5779, 0.6203, 0.6377, 0.6391, 0.6683, 0.7212, 0.7465, 0.6958, 0.7891, 0.7989]

Sputtering_Rotation_W_W['Theta0_star'] = [116.47, 114.98, 112.62, 110.82, 109.40, 108.23, 106.42, 105.05, 103.98, 101.77, 99.65, 98.95, 98.38, 97.51, 95.95, 95.32, 93.77, 93.37, 92.39]





#  Carbon to Tungsten

Sputtering_Rotation_C_W = dict.fromkeys(List_Paramters)

Sputtering_Rotation_C_W['E0'] = [70, 100, 150, 200, 250, 300, 350, 400, 500, 600, 1000]

Sputtering_Rotation_C_W['f'] = [8.64058, 4.40983, 3.10475, 2.38205, 2.15983, 2.24701, 2.30960, 2.44819, 2.63258, 2.76898, 2.90335]

Sputtering_Rotation_C_W['b'] = [7.82604, 4.25497, 2.68550, 1.84344, 1.48361, 1.41341, 1.36178, 1.38308, 1.40376, 1.42265, 1.34780]

Sputtering_Rotation_C_W['c'] = [0.41046, 0.56592, 0.62580, 0.72769, 0.80484, 0.83535, 0.85930, 0.86511, 0.87170, 0.87466, 0.88860]

temp_array = (np.pi - np.arccos(np.sqrt(1/(1+np.array(Sputtering_Rotation_C_W['E0'][:])/8.68)))) # E_sp from Carolina Bjo ̈rkas

Sputtering_Rotation_C_W['Theta0_star'] = temp_array.tolist()




#  Tungsten to Carbon

Sputtering_Rotation_W_C = dict.fromkeys(List_Paramters)

Sputtering_Rotation_W_C['E0'] = [150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000]

Sputtering_Rotation_W_C['f'] = [42.17052, 33.54171, 28.41065, 24.37823, 21.02544, 18.35384, 16.24075, 14.54318, 12.33922, 11.08408, 10.38992, 9.69051]

Sputtering_Rotation_W_C['b'] = [18.94823, 15.04655, 12.62815, 10.66071, 8.99644, 7.68049, 6.65058, 5.82076, 4.76094, 4.14915, 3.80177, 3.42556]

Sputtering_Rotation_W_C['c'] = [0.49332, 0.54302, 0.57918, 0.61419, 0.65045, 0.68349, 0.71204, 0.73815, 0.77527, 0.79950, 0.81388, 0.82930]

temp_array = (np.pi - np.arccos(np.sqrt(1/(1+np.array(Sputtering_Rotation_W_C['E0'][:])/7.42)))) # E_sp from Carolina Bjo ̈rkas

Sputtering_Rotation_W_C['Theta0_star'] = temp_array.tolist()








































































