#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 13:46:16 2023

@author: de
"""

import numpy as np
import matplotlib.pyplot as plt
from pyGITR.math_helper import *

def Gaussian(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 5.0, mu: float = 10.0, beta: float = 1.0, Normalized=True):
    f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
    
    if Normalized:
        f = f/Integrale(f, x, Array=False)
    return f

def Gaussian_test(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 15.0, mu: float = 30.0, beta: float = 1.0, Normalized=True):
    f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
    
    if Normalized:
        f = f/Integrale(f, x, Array=False)
    return f


def Gaussian_Jerome(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 20.0, mu: float = 120, beta: float = 0.0, Normalized=True):
    f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
    # if beta > 0:
    #     f[np.argwhere(x<0)] = 0
    if Normalized:
        f = f/Integrale(f, x, Array=False)
    return f

def Thomson(x: np.ndarray = np.linspace(0, 300, 10000), xb: float = 8.64, xc: float = 200, Normalized=True):
    assert not (xc <= xb), "xc cannot be <= xb"
    f = x/(x + xb) ** 3*(1.0-np.sqrt((x+xb)/(xc+xb)))
    
    f[np.argwhere(x > xc)] = 0.0
    if Normalized:
        f = f/Integrale(f, x, Array=False)
    return f


rng = np.random.default_rng()
sigma_x = 10
x = np.linspace(0, 5*sigma_x, 100000)

g_ = Gaussian(x)   #
gw_ = Gaussian_test(x)  #

G_ = Gaussian(x) * Gaussian_test(x)  #

I = Integrale(g_, x, Array=True)
Iw = Integrale(G_, x, Array=True)
G_ = G_/Iw[-1]
Iw = Iw/Iw[-1]


n= 100000
xi = rng.random(n)
# xi = np.linspace(0.1,0.9,n)
v = np.interp(xi, I, x)
vw = np.interp(xi, Iw, x)

sampled_x_physics = x[np.searchsorted(Integrale(g_, x, Normalized=True), xi, side='left')]
sampled_x = x[np.searchsorted(Integrale(G_, x, Normalized=True), xi, side='left')]

weights = 1/gw_[np.searchsorted(Integrale(G_, x, Normalized=True), xi, side='left')]

wg = 1/Gaussian_test(vw)     #

# plt.figure()
# plt.plot(xi,wg)
# plt.show()
#%%
plt.figure()
plt.hist(sampled_x_physics,density=True,bins=100,alpha=0.5,label="physics")
plt.hist(sampled_x ,density=True,bins=100, weights=weights,label="weighted samples")
plt.legend()
plt.show()
#%%


# plt.figure()
# #plt.plot(x,g_,label="p distribution")
# #plt.plot(x,gw_,label="q distribution")
# plt.plot(x,G_)
# plt.legend()
# plt.show()
#%%
plt.figure()
plt.hist(wg,bins=100)
plt.hist(weights,bins=100,alpha=0.5)
plt.show()


