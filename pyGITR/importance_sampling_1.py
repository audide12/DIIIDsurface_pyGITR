#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Fri Mar 10 12:54:45 2023
"""
import numpy as np
from pyGITR.math_helper import *
from typing import Callable
import matplotlib.pyplot as plt
import pydoc
import netCDF4
import os

def Gaussian(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 10.0, mu: float = 120.0, beta: float = 1.0, Normalized=True):
     f = np.exp(-((x-mu)/sigma)**2)
     
     f[np.argwhere(x<0)] = 0
     
     if Normalized:
         f = f/Integrale(f, x, Array=False)
     return f

def Gaussian_test(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 10.0, mu: float = 130, beta: float = 1.0, Normalized=True):
     f = np.exp(-((x-mu)/sigma)**2)
      
     f[np.argwhere(x<0)] = 0
      
     if Normalized:
         f = f/Integrale(f, x, Array=False)
     return f



def Thomson(x: np.ndarray = np.linspace(0, 300, 10000), xb: float = 8.64, xc: float = 100, Normalized=True):
    assert not (xc <= xb), "xc cannot be <= xb"
    f = x/(x + xb) ** 3*(1.0-np.sqrt((x+xb)/(xc+xb)))
    f[np.argwhere(x > xc)] = 0.0
    # if Normalized:
    #     f = f/Integrale(f, x, Array=False)
    return f


def Levy(x=np.linspace(0.1,10,10000), c=1, mu=0):
    return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)

x = np.linspace(1, 200, 5005)


pdf_p = Gaussian(x)
pdf_q = Gaussian_test(x)
#print(pdf_q)
plt.figure()
#plt.plot(x,pdf_p*pdf_q,label="p times q distribution")
plt.plot(x,Integrale(np.multiply(pdf_p,1.0), x, Array=True),label="cumulative of p times q")
# plt.plot(x,pdf_q,label="sampling distribution")
plt.legend()
plt.show()


# Gaussian(np.array([x]))
# Gaussian_Jerome(np.array([x]))





