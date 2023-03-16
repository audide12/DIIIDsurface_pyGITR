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

def Gaussian(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 5.0, mu: float = 120.0, beta: float = 1.0, Normalized=True):
    f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
    
    f[np.argwhere(x<0)] = 0
    
    if Normalized:
        f = f/Integrale(f, x, Array=False)
    return f

def Gaussian_test(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 20.0, mu: float = 130.0, beta: float = 1.0, Normalized=True):
    f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
    
    f[np.argwhere(x<0)] = 0
    
    # if Normalized:
    #     f = f/Integrale(f, x, Array=False)
    return f


def Gaussian_Jerome(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 5.0, mu: float = 120, beta: float = 0.0, Normalized=True):
    f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
    if beta > 0:
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


pdf_p = Gaussian_Jerome(x)
pdf_q = Gaussian_test(x)

pdf_p_times_q = np.multiply(pdf_p,pdf_q)

#Normalization = 1.0
Normalization = Integrale(pdf_p_times_q, x, Array=False)
#Normalization = 1.0

pdf_p_times_q = np.divide(pdf_p_times_q,Normalization)


pdf_p_times_q_divide_q = np.divide(pdf_p_times_q,pdf_q)

#Normalization = Integrale(pdf_p_times_q_divide_q, x, Array=False)
#Normalization = 1.0
#print(Normalization)

pdf_p_times_q_divide_q = np.divide(pdf_p_times_q_divide_q,Normalization)


# #print(pdf_q)
plt.figure()
plt.plot(x,pdf_p,label="p distribution")
plt.plot(x,pdf_q,label="q distribution")
#plt.plot(x,pdf_p_times_q,label="multiplied distribution")
#plt.plot(x,Integrale(np.multiply(pdf_p,1.0), x, Array=True),label="cumulative of p times q")
#plt.plot(x,pdf_p_times_q_divide_q,label="multiplied distribution")
plt.xlim(100,150)
plt.legend()
plt.show()


# Gaussian(np.array([x]))
# Gaussian_Jerome(np.array([x]))





