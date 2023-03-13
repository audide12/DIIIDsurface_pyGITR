#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 13:46:16 2023

@author: de
"""

import matplotlib.pyplot as plt
import numpy as np
from pyGITR.Particles import *

def Levy(x=np.linspace(0.01,10,10000), c=1, mu=0):
    return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)


def weighted_particles(initial_distrib, sampling_distrib):
    
    weights = initial_distrib/sampling_distrib
    
    # plt.figure()
    # plt.hist(initial_distrib, bins=100, label="actual physics distribution")
    # plt.hist(desired_distrib, bins=100, label="sampling distribution")
    # plt.title(label="Comparison of distributions")
    # plt.xlim(90,110)
    # plt.legend()
    # plt.show()
    
    
    # plt.figure()
    # plt.hist(weights, bins=100, label="weights")
    # plt.title(label="Comparison of weights")
    # plt.legend()
    # plt.show()
    
    multiplied_distribution = np.multiply(weights,sampling_distrib)
    
    
    # plt.figure()
    # plt.hist(initial_distrib-multiplied_distribution, bins=10, label="actual physics distribution")
    # plt.title(label="Comparison of initial distributions")
    # plt.legend()
    # plt.show()
    
    # plt.figure()
    # plt.hist(multiplied_distribution, bins=100, label="multiplied distribution")
    # plt.hist(initial_distrib, bins=100, label="actual physics distribution",alpha=0.5)
    # plt.title(label="Overlapping distributions")
    # plt.legend()
    # plt.show()
    
    
    return sampling_distrib,weights
    
# weights = weighted_particles( gaussian(np.linspace(-10,10,100), mean=0, stddev=1), gaussian(np.linspace(-10,10,100), mean=3, stddev=1))

mean_initial = 100
stddev_initial = 10.0

mean_desired = 110
stddev_desired = 25.0

no_particles = 10000

initial_distrib = np.random.normal(mean_initial,stddev_initial,no_particles)
sampling_distrib = np.random.normal(mean_desired,stddev_desired,no_particles)

s_distrib,w = weighted_particles(initial_distrib,sampling_distrib)


# print(Levy(x=np.linspace(0.1,10,10000)))

# plt.figure()
# plt.hist(initial_distrib, bins=100, label="actual physics distribution")
# plt.hist(sampling_distrib, bins=100, label="sampling distribution",alpha=0.5)
# plt.title(label="Comparison of distributions")
# #plt.xlim(90,110)
# plt.legend()
# plt.show()


# plt.figure()
# plt.hist(initial_distrib, bins=100, label="actual physics distribution")
# plt.hist(s_distrib, bins=100, label="s distribution",alpha=0.5,weights=w)
# plt.title(label="of distributions")
# #plt.xlim(90,110)
# plt.legend()
# plt.show()



# plt.figure()
# plt.scatter(s_distrib,w)
# plt.show()