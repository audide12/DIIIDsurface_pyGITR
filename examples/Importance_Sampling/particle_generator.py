#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:41:13 2023

@author: de
"""


from pyGITR.particleSource_functions import *
from pyGITR.Particles import *

nP = 1000000

p = ParticleDistribution()
p.SetAttr('Np', nP)

#p.SetAttr(['vz'],'Levy', x=np.linspace(0.001,10,1000), c=2, mu=0)   # without importance sampling
    
mean = 2.0
sigma = 0.1


p.SetAttr_weighted(['vz'],'Thomson','Gaussian')  # first is p and second as q
p.SetAttr(['vx'],'Gaussian')
##  p is actual distribution and q is sampling distribution

q_distribution = p.Particles['vx'] # q distribution

sampled_distribution = p.Particles['vz']
weights = p.Particles['weights']
weights = nP*weights/np.sum(weights)

p_distribution = p.Particles['actual']      # p distribution

#%%

plt.figure()
#plt.hist(sampled_distribution, bins=100, weights=np.ones(len(sampled_distribution)),label="sampled with weights --")
plt.hist(sampled_distribution, bins=50, weights=weights,label="sampled with weights")
#plt.hist(sampled_distribution, bins=100,label="sampled without weights")
plt.hist(p_distribution,bins=50, label="physics",alpha=0.5)
plt.legend()
plt.show()

#%%
plt.figure()
plt.hist(q_distribution, bins=100,label="q distribution")
#plt.hist(sampled_distribution, bins=50, weights=weights,label="sampled with weights")
#plt.hist(sampled_distribution, bins=100,label="sampled without weights")
#plt.hist(p_distribution,bins=50, label="physics",alpha=0.5)
plt.legend()
plt.show()