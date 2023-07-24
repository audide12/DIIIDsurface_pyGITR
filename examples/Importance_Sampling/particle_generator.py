#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:41:13 2023

@author: de
"""


from pyGITR.particleSource_functions import *
from pyGITR.Particles import *

nP = 100000

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

p_distribution = p.Particles['actual']      # p distribution


plt.figure()
plt.hist(p.Particles['actual'],density=True,bins=100,alpha=0.5,label="Distribution")
plt.hist(p.Particles['vz'] ,density=True,bins=100, weights=p.Particles['weights'],label="weighted samples")
plt.legend()
plt.show()

#%%
#weights1 = nP*weights/(np.sum(weights))
weights1 = weights/(np.sum(weights))


plt.figure()
plt.hist(p_distribution,bins=100,label="physics")
plt.hist(sampled_distribution,bins=100, weights=weights1,label="weighted samples",alpha=0.5)
plt.legend()
plt.show()

#%%
plt.figure()
plt.hist(q_distribution, bins=100,label="Sampling distribution")
#plt.hist(sampled_distribution, bins=50, weights=weights,label="sampled with weights")
#plt.hist(sampled_distribution, bins=100,label="sampled without weights")
plt.hist(p_distribution,bins=50, label="Underlying physics distribution",alpha=0.3)
plt.legend()
plt.show()