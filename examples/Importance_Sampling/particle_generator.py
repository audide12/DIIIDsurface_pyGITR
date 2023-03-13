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


p.SetAttr_weighted(['vz'],'Gaussian','Gaussian_test')  # first is p and second as q

p.SetAttr(['vx'],'Gaussian')

##  p is actual distribution and q is sampling distribution

actual_distribution = p.Particles['actual']

actual_distribution_1 = p.Particles['vx']



sampled_distribution = p.Particles['vz']
weights = p.Particles['weights']


multiplied = np.multiply(sampled_distribution,weights)

plt.figure()
plt.hist(actual_distribution, bins=100, label="actual physics distribution")
plt.hist(sampled_distribution, bins=100, label="sampling physics distribution",alpha=0.5)
plt.title(label="Comparison of distributions")
#plt.xlim(90,110)
plt.legend()
plt.show()




plt.figure()
plt.hist(multiplied, bins=100, label="Sampling distribution multiplied by weights")
# plt.hist(sampled_distribution, bins=100, label="sampling distribution")
# plt.hist(weights, bins = 100, label = "weights")
plt.title(label="Multiplied")
# plt.xlim(0,10)
plt.legend()
plt.show()



# plt.figure()
# plt.scatter(sampled_distribution,weights)
# plt.xlabel("Sampled Distributions")
# plt.ylabel("Weights")
# plt.title(label="Sampled Distribution vs weights")
# plt.show()



#p.PlotAttr('vx')

#p.WriteParticleFile(particleFile1)


