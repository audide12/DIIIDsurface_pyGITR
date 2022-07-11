#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 13:42:31 2021

@author: jguterl
"""

for S in Data.Simulations.Parameters:
    Data.Simulations[0].Parameters[0]




Parameters = [k.split('.')[-1]for k in Run.ParamScan.keys()]
Nvalues = [len(v['Values']) for v in Run.ParamScan.values()]

Parameters = Post.Simulations[0].Parameters
GroupBy='biasPotential'
GroupAxis = Parameters.index(GroupBy)
PlotBy='lambda'
PlotAxis = Parameters.index(PlotBy)
GroupBy='biasPotential'
GroupAxis = Parameters[0].index(GroupBy)
PlotBy='lambda'
PlotAxis = Parameters.index(PlotBy)

for idx in np.ndindex(Post.IdxArray.shape):
    ValueParameter[idx] = Post.Simulations[Post.IdxArray[idx]].Values[]
    Data[idx] =

def GetData(ax,Sim):
            E= rget(Sim,['SurfaceData','Data','E'])
            gridE= rget(Sim,['SurfaceData','Data','E'])
            ax.plot(gridE,E)
            ax.set_xlabel('E[eV]')

def PlotArray(self,GetData):
    GroupBy='biasPotential'
GroupAxis = Parameters[0].index(GroupBy)
PlotBy='lambda'
PlotAxis = Parameters.index(PlotBy)

fig, ax = plt.subplots(*Post.IdxDim)
for j in range(Post.IdxDim[GroupAxis]):
    for i in range(Post.IdxDim[GroupPlot]):
        ax[i,j].plot(GetData(self.))
    x = [Post.Simulations[k].Values[PlotAxis] for k in Post.IdxArray.take(j,axis=GroupAxis)]
    p = [Post.Simulations[k].Values[GroupAxis] for k in Post.IdxArray.take(j,axis=GroupAxis)]
    assert all([pp == p[0] for pp in p]) 'Something is wrong with the parameters'
    y = rget(Sim,['SurfaceData','Data','E'])


def Data(Sim):
    E= rget(Sim,['SurfaceData','Data','E'])
    gridE= rget(Sim,['SurfaceData','Data','E'])
plt.plot()

fig, ax = plt.subplots(5,2)
