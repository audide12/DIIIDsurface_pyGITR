#import electronvolt_num as units
import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
plt.ion()
FileNameSurface='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/surface.nc'
FileNameParticle='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/particleSource.nc'
import netCDF4
import os
#from SimManager import rget

class PostProcess():
    def __init__(self, Simulations = None, **kwargs):
        if type(Simulations) == str:
            Simulations = np.load(Simulations,allow_pickle=True)
        if type(Simulations) == np.ndarray:
            Simulations = Simulations.tolist()
        if type(Simulations) != list:
            Simulations = [Simulations]


        self.Simulations = Simulations
        self.SurfaceFile = 'output/surface.nc'
        self.ParticleStartFile = 'output/particleSource.nc'
        self.ParticleEndFile = 'output/particleSource_end.nc'
        self.HistoryFile = 'output/history.nc'
        self.SurfacePath = None
        self.Verbose = False
        for k,v in kwargs:
            if hasattr(self,k):
                setattr(self,k,v)
        print(self.Simulations)
        Idx = np.array(tuple(S.Idx for S in self.Simulations))

        self.NDim = Idx.shape[1]
        self.IdxDim = tuple(np.max(Idx,axis=0)+1)
        if len(self.IdxDim) == 1:
            self.IdxDim = (np.max(Idx)+1,1)
        print('IdxDim:',self.IdxDim)
        self.IdxArray = np.empty(self.IdxDim,dtype=int)
        for S in self.Simulations:
            self.IdxArray[S.Idx] = S.Id
        if len(self.Simulations)>0:
            assert all([S.Parameters == self.Simulations[0].Parameters for S in Simulations] )
            for S in self.Simulations:
                S.DicParameters = dict((k,v) for k,v in zip(S.Parameters,S.Values))
            
            self.Parameters = self.Simulations[0].Parameters
            self.DicParameters = dict((k,list(set([S.DicParameters[k] for S in self.Simulations]))) for k in self.Parameters)

            
        else:
            self.Parameters = []

        self.GetSurfaceData()
        self.GetParticleData()


    def GetSurfaceData(self, Format='NETCDF4'):
        for S in self.Simulations:
            FilePath = os.path.join(S.Directory,self.SurfaceFile)
            if self.Verbose:
                print('Importing surface data from {}'.format(FilePath))
            if not hasattr(S,'Data'):
                S.Data = {}
            S.Data['SurfaceData'] = ImportSurfaceData(FilePath, Format = Format)

    def GetHistoryData(self, Format='NETCDF4'):
        for S in self.Simulations:
            FilePath = os.path.join(S.Directory,self.HistoryFile)
            if self.Verbose:
                print('Importing history data from {}'.format(FilePath))
            if not hasattr(S,'Data'):
                S.Data = {}
            S.Data['HistoryData'] = ImportHistoryData(FilePath, Format = Format)

    def GetParticleData(self, Format='NETCDF4'):
        for S in self.Simulations:
            FilePath = os.path.join(S.Directory,self.ParticleStartFile)
            if self.Verbose:
                print('Importing particle data from {}'.format(FilePath))
            if not hasattr(S,'Data'):
                S.Data = {}
            S.Data['ParticleStartData'] = ImportParticleData(FilePath, Format = Format)

            FilePath = os.path.join(S.Directory,self.ParticleEndFile)
            if self.Verbose:
                print('Importing particle data from {}'.format(FilePath))
            if not hasattr(S,'Data'):
                S.Data = {}
            S.Data['ParticleEndData'] = ImportParticleData(FilePath, Format = Format)

    def PlotArray(self,GetData, **kwargs):
        self.fig, self.ax = plt.subplots(*self.IdxDim)
        if type(self.ax) != np.ndarray:
            self.ax = np.array([self.ax])
        if len(self.ax.shape)<2:
            self.ax = self.ax.reshape(*self.IdxDim)
        for j in range(self.IdxDim[1]):
            for i in range(self.IdxDim[0]):
                Sim = self.Simulations[self.IdxArray[i,j]]
                GetData(self.ax[i,j],Sim, **kwargs)
                Str = ' \n '.join(["{}={}".format(p,v) for (p,v) in zip(Sim.Parameters, Sim.Values)])
                self.ax[i,j].set_title(Str, loc='center', wrap=True)
                
    def CumulateParticleData(self):
        if not hasattr(self, 'CumulativeData'):
            self.CumulativeData = {}
        for p in ['ParticleStartData','ParticleEndData']:
            self.CumulativeData[p] = {'Data':dict((k,np.concatenate(tuple([S.Data[p]['Data'][k] for S in self.Simulations]))) for k in self.Simulations[0].Data[p]['Data'].keys())}

        # x = [Post.Simulations[k].Values[PlotAxis] for k in Post.IdxArray.take(j,axis=GroupAxis)]
        # p = [Post.Simulations[k].Values[GroupAxis] for k in Post.IdxArray.take(j,axis=GroupAxis)]
    # assert all([pp == p[0] for pp in p]) 'Something is wrong with the parameters'
    # y = rget(Sim,['SurfaceData','Data','E'])
    # plot()


def ImportSurfaceData(FilePath='output/surface.nc', Format='NETCDF4'):
    SurfaceData = {}
    SurfaceData['FilePath'] = FilePath
    SurfaceData['RawData'] = LoadData( FilePath, Format)
    SurfaceData['Data'] = GetDistribSurface(SurfaceData['RawData'])
    return SurfaceData


def ImportHistoryData(FilePath='output/history.nc', Format='NETCDF4'):
    Data = {}
    Data['FilePath'] = FilePath
    Data['RawData'] = LoadData( FilePath, Format)
    Data['Data'] = GetHistoryParticle(Data['RawData'])
    return Data

def ImportParticleData(FilePath='output/particleSource_end.nc', Format='NETCDF4'):
    ParticleData = {}
    ParticleData['FilePath'] = FilePath
    ParticleData['RawData'] = LoadData( FilePath, Format)
    ParticleData['Data'] = GetDistribParticle(ParticleData['RawData'])
    return ParticleData

def LoadData(FilePath, Format='NETCDF4'):
        FilePath = os.path.expanduser(FilePath)
        if Format == 'NETCDF4':
            return Dataset(FilePath, "r", format="NETCDF4")
        else:
            raise ValueError('Unknown format.')

def GetDistribSurface(SurfaceData):
    Data = {}
    Data['Distrib'] = np.array(SurfaceData.variables['surfEDist'])
    Data['GridE'] = np.array(SurfaceData.variables['gridE'])
    Data['GridA'] = np.array(SurfaceData.variables['gridA'])
    Data['E'] = np.sum(Data['Distrib'][:,:,:],(0,2))
    Data['A'] = np.sum(Data['Distrib'][:,:,:],(0,1))
    Data['S'] = np.sum(Data['Distrib'][:,:,:],(1,2))
    Data['Tot'] = np.sum(Data['Distrib'][:,:,:],(0,1,2))
    return Data

def GetDistribParticle(ParticleData):
    Data={}
    for k in ParticleData.variables.keys():
        Data[k] = np.array(ParticleData.variables.get(k))
    v = np.sqrt(Data['vz']**2+Data['vy']**2+Data['vx']**2)
    vp = np.sqrt(Data['vy']**2+Data['vx']**2)
    Data['E'] = 1/2*Data['amu']*(v**2)*units.mp/units.eV*(units.V)**2
    Data['theta'] = np.arccos(np.abs(Data['vz'])/v)
    Data['phi'] = np.arccos(np.abs(Data['vx'])/vp)
    return Data

def GetHistoryParticle(ParticleData):
    Data={}
    for k in ParticleData.variables.keys():
        Data[k] = np.array(ParticleData.variables.get(k))
    return Data

# class PostProcess():
#     def __init__(self):
#         pass
#     def LoadData(self,FilePath, Format='NETCDF4'):
#         if Format == 'NETCDF4':
#             return Dataset(FilePath, "r", format="NETCDF4")
#         else:
#             raise ValueError('Unknown format.')
#     def ImportData(self,FileName,Folder='output',Path=None):
#         if Path is None:


# class DistributionAnalysis()

# SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
# ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")
# D={}
# vx = np.array(ParticleData.variables['vx'])
# vy = np.array(ParticleData.variables['vy'])
# vz = np.array(ParticleData.variables['vz'])
# amu = np.array(ParticleData.variables['amu'])
# E0=1/2*amu*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2

# Distrib = np.array(SurfaceData.variables['surfEDist'])
# GridE = np.array(SurfaceData.variables['gridE'])
# GridA = np.array(SurfaceData.variables['gridA'])
# E=np.sum(Distrib[:,:,:],(0,2))
# A=np.sum(Distrib[:,:,:],(0,1))
# S=np.sum(Distrib[:,:,:],(1,2))
# print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))
# fig = plt.figure()
# ax = fig.subplots(1,4)
# ax[0].plot(GridE,E)
# ax[0].hist(E0,GridE)
# ax[1].plot(GridA,A)
# ax[2].hist(vx)
# ax[2].hist(vy,alpha=0.2)
# ax[2].hist(vz,alpha=0.2)
# ax[3].plot(np.arange(0,S.shape[0]),S)
# #%%
# FileNameHistory='/home/jguterl/Dropbox/python/pyGITR/examples/large_box4/output/history.nc'
# HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")
# x = np.array(HistoryData.variables['x'])
# z = np.array(HistoryData.variables['z'])
# y = np.array(HistoryData.variables['y'])

# nT = HistoryData.dimensions['nT'].size
# nP = HistoryData.dimensions['nP'].size
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for i in range(0,nP):
#     ax.plot(x[i,:],y[i,:],z[i,:])
def PlotEStart(ax,Sim):
    E= rget(Sim.Data,['ParticleStartData','Data','E'])
    ax.hist(E)
    ax.set_xlabel('E[eV]')

def PlotEEnd(ax,Sim):
    E= rget(Sim.Data,['ParticleStartData','Data','E'])
    ax.hist(E)
    ax.set_xlabel('E[eV]')

def PlotEStartEnd(ax,Sim,**kwargs):
    Estart= rget(Sim.Data,['ParticleStartData','Data','E'])
    Eend= rget(Sim.Data,['ParticleEndData','Data','E'])
    ax.hist(Estart,**kwargs)
    ax.hist(Eend,**kwargs)
    ax.set_xlabel('E[eV]')

def PlotEStartEnd(ax,Sim,**kwargs):
    Estart = rget(Sim.Data,['ParticleStartData','Data','E'])
    Eend= rget(Sim.Data,['ParticleEndData','Data','E'])
    ax.hist(Estart,**kwargs)
    ax.hist(Eend,**kwargs)
    ax.set_xlabel('E[eV]')

def SurfaceAngle(ax,Sim):
    A = rget(Sim.Data,['SurfaceData','Data','A'])
    gridA = rget(Sim.Data,['SurfaceData','Data','GridA'])
    ax.plot(gridA,A)
    ax.set_xlabel('angle [degree]')

def PlotHistory(Sim):
    x = Sim.Data['HistoryData']['Data']['x']
    y = Sim.Data['HistoryData']['Data']['y']
    z = Sim.Data['HistoryData']['Data']['z']
    #fig = plt.figure()
    ax = plt.axes(projection='3d')
    for i in range(x.shape[0]):
        ax.plot(x[i,:],y[i,:],z[i,:])
        ax.set_xlabel('x')
        ax.set_ylabel('y')