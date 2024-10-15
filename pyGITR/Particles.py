# -*- coding: utf-8 -*-


import numpy as np
from pyGITR.math_helper import *
from typing import Callable
import matplotlib.pyplot as plt
import pydoc
import netCDF4
import os

class Distribs():

    #def Gaussian(x: np.ndarray = np.linspace(-10, 10, 10000), sigma: float = 1.0, mu: float = 0.0, beta: float = 0.0, Normalized=True):
    #    f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
    #    #if beta > 0:
    #    #    f[np.argwhere(x<0)] = 0
    #    if Normalized:
    #        f = f/Integrale(f, x, Array=False)
    #    return f

    def Gaussian(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 30.0, mu: float = 70.0, Normalized=True):
        f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
        #f[np.argwhere(x<0)] = 0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f

    def Gaussian_test(x: np.ndarray = np.linspace(-15000, 15000, 100000), sigma: float = 15.0, mu: float = 30.0, Normalized=True):
        f = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)/sigma)**2)
        #f[np.argwhere(x<0)] = 0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f


    def Gaussian_Jerome(x: np.ndarray = np.linspace(0, 15000, 100000), sigma: float = 20.0, mu: float = 20, beta: float = 0.0, Normalized=True):
        f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
        if beta > 0:
            f[np.argwhere(x<0)] = 0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f    
    
    def Gaussian_Zack(x: np.ndarray = np.linspace(-10, 10, 10000), sigma: float = 1.0, mu: float = 0.0, beta: float = 0.0, Normalized=True):
        f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
        if beta > 0:
            f[np.argwhere(x<0)] = 0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f
    
    # def Maxwellian(x: np.ndarray, sigma: float = 1, mu: float = 0, beta: float = 0, Normalized=True):
    #     f = np.abs(x)**beta*np.exp(-1.0/2.0*((x-mu)/sigma)**2)
    #     if Normalized:
    #         f = f/Integrale(f, x, Array=False)
    #     return f    
    

    def Thomson(x: np.ndarray = np.linspace(0, 300, 10000), xb: float = 6.67, xc: float = 200, Normalized=True):
        #xb = 6.67 for C, xb = 8.64 for W, xb = 4.73 for Si
        
        assert not (xc <= xb), "xc cannot be <= xb"
        f = x/(x + xb) ** 3*(1.0-np.sqrt((x+xb)/(xc+xb)))
        #f[np.argwhere(x > xc)] = 0.0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f

    def Uniform(x=np.linspace(0, 1, 10), xmin=None, xmax=None, Normalized=True):
        assert not (xmin is not None and xmax is not None and xmin >= xmax), "xmin cannot be <= xmax"
        f = np.full(x.shape, 1)
        if xmin is not None:
            f[np.argwhere(x < xmin)] = 0
        if xmax is not None:
            f[np.argwhere(x > xmax)] = 0

        if Normalized:
            f = f/Integrale(f, x, Array=False)

        return f

    def SinCos(x=np.linspace(0, np.pi/2, 10000), xmin=None, xmax=None, Normalized=True):
        assert not (xmin is not None and xmax is not None and xmin >= xmax), "xmin cannot be <= xmax"

        f = np.sin(x)*np.cos(x)
        if xmin is not None:
            f[np.argwhere(x < xmin)] = 0
        if xmax is not None:
            f[np.argwhere(x > xmax)] = 0
        if Normalized:
            f = f/Integrale(f, x, Array=False)
        return f
    
    def Levy(x=np.linspace(0.1,10,10000), c=1, mu=0):
        return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)

    @classmethod
    def GetPdf(cls, f:str or Callable[np.array, np.ndarray], **kwargs) -> (np.ndarray,np.ndarray):

        if type(f) == str:
            if hasattr(cls, f):
                f = getattr(cls,f)
                if kwargs.get('x') is None:
                    x = f.__defaults__[0]
                else:
                    x = kwargs.pop('x')

            else:
                raise KeywordError('Cannot find the function "{}" in attributes.'.format(f))
        else:
            if kwargs.get('x') is None:
                raise KeywordError('Must provide a vector of x values for the function f to get f(x)')
            else:
                x = kwargs.pop('x')

        return x, f(x,**kwargs)
    
    @classmethod
    def GetWeightedPdf(cls, f_q:str or Callable[np.array, np.ndarray], f_p:str or Callable[np.array, np.ndarray], **kwargs) -> (np.ndarray,np.ndarray,np.ndarray):
        
        if type(f_q) == str:
            if hasattr(cls, f_q):
                f_q = getattr(cls,f_q)
                if kwargs.get('x') is None:
                    x = np.linspace(0, 5000, 100000)   #specifically for vz
                else:
                    x = kwargs.pop('x')
            else:
                raise KeywordError('Cannot find the function "{}" in attributes.'.format(f_q))
        else:
            if kwargs.get('x') is None:
                raise KeywordError('Must provide a vector of x values for the function f to get f(x)')
            else:
                x = kwargs.pop('x')
                
        if type(f_p) == str:
            if hasattr(cls, f_p):
                f_p = getattr(cls,f_p)

            else:
                raise KeywordError('Cannot find the function "{}" in attributes.'.format(f_p))
            
        # print(len(x))
        # print(f_p(x,**kwargs))
        # print(f_q(x,sigma,mu,**kwargs))

        return x, f_q(x,**kwargs), f_p(x,**kwargs)
    
    
    
    @classmethod
    def PlotDistrib(cls, DistribName, **kwargs):
        """


        Args:
            cls (TYPE): DESCRIPTION.
            DistribName (TYPE): DESCRIPTION.
            **kwargs (TYPE): DESCRIPTION.

        Returns:
            None.

        """
        if hasattr(cls, DistribName):
            f = getattr(cls,DistribName)
            if kwargs.get('x') is None:
                kwargs['x'] = f.__defaults__[0]
                Name = ';' + DistribName+(';').join(['{}={}'.format(v,d) for v,d in zip(f.__code__.co_varnames[1:-2],f.__defaults__[1:])])

            plt.plot(kwargs['x'], f(**kwargs), label=Name)
            plt.legend()
        else:
            print('Cannot find the distribution function "{}" in BaseDistrib.'.format(DistribName))

    @classmethod
    def ShowAvailablePdfs(cls):
        for D in ['Gaussian', 'Uniform', 'SinCos', 'Thomson']:
            print('Distribution: {}'.format(pydoc.render_doc(getattr(cls,D), "%s")))

    def GenerateDistribution(x:np.ndarray ,pdf:np.ndarray, N:int=10000):
        assert type(N) == int and N>0, "Nsamples must be an integer > 0. N={}".format(N)
        assert type(x) == np.ndarray and type(pdf) == np.ndarray, " x and pdf must be a numpy array"
        assert x.shape == pdf.shape, " x and pdf must have the same shape.:x.shape={}; pdf.shape={}".format(x.shape,pdf,shape)
        
        #print(len(x[np.searchsorted(Integrale(pdf, x, Normalized=True), np.random.rand(N), side='left')]))
        return x[np.searchsorted(Integrale(pdf, x, Normalized=True), np.random.rand(N), side='left')]
    

        
    
    def GenerateDistributionWeighted(x:np.ndarray ,pdf:np.ndarray, pdf_p:np.ndarray,pdf_q:np.ndarray, N:int=10000):
        assert type(N) == int and N>0, "Nsamples must be an integer > 0. N={}".format(N)
        assert type(x) == np.ndarray and type(pdf) == np.ndarray, " x and pdf must be a numpy array"
        assert x.shape == pdf.shape, " x and pdf must have the same shape.:x.shape={}; pdf.shape={}".format(x.shape,pdf,shape)
        
        rng = np.random.default_rng()
        
        xi = rng.random(N)
        
        sampled_x_physics = x[np.searchsorted(Integrale(pdf_p, x, Normalized=True), xi, side='left')]
        sampled_x = x[np.searchsorted(Integrale(pdf, x, Normalized=True), xi, side='left')]

        weights = 1/pdf_q[np.searchsorted(Integrale(pdf, x, Normalized=True), xi, side='left')]
        
        # plt.figure()
        # plt.hist(sampled_x_physics,density=True,bins=100,alpha=0.5,label="physics")
        # plt.hist(sampled_x ,density=True,bins=100, weights=weights,label="weighted samples")
        # plt.legend()
        # plt.show()
        
        return sampled_x,weights,sampled_x_physics
            
    
    

class ParticleDistribution():

    def __init__(self, Np=1000, ListAttr=['x','y','z','vx','vy','vz']):
        self.ListAttr = ListAttr
        self.Particles = {'Np':Np}
        for Attr in self.ListAttr:
            self.AddAttr(Attr)


    def ShowAvailablePdfs(self):
        Distribs.ShowAvailablePdfs()

    def PlotAttr(self, Attr=[], NewFig=True):
        if type(Attr) !=list:
            Attr = [Attr]
        if Attr == []:
            Attr = [A for A,V in self.Particles.items() if type(V) == np.ndarray]
        if NewFig:
            plt.figure()
        for A in Attr:
            assert A in list(self.Particles.keys()), 'Cannot find "{}" among list of attributes available:{}'.format(A, list(self.Particles.keys()))
            v = self.Particles[A]

            plt.hist(v, label=A, density=True)
            plt.legend()




    def AddAttr(self, Attr):
        assert self.Particles.get(Attr) is None, 'Attribute {} already exists.'.format(Attr)
        self.Particles[Attr] = np.full((self.Particles['Np']),0)
        if Attr not in self.ListAttr:
            self.ListAttr.append(Attr)

    def ShowAttr(self):
        print('--- Particle attributes ---')
        for A,V in self.Particles.items():
            Type = V.__class__.__name__
            if type(V) == np.ndarray:
                Type = "{}:{}".format(Type,V.shape)
            print('-- {} = {} ({})'.format(A,V, Type))
        print('--- End of Particle attributes ---')

    def SetAttr(self, Attr: list or str, f: float or np.ndarray or str or int, **kwargs):
        if type(Attr) == str:
            Attr = [Attr]
        for A in Attr:
            assert A in list(self.Particles.keys()), 'Cannot find "{}" among list of attributes available:{}'.format(A, list(self.Particles.keys()))
            if type(f) == float or type(f) == int:
                if A != 'Np':
                    self.Particles[A] = np.full((self.Particles['Np']),f)
                else:
                    self.Particles[A] = f

            elif type(f) == np.ndarray:
                self.Particles[A] = f
            else:
                self.Particles[A] = self.Generate(self.Particles['Np'], f, **kwargs)
                
    def SetAttr_weighted(self, Attr: list or str, f_p: str, f_q: str, mu: float = 5, sigma: float = 1.0  , **kwargs):
        if type(Attr) == str:
            Attr = [Attr]
            
            
        for A in Attr:
            assert A in list(self.Particles.keys()), 'Cannot find "{}" among list of attributes available:{}'.format(A, list(self.Particles.keys()))
            self.Particles[A],self.Particles['weights'],self.Particles['actual'] = self.Generate_Weighted(self.Particles['Np'], f_q, f_p , **kwargs)        
            
            # self.Particles[A],self.Particles['weights'], self.Particles['actual'] = self.Generate_Weighted(self.Particles['Np'], f_q, f_p ,mu, sigma, **kwargs)        
        
                        

    def ScaleAttr(self,Attr, ScaleFactor):
        if type(Attr) == str:
            Attr = [Attr]
        for A in Attr:
            assert A in list(self.Particles.keys()), 'Cannot find "{}" among list of attributes available:{}'.format(A, list(self.Particles.keys()))
            self.Particles[A] = self.Particles[A]*ScaleFactor

    def RotateAngle(self, Field:str, theta:float, phi:float, Degree=True):
        self.Rotate(Field, [0,1,0], theta, Degree)
        self.Rotate(Field, [0,0,1], phi, Degree)

    def Rotate(self, Field:str, AxisVector:np.ndarray or list, Angle: float, Degree=True) -> None:
        """
        Rotate the vector field (x,y,z) or (vx,vy,vz) around the axis AxisVector by angle Angle

        Args:
            Field (str): DESCRIPTION.
            AxisVector (np.ndarray or list): DESCRIPTION.
            Angle (float): DESCRIPTION.
            Degree (TYPE, optional): if True, angle are in degree. Defaults to True.

        Returns:
            None: DESCRIPTION.

        """
        assert Field == 'v' or Field == 'x', 'Field must be either "v" or "x"'
        assert type(AxisVector) == np.ndarray or list, 'AxisVector must be a list or numpy array'

        if Field == 'x':
            Fields =  ['x','y','z']
        else:
            Fields = ['vx','vy','vz']

        Vs = [self.Particles[k] for k in  Fields]
        V  = RotateCoordinates(*Vs, AxisVector, Angle, Degree)
        for i,F in enumerate(Fields):
            self.Particles[F] = V[:,i].squeeze()

    def Generate(self, Np, DistribName, **kwargs):
        """


        Args:
            Np (int, optional): Number fo samples. Defaults to 10000.
            DistribName (TYPE): Probability distribution function.

        Returns:
            Distribution (np.ndarray): distribution

        """

        x, pdf = Distribs.GetPdf(DistribName, **kwargs)
        
        return Distribs.GenerateDistribution( x, pdf, Np)
    
    
    def Generate_Weighted(self, Np, q_DistribName, p_DistribName, **kwargs):
        """
        Args:
            Np (int, optional): Number fo samples. Defaults to 10000.
            DistribName (TYPE): Probability distribution function.

        Returns:
            Distribution (np.ndarray): distribution

        """
        
        x, pdf_q, pdf_p = Distribs.GetWeightedPdf(q_DistribName,p_DistribName, **kwargs)
        
        pdf_p_times_q = np.multiply(pdf_p,pdf_q)
        
        Normalization = Integrale(pdf_p_times_q, x, Array=False)
        #Normalization = 1.0
        
        pdf_p_times_q = np.divide(pdf_p_times_q,Normalization)
        
        sampled_x,weights,sampled_x_physics = Distribs.GenerateDistributionWeighted(x, pdf_p_times_q,pdf_p, pdf_q, Np)
        
        #weights = np.divide(1.0,pdf_q[np.searchsorted(x,sampled_x)])
                
        return sampled_x, weights,sampled_x_physics

    def CheckParticles(self):
        for k,v in self.Particles.items():
            if k != 'Np' and k in self.ListAttr and type(v) == np.ndarray:
                assert len(v.shape)<2 and v.shape[0] == self.Particles['Np'], 'Wrong dimension for particle attribute {}  with shape {}. Np={}'.format(k,v.shape,self.Particles['Np'])

    def WriteParticleFile(self,FileName, Format='f8', Folder='input'):
        self.CheckParticles()
        File = netCDF4.Dataset(os.path.join(Folder,FileName), 'w', format='NETCDF4')

        np = File.createDimension('nP', self.Particles['Np'])
        Var = {}
        for k in self.ListAttr:
            Var[k] = File.createVariable(k,'f8','nP')
            Var[k][:] = self.Particles[k]
        File.close()










