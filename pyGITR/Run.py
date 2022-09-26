#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 20:37:20 2021

@author: jguterl
"""


import libconf
import os,io
import shutil
from SimManager import SimulationManager, MakeSimFolder, MakeParameterArray, MakeParamInfo, UpdateInputFile



class Run(SimulationManager):
    def __init__(self):
        self.GITRExecPath = os.path.expanduser('/Users/de/Research/GITR/build/GITR')
        self.ParamScan = {}
        self.ParamModif = {}
        self.ReferenceDirectory = None
        self.SimRootPath = None
        self.Verbose = False
        self.ListDirectory = []
        self.Init()

    def ShowSummary(self, ShowAll=False):
        print('GITRRun instance attributes:')
        for k in list(self.__dict__.keys()):
            if k == 'ListScan' and not ShowAll:
                continue
            else:
                print(' - {} : {}'.format(k,getattr(self,k)))

    def SetGITRExec(self, ExecPath:str):
        ExecPath = os.path.expanduser(ExecPath)
        assert os.path.exists(ExecPath) == True, 'Cannot find GITR executable: {}'.format(ExecPath)
        self.GITRExecPath = ExecPath



    def SetParamScan(self,  ConfigFile,Params, Values=None):
        self.ParamScan = {}
        if type(Params) == dict:
            self.ParamScan=dict((k,{'ConfigFile':ConfigFile,'Values':v}) for k,v in Params.items())
        else:
            if type(Values) != list:
                raise KeyError('Values of parameters must be given as a a list')
            self.ParamScan[Params] = {'ConfigFile':ConfigFile,'Values':Values}



    def AddParamScan(self, ConfigFile, Params, Values=None):
        if type(Params) == dict:
            self.ParamScan.update(dict((k,{'ConfigFile':ConfigFile,'Values':v}) for k,v in Params.items()))
        else:
            if type(Values) != list:
                raise KeyError('Values of parameters must be given as a a list')
            self.ParamScan[Params] = {'ConfigFile':ConfigFile,'Values':Values}


    def SetModifParam(self, ConfigFile, Params, Value=None, AddParam=False):
        self.ParamModif = {}
        if type(Params) == dict:
            self.ParamModif= dict((k,{'Value':Value,'ConfigFile':ConfigFile,'AddParam':AddParam}) for k,v in Params.items())
        else:
            assert Value is not None
            self.ParamModif[Params]={'Value':Value,'ConfigFile':ConfigFile,'AddParam':AddParam}

    def ModifParam(self, ConfigFile, Params, Value=None):

        if type(Params) == dict:
            self.ParamModif.update(dict((k,{'Value':Value,'ConfigFile':ConfigFile}) for k,v in Params.items()))
        else:
            assert Value is not None
            self.ParamModif[Params]={'Value':Value,'ConfigFile':ConfigFile}

    def ApplyModifyParam(self, AddParam=False):
        Dic={}
        Dic['ParameterInfo'] = MakeParamInfo(self.ParamModif)
        for D in self.ListDirectory:
            Dic.update({'Directory': D})
            UpdateInputFile(Dic, self.LoadMethod, self.DumpMethod, AddParam, self.Verbose)










    def SetReferenceDirectory(self, ReferenceDirectory:str)->None:
        ReferenceDirectory = os.path.expanduser(ReferenceDirectory)
        assert os.path.exists(ReferenceDirectory) == True, 'Cannot find ReferenceDirectory: {}'.format(ReferenceDirectory)
        self.ReferenceDirectory = ReferenceDirectory

    def SetSimRootPath(self, SimRootPath:str)->None:
        SimRootPath = os.path.expanduser(SimRootPath)
        self.SimRootPath = SimRootPath
    @staticmethod
    def DumpMethod(FilePath:str, Config:dict)->None:
            FilePath =  os.path.expanduser(FilePath)
            assert os.path.exists(FilePath), "Cannot find the configfile {}".format(FilePath)
            f = io.open(FilePath,'w')
            libconf.dump(Config,f)
            f.close()
    @staticmethod
    def LoadMethod(FilePath:str)->dict:
            FilePath =  os.path.expanduser(FilePath)
            assert os.path.exists(FilePath), "Cannot find the configfile {}".format(FilePath)
            f = io.open(FilePath,'r')
            Config = libconf.load(f)
            f.close()
            return Config

    def SetupScan(self, OverWrite=False, Format='value', AddParam=False):
        assert type(self.ReferenceDirectory)==str and os.path.exists(self.ReferenceDirectory), "Cannot find base folder path {}".format(self.ReferenceDirectory)
        assert type(self.SimRootPath)==str
        assert os.path.exists(self.GITRExecPath)
        print('>>>> Setting up scanning of parameters {} from reference folder {} with SimRootPath: {}'.format(list(self.ParamScan.keys()),self.ReferenceDirectory,self.SimRootPath))
        print('Executable: {}'.format(self.GITRExecPath))
        self.ParameterArray = MakeParameterArray(self.ParamScan, Format, self.Verbose)

        for SimInfo in self.ParameterArray.flat:
            Directory = MakeSimFolder(SimInfo['Suffix'], self.ReferenceDirectory, self.SimRootPath, OverWrite, self.Verbose)
            self.SetSimulation(SimInfo,Directory, self.GITRExecPath, self.LoadMethod, self.DumpMethod, AddParam, self.Verbose)
            self.ListDirectory.append(Directory)
        print('>>>> Applying modification of fixed parameters')
        self.ApplyModifyParam()
        print('>>>> Modification of fixed parameters applied.')
        self.DumpInfo(Folder=self.SimRootPath)



    def Clean(self,OutputDirectory='output'):
        for L in self.ListDirectory:
            if os.path.exists(os.path.join(L,OutputDirectory)):
                print('Cleaning {}'.format(os.path.join(L,OutputDirectory)))
                shutil.rmtree(os.path.join(L,OutputDirectory))
            os.mkdir(os.path.join(L,OutputDirectory))











        #         ListNewFolder =["{}_{}_{}".format(self.BaseFolder,ParamName, Value)
        #     for k,v in self.ParamScan.items():
        #     ParamName = k.split('.')[-1]
        #     for a in itertools.product(v1,v2,v3):
        #     for V in v:
        #         NewFolder = "{}_{}_{}".format(self.BaseFolder,ParamName, Value)



        # NewFolder = "{}_{}_{}".format(BaseFolder,'_'.join(ListSuffix))
        # print("Copying {} to {}".format(BaseFolder,NewFolder))
        # CopyFolder(BaseFolder,NewFolder)



