#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:53:18 2021

@author: jguterl
"""
import numpy as np
FileName='/home/jguterl/Dropbox/python/pyGITR/examples/large_box/output/particleSource.nc'
#FileName='/home/jguterl/Dropbox/python/pyGITR/examples/large_box/input/particleConf.nc'
import netCDF4
from netCDF4 import Dataset
rootgrp = Dataset(FileName, "r", format="NETCDF4")
rootgrp.variables['vz'][:]

FileName='/home/jguterl/Dropbox/python/pyGITR/examples/large_box/output/surface.nc'
SurfaceData = Dataset(FileName, "r", format="NETCDF4")
Distrib=np.array(SurfaceData.variables['surfEDist'])