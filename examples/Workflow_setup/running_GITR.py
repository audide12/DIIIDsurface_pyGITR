#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 14:14:20 2022

@author: de
"""

#%%

Run = pyGITR.Run()
Run.Verbose = True

Run.SetReferenceDirectory('.')
Run.SetSimRootPath('/Users/de/Research/DIIIDsurface_pyGITR/examples/Workflow_setup/')
Run.SetupScan(OverWrite=True)

Run.Clean()
Run.LaunchBatch()
