#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for small/large dots experiments on DiMES
# Reference paper: https://iopscience.iop.org/article/10.1088/1361-6587/ab5144/meta
@author: Jerome Guterl (guterlj@fusion.gat.com)
"""

import gmsh
from pyGITR.gmsh_helper import SetGroups

# Initialize gmsh session
gmsh.initialize()

# Define DiMES cap surface
xDiMES=0.0
yDiMES=0.0
zDiMES=0.0
rDiMES=0.025

# Define small dot on DiMES cap surface
xsDot=-0.01
ysDot=0.0
zsDot=0.0
rsDot=0.0005

# Define large dot on DiMES cap surface
xlDot=0.0
ylDot=0.0
zlDot=0.0
rlDot=0.005

# Define metal rings
xBlock=xDiMES
yBlock=yDiMES
zBlock=zDiMES
dxBlock = 0.01
dyBlock = 0.01

# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
dxBox = rDiMES*12.5
dyBox = rDiMES*12.5
dzBox = rDiMES*12.5

# Tags
TagDiMES0=0
TagDiMES=10
TaglDot=20
TagsDot=30
TagBox0=40

# Create metal ring
gmsh.model.occ.addRectangle(xBlock, yBlock, zBlock, dxBlock, dyBlock)

# Create 2D surfaces
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0) # the second r is used for elipse definition

# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
TagBox0 = gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0)
# print(type(s0))
# print(s0)
gmsh.model.occ.remove([(3, TagBox0)]) # removes the solid volume?
gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0)], -1, removeTool=False)
s1 = gmsh.model.occ.getEntities(2)
TagBoundBox = list(set(s1) - set(s0))
# print(type(s1))
# print(s1)

# Create dots
gmsh.model.occ.addDisk(xsDot,ysDot,zsDot,rsDot,rsDot,TagsDot)
gmsh.model.occ.addDisk(xlDot,ylDot,zlDot,rlDot,rlDot,TaglDot)
DiMESTag = gmsh.model.occ.cut( [(2,TagDiMES0)],[(2,TagsDot),(2,TaglDot)], TagDiMES,removeTool=False) # cuts the end from the first in the list, assigns new tag

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Set number of elements on the boundary of each dots and DiMES cap
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 50)



# Prevent very small elements in small dots
gmsh.option.setNumber("Mesh.MeshSizeMin", rsDot/3)

# Define groups to allow setting properties of elements when generating geometry input with GeomSetup
TaglDot = [(2,20)]
TagsDot = [(2,30)]
SetGroups(gmsh.model, DiMESTag, "DiMES", [255, 0, 0])
SetGroups(gmsh.model, TaglDot, "Large Dot", [0, 255, 0])
SetGroups(gmsh.model, TagsDot, "Small Dot", [0, 0, 255])
SetGroups(gmsh.model, TagBoundBox, "BoundBox", [230, 230, 250])

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("small_large_dots_DiMES.msh")

# Close gmsh session
gmsh.finalize()