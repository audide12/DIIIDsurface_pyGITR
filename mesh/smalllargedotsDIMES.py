#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for small/large dots experiments on DiMES
# Reference paper: https://iopscience.iop.org/article/10.1088/1361-6587/ab5144/meta
@author: Jerome Guterl (guterlj@fusion.gat.com)
"""

import gmsh

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

# Tags
TagDiMES0=0
TagDiMES=10
TaglDot=20
TagsDot=30

# Create 2D surfaces
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0)
gmsh.model.occ.addDisk(xsDot,ysDot,zsDot,rsDot,rsDot,TagsDot)
gmsh.model.occ.addDisk(xlDot,ylDot,zlDot,rlDot,rlDot,TaglDot)
gmsh.model.occ.cut([(2,TagDiMES0)],[(2,TagsDot),(2,TaglDot)],TagDiMES,removeTool=False)

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Set number of elements on the boundary of each dots and DiMES cap
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 50)

# Prevent very small elements in small dots
gmsh.option.setNumber("Mesh.MeshSizeMin", rsDot/3)

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("small_large_dots_DiMES.geo")

# Close gmsh session
gmsh.finalize()