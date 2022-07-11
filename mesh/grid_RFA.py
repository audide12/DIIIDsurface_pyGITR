#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for micro-trenches experiments on DiMES
# Reference paper:
@author: guterl
"""
import gmsh
import math

# (Re)initialize gsmh session
try:
    gmsh.finalize()
except:
    pass
finally:
    gmsh.initialize()

# Define grid boundary
xGrid = 0
yGrid = 0
zGrid = 0

dxGrid = 6.5e-3
dyGrid = 5e-3
dzGrid = 50e-6

dxBlock = 0.2e-3
dyBlock = 0.04e-3

dxInter = 0.025e-3
dyInter = 0.025e-3

# Mesh parameters
LMesh= 0.00001

# Make grid geometry
x = xGrid
y = yGrid
i=0
j=0
TagRec=[]
while x <= xGrid+dxGrid:
    TagRec.append(gmsh.model.occ.addRectangle(x,yGrid,zGrid,dxInter,dyGrid+dyInter))
    x = x + dxBlock+dxInter
    i += 1
TagRec.append(gmsh.model.occ.addRectangle(xGrid+dxGrid,yGrid,zGrid,dxInter,dyGrid+dyInter))
while y <= yGrid+dyGrid:
    TagRec.append(gmsh.model.occ.addRectangle(xGrid+dxInter,y,zGrid,dxGrid,dyInter))
    y = y + dyBlock+dyInter
    j += 1
TagRec.append(gmsh.model.occ.addRectangle(xGrid,yGrid+dyGrid,zGrid,dxGrid+dxInter,dyInter))
TagTemp=TagRec.pop(0)
TagGrid = gmsh.model.occ.fuse([(2,TagTemp)],[(2,idx) for idx in TagRec])
#TagGrid =gmsh.model.occ.extrude(TagGrid[0], 0, 0, dzGrid)

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# SetupMesh
gmsh.option.setNumber("Mesh.MeshSizeMin", LMesh)
gmsh.option.setNumber("Mesh.MeshSizeMax", LMesh)

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
#gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("leding_edges_DiMES.msh")

# close gmsh session
#gmsh.finalize()

