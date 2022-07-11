#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for micro-trenches experiments on DiMES
# Reference paper:
@author: guterl
"""
import gmsh
from pyGITR.gmsh_helper import SetGroups

# Clean up gmsh session if script already executed
# (Re)initialize gsmh session
try:
    gmsh.finalize()
except:
    pass
finally:
    gmsh.initialize()

# Define DiMES cap surface
xDiMES = 0.0
yDiMES = 0.0
zDiMES = 0.0
rDiMES = 0.025

# Define 4 identical blocks on center of DiMES cap surface
NBlocks = 4
xBlock = [-0.01, 0.0, 0.0, -0.01]
yBlock = [-0.01, -0.01, 0.0, 0.0]
zBlock = [0.0, 0.0, 0.0, 0.0]
zBlock = [0.0, 0.0, 0.0, 0.0]
dzBlock = [0.001, 0.0003, 0.0, 0.0]
dxBlock = 0.01
dyBlock = 0.01

# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
dxBox = rDiMES*2.5
dyBox = rDiMES*2.5
dzBox = rDiMES*2.5

# Mesh parameters
LmeshDiMES = 0.003
LmeshBlock = 0.0005
LmeshBlockRefine = 0.0001  # increase mesh resolution at leadign edges
LmeshBound = 0.005
ddxBox = 0.0001
ddyBox = 0.0001
dyBoxRefine = 0.002
dxBoxRefine = 0.0005

# Tags
TagDiMES0 = 1
TagBox0 = 10
TagDiMES = 2000
TagBlock = [100*(i+1) for i in range(NBlocks)]

# Create DiMES surface
gmsh.model.occ.addDisk(xDiMES, yDiMES, zDiMES, rDiMES, rDiMES, TagDiMES0)

# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
TagBox0 = gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0)
gmsh.model.occ.remove([(3, TagBox0)])
gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0)], -1, removeTool=False)
s1 = gmsh.model.occ.getEntities(2)
TagBoundBox = list(set(s1) - set(s0))

# Create micro-trenches
for i in range(NBlocks):
    TagBlock[i] = gmsh.model.occ.addRectangle(
        xBlock[i], yBlock[i], zBlock[i], dxBlock, dyBlock)

TagDiMES = gmsh.model.occ.cut([(2, TagDiMES0)], [(
    2, TagBlock[i]) for i in range(NBlocks)], TagDiMES, removeTool=False)
TagExtBlock = TagBlock.copy()


for i in range(NBlocks):
    if dzBlock[i] > 0:
        TagExtBlock[i] = gmsh.model.occ.extrude(
            [(2, TagBlock[i])], 0, 0, dzBlock[i])

# Fusion two extruded blocks to remove duplicate neighbor surfaces
Cube1, a = gmsh.model.occ.fuse([(3, 2)], [(3, 1)])
gmsh.model.occ.remove(Cube1)

# Remove bottom surface from fused extruded cubes
gmsh.model.occ.remove([(2, 2003)])

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

### Setup mesh size
gmsh.model.mesh.field.add("Box", 1)
gmsh.model.mesh.field.setNumber(1, "VIn", LmeshBlock)
gmsh.model.mesh.field.setNumber(1, "VOut", 1)
gmsh.model.mesh.field.setNumber(1, "XMin", min(xBlock)-ddxBox)
gmsh.model.mesh.field.setNumber(1, "XMax", max(xBlock)+dxBlock+ddxBox)
gmsh.model.mesh.field.setNumber(1, "YMin", min(yBlock)-ddyBox)
gmsh.model.mesh.field.setNumber(1, "YMax", max(yBlock)+dyBlock+ddyBox)
gmsh.model.mesh.field.setNumber(1, "ZMin", 0)
gmsh.model.mesh.field.setNumber(1, "ZMax", max(dzBlock)*1.1)
gmsh.model.mesh.field.setNumber(1, "Thickness", dxBox)

gmsh.model.mesh.field.add("Box", 2)
gmsh.model.mesh.field.setNumber(2, "VIn", LmeshBlockRefine)
gmsh.model.mesh.field.setNumber(2, "VOut", 1)
gmsh.model.mesh.field.setNumber(2, "YMin", 0-dyBoxRefine)
gmsh.model.mesh.field.setNumber(2, "YMax", 0+dyBoxRefine)
gmsh.model.mesh.field.setNumber(2, "XMin", min(xBlock)-dxBoxRefine)
gmsh.model.mesh.field.setNumber(2, "XMax", max(xBlock)+dxBlock+dxBoxRefine)
gmsh.model.mesh.field.setNumber(2, "ZMin", 0)
gmsh.model.mesh.field.setNumber(2, "ZMax", max(dzBlock)*1.1)
gmsh.model.mesh.field.setNumber(2, "Thickness", dyBoxRefine)

gmsh.model.mesh.field.add("Ball", 3)
gmsh.model.mesh.field.setNumber(3, "VIn", LmeshDiMES)
gmsh.model.mesh.field.setNumber(3, "VOut", LmeshBound)
gmsh.model.mesh.field.setNumber(3, "XCenter", xDiMES)
gmsh.model.mesh.field.setNumber(3, "YCenter", yDiMES)
gmsh.model.mesh.field.setNumber(3, "ZCenter", zDiMES)
gmsh.model.mesh.field.setNumber(3, "Radius", rDiMES)

# Set mesh size fieldslist
gmsh.model.mesh.field.add("Min", 4)
gmsh.model.mesh.field.setNumbers(4, "FieldsList", [1, 2, 3])

# Turn off background(global) mesh size rules and set fieldslist as mesh size rules.
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.model.mesh.field.setAsBackgroundMesh(4)

# We can constraint the min and max element sizes to stay within reasonnable values
#gmsh.option.setNumber("Mesh.MeshSizeMin", 0.00001)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)

# Define groups to allow setting properties of elements when generating geometry input with GeomSetup
SetGroups(gmsh.model, TagDiMES, "DiMES", [255, 0, 0])
SetGroups(gmsh.model, TagExtBlock, "Blocks", [0, 0, 255])
SetGroups(gmsh.model, TagBoundBox, "BoundBox", [0, 0, 255])

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("leading_edges_DiMES.msh")

# close gmsh session. Comment when working on the geometry.
gmsh.finalize()
