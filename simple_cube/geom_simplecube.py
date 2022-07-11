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

# Define center cube and length
xBox = 0.0
yBox = 0.0
zBox = 0.0
L = 0.001



# Mesh parameters
Lmesh = 0.00

# Tags
TagDiMES0 = 1
TagBox0 = 10
TagDiMES = 2000


# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
TagBox0 = gmsh.model.occ.addBox(xBox-L, yBox-L, zBox-L, L, L, L, TagBox0)
gmsh.model.occ.remove([(3, TagBox0)])

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# We can constraint the min and max element sizes to stay within reasonnable values
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("simple_cube.msh")

# close gmsh session. Comment when working on the geometry.
gmsh.finalize()
