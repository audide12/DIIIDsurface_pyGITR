#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:14:06 2022

@author: de
"""

#!/usr/bin/env python3

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gmsh
from pyGITR.gmsh_helper import SetGroups

# Initialize gmsh session
gmsh.initialize()

# shift entire mesh to align with center of DiMES
r_shift = 1.485
z_shift = -1.250

# Define DiMES cap surface
xDiMES=0.0+r_shift
yDiMES=0.0
zDiMES=0.0
# zDiMES=z_shift
rDiMES=0.025



# Define inside rectangle
xRectDimes = xDiMES
yRectDimes = 0.0
zRectDimes = 0.0
L1_rect = 0.02
L2_rect = 0.025
Strip_distance = 0.006   # verify this
Strip_width = 0.001



# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
# dxBox = rDiMES*8.5
# dyBox = rDiMES*8.5
# dzBox = rDiMES*8.5


dxBox = rDiMES*3
dyBox = rDiMES*3
dzBox = rDiMES*3

dxBox = rDiMES*10
dyBox = rDiMES*10
dzBox = rDiMES*10


# Define metal rings
dxBlock = 0.05
dyBlock = dyBox
xBlock=xDiMES-0.01-rDiMES-dxBlock
yBlock=yDiMES
zBlock=zDiMES

# Tags
TagDiMES0=0
TagDiMES=10
TagRectDimes=20
TagStrip_L_Dimes=50
TagStrip_R_Dimes=60
TagBox0=40



# When add a shape, return an int tag, either automatically, or the tag given
# Get entities returns the dimensionality and tag of the current entities
# Cut removes a list of specified entities (via their tags), and returns ???

# there are 3d entities, and 2d entities
# when the 3d box is made, 1 3d entitiy with the assigned tag is made
# along with the 2d entities that comprise the 3d shape


# Create 2D surfaces
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0) # the second r is used for elipse definition

gmsh.model.occ.addRectangle(xRectDimes - L1_rect/2, yRectDimes - L2_rect/2, zRectDimes, L1_rect, L2_rect, tag=TagRectDimes, roundedRadius=0.)

#Left strip
gmsh.model.occ.addRectangle(xRectDimes - L1_rect/2 - Strip_distance - Strip_width , yRectDimes - L2_rect/2, zRectDimes, Strip_width, L2_rect, tag=TagStrip_L_Dimes, roundedRadius=0.)

#Right strip
gmsh.model.occ.addRectangle(xRectDimes + L1_rect/2 + Strip_distance , yRectDimes - L2_rect/2, zRectDimes, Strip_width, L2_rect, tag=TagStrip_R_Dimes, roundedRadius=0.)


# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0)
gmsh.model.occ.remove([(3, TagBox0)]) # removes the 3d shape from the 3d shape list, leaving only 2d shapes in the 2d shape list

gmsh.model.occ.cut([(2, TagDiMES0)], [(2, TagRectDimes),(2,TagStrip_R_Dimes),(2,TagStrip_L_Dimes)], -1, removeTool=False)


s1 = gmsh.model.occ.getEntities(2)
TagBoundBox = list(set(s1) - set(s0))
print(type(s1))
print('s1',s1)
print(TagBoundBox)



DiMESTag = gmsh.model.occ.getEntities(2)

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Set number of elements on the boundary of each dots and DiMES cap
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 50)

#gmsh.option.setNumber("Mesh.MeshSizeMax", L1_rect) # to increase the mesh sizes


# Define groups to allow setting properties of elements when generating geometry input with GeomSetup

TagRect = [(2,20)]
TagStrip_L = [(2,50)]
TagStrip_R = [(2,60)]

SetGroups(gmsh.model, DiMESTag, "DiMES", [255, 0, 0])
SetGroups(gmsh.model, TagBoundBox, "BoundBox", [230, 230, 250])
SetGroups(gmsh.model, TagRect, "LargeRectangle", [0, 255, 0])
SetGroups(gmsh.model, TagStrip_L, "LeftStrip", [0, 0, 250])
SetGroups(gmsh.model, TagStrip_R, "RightStrip", [0, 255, 230])


# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("DiMES_SiC.msh")

# Close gmsh session
gmsh.finalize()