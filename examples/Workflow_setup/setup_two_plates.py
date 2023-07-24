#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:54:17 2022

@author: audide
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
height = 3.0
width = 3.0
z_down = 0.0
z_up = 0.5
x_center = 0.0
y_center = 0.0
tag_down = 1
tag_up = 2



# Create DiMES surface
gmsh.model.occ.addRectangle(x_center, y_center, z_down, height, width, tag_down)
gmsh.model.occ.addRectangle(x_center, y_center, z_up, height, width, tag_up)





gmsh.model.occ.synchronize()



# We can constraint the min and max element sizes to stay within reasonnable values
gmsh.option.setNumber("Mesh.MeshSizeMin", 10)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 12)



# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("two_plates.msh")

# close gmsh session. Comment when working on the geometry.
gmsh.finalize()
