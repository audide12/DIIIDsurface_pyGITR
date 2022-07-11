#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for micro-trenches experiments on DiMES
# Reference paper:
@author: guterl
"""
import gmsh
gmsh.initialize()

import math

# Clean up gmsh session if script already executed
try: gmsh.finalize()
except:
    pass
# initialize gsmh session
gmsh.initialize()


#generate Dimes surface
xDiMES=0.0
yDiMES=0.0
zDiMES=0.0
rDiMES=0.05
TagDiMES0=1000
TagDiMES=1111
TagTrench=2222
MeshSizeDiMES=0.01
#2D surface for DiMES
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0)


# Add elements
axTrench,ayTrench,azTrench=0,0,1
xTrench,yTrench,zTrench=0,0,0
alphaTrench=math.pi/4
lTrench,LTrench=0.001,0.001
dxTrench,dyTrench,dzTrench=0.01,0.01,0.01
TagObject=2000
g=gmsh.model.occ
#gmsh.model.occ.addRectangle(xTrench,yTrench,zTrench,lTrench,LTrench,TagObject)
gmsh.model.occ.addBox(xTrench, yTrench, zTrench, dxTrench, dyTrench, dzTrench,TagTrench)
gmsh.model.occ.cut([(2,TagDiMES0)],[(3,TagTrench)],TagDiMES,removeTool=False)
print(gmsh.model.occ.getEntities())
#gmsh.model.occ.addBox(xTrench, yTrench, zTrench, dxTrench, dyTrench, dzTrench,TagTrench)
#gmsh.model.occ.remove([(2,1111)])
#gmsh.model.occ.remove([(2,1112)])
#gmsh.model.occ.remove([(2,1113)])
gmsh.model.occ.remove([(3,TagTrench)])
print(gmsh.model.occ.getEntities())
#gmsh.model.occ.remove([(2,1113)])
#gmsh.model.occ.remove([(2,1114)])
#gmsh.model.occ.remove([(2,1115)])
gmsh.model.occ.remove([(2,1005)])
gmsh.model.occ.fuse([(2,TagDiMES)], [(2,1003), (2,1004)])
print(gmsh.model.occ.getEntities())

gmsh.model.occ.mesh.setSize([(0, i) for i in [2,3,4,5,6,7,8,9]]
                             , 0.001)
gmsh.model.occ.synchronize()
#gmsh.model.occ.extrude(dimTags, dx, dy, dz)
#gmsh.model.occ.rotate([(2,TagObject)],xTrench,yTrench,zTrench,axTrench,ayTrench,azTrench,alphaTrench)
#gmsh.model.occ.cut([(2,TagDiMES0)],[(2,TagObject)],TagDiMES,removeTool=False)
# Synchronize necessary before mesh setup and generation

#%%
# Setup mesh
gmsh.option.setNumber("Geometry.NumSubEdges", 100)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
# And we set the minimum number of elements per 2*Pi radians:
#gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)

# We can constraint the min and max element sizes to stay within reasonnable
# values (see `t10.py' for more details):
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.00001)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)
# Generate mesh
gmsh.model.mesh.generate(2)

# Save and Display
gmsh.write("t1000.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#gmsh.finalize()
# microtrench generation
#gmsh.model.geo.addPoint(0.0, 0.0, 0, lc, 1)