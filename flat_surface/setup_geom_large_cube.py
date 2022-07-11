# -*- coding: utf-8 -*-


from pyGITR.Geom import GeomSetup
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
L = 0.05



# Mesh parameters
Lmesh = 2*L

# Tags
TagDiMES0 = 1
TagBox0 = 10
TagBox1 = 11
TagDiMES = 2000


# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
TagBox0 = gmsh.model.occ.addRectangle(xBox-L, yBox-L, zBox, 2*L, 2*L,TagBox0)
TagBox1 = gmsh.model.occ.addRectangle(xBox-L, yBox-L, zBox+L, 2*L, 2*L, TagBox1)
#gmsh.model.occ.remove([(3, TagBox0)])

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# We can constraint the min and max element sizes to stay within reasonnable values
gmsh.option.setNumber("Mesh.MeshSizeMin", Lmesh)
gmsh.option.setNumber("Mesh.MeshSizeMax", Lmesh)

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
#gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("large_box.msh")

# close gmsh session. Comment when working on the geometry.
gmsh.finalize()
# Load geometry
g = GeomSetup('large_box.msh', Verbose=True)

# Show existing groups
g.ShowGroups()

# Plot mesh for some groups
g.Plot()

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
g.SetAxisLim(-0.10, 0.10)

# Show normals for the boundbox
g.ShowNormals(L=0.02)

# Set properties for each group
g.SetElemAttr([],'surface', 1)


# set inDir. Empty list = all elements
g.SetElemAttr([], 'inDir',-1)

# Set Z for material
g.SetElemAttr([], 'Z', 0)

# Set potential for biasing
g.SetElemAttr([],'potential',0)
g.SetAttr('lambda',0.001)
g.SetAttr('zref',0)
g.SetAttr('periodic_bc_x',1)
g.SetAttr('periodic_bc_y',1)
g.SetAttr('periodic_bc_x0',-L)
g.SetAttr('periodic_bc_x1',+L)
g.SetAttr('periodic_bc_y0',-L)
g.SetAttr('periodic_bc_y1',+L)
# Plot geometry showing values of Z with color
g.Plot(ElemAttr='Z', Alpha=0.1)


# Write the geometry file
g.WriteGeomFile(Folder='input',OverWrite=True)



