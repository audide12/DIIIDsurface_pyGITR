# -*- coding: utf-8 -*-


from pyGITR.Geom import GeomSetup

# Load geometry
g = GeomSetup('simple_cube.msh', Verbose=True)

# Show existing groups
g.ShowGroups()

# Plot mesh for some groups
g.Plot()

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
g.SetAxisLim(-0.02, 0.02)

# Show normals for the boundbox
g.ShowNormals()

# Set properties for each group
g.SetElemAttr([],'surface', 1)


# set inDir. Empty list = all elements
g.SetElemAttr([], 'inDir',1)

# Set Z for material
g.SetElemAttr([], 'Z', 0)

# Set potential for biasing
g.SetElemAttr([],'potential',0)

# Plot geometry showing values of Z with color
g.Plot(ElemAttr='Z', Alpha=0.1)

# Write the geometry file
g.WriteGeomFile(Folder='input')



