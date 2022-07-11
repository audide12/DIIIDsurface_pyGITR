# -*- coding: utf-8 -*-


from pyGITR.Geom import GeomSetup

# Load geometry
g = GeomSetup('/Users/guterl/Dropbox/python/pyGITR/examples/mesh/leading_edges_DiMES.msh', Verbose=True)

# Show existing groups
g.ShowGroups()

# Plot mesh for some groups
g.Plot(["DiMES","Blocks"])

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
g.SetAxisLim(-0.02, 0.02)

# Plot mesh for some groups
g.Plot(["BoundBox"])

# Show normals for the boundbox
g.ShowNormals(["BoundBox"])

# Set properties for each group
g.SetElemAttr(["BoundBox"], 'surface', 0)

# set 'surface' property. Default is 0.
g.SetElemAttr(["Blocks"], 'surface', 1)

# set inDir. Empty list = all elements
g.SetElemAttr([], 'inDir',1)

# Set Z for material
g.SetElemAttr(["DiMES"], 'Z', 6)
g.SetElemAttr(["Blocks"], 'Z', 74)

# Set potential for biasing
g.SetElemAttr([],'potential',0)

# Plot geometry showing values of Z with color
g.Plot(ElemAttr='Z', Alpha=0.1)

# Write the geometry file
g.WriteGeomFile()



