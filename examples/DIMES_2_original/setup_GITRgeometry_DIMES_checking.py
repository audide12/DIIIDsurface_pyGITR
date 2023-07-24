#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 18:12:10 2023

@author: de
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:56:05 2022

@author: de
"""

# -*- coding: utf-8 -*-

import os
 


# Load geometry
g = GeomSetup('DiMES_SiC_checking.msh', Verbose=True)

# Show existing groups
g.ShowGroups()

# Plot mesh for some groups
g.Plot_Geom(["DiMES"])
#g.Plot(["BoundBox"])
g.Plot_Geom(["LargeRectangle"])
g.Plot_Geom(["LeftStrip"])
g.Plot_Geom(["RightStrip"])

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
# g.SetAxisLim(-1.25, 1.25)

# Show Centroids
# g.ShowCentroids()

# Show normals for the DiMES
# g.ShowNormals("DiMES")
# g.ShowNormals("LargeRectangle")
# g.ShowNormals("LeftStrip")
# g.ShowNormals("RightStrip")


# Set properties for each group
# set 'surface' property. Default is 0.
g.SetElemAttr(["DiMES"], 'surface', 1)
g.SetElemAttr(["LargeRectangle"], 'surface', 1)
g.SetElemAttr(["LeftStrip"], 'surface', 1)
g.SetElemAttr(["RightStrip"], 'surface', 1)
g.SetElemAttr(["BoundBox"], 'surface', 0)

# set inDir. Empty list = all elements
g.SetElemAttr([],'inDir',1)
g.SetElemAttr(["BoundBox"], 'inDir',-1) # changed this from +1
#g.ShowInDir(["BoundBox"])
#g.ShowInDir(["DiMES","LargeRectangle","LeftStrip","RightStrip"])


# Set Z for material
g.SetElemAttr(["DiMES"], 'Z', 6) # for Carbon head
g.SetElemAttr(["LargeRectangle"], 'Z', 20) # for Silicon Carbide coating
g.SetElemAttr(["LeftStrip"], 'Z', 20) # for Silicon Carbide coating
g.SetElemAttr(["RightStrip"], 'Z', 20) # for Silicon Carbide coating

#g.SetElemAttr(["Large Dot","Small Dot","Metal Ring"], 'Z', 74)
g.SetElemAttr(["BoundBox"], 'Z', 0)

# Set potential for biasing. Empty list = all elements
g.SetElemAttr([],'potential',0)

# Plot geometry showing values of Z with color
#g.Plot(ElemAttr='Z', Alpha=0.1)

# Write the geometry file
g.WriteGeomFile(Folder="input",OverWrite=True)

# os.system("mv gitrGeom.cfg input")



