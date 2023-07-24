"""
Generation of particles distribution for GITR.
@author: zack bergstrom

To Do:
displace point inward
rotate parallel to normal vector
points are close to surface, but not always within the surface
scale velocities: to energy, vpara, vperp
rotate to parallel with surface norm
rotate to magnetic field reference frame?
need to more velocity rotations into the surface loop
"""
from particleSource_functions import *

#x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,indir = loadCFG()
x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,indir = loadCFG("gitrGeom.cfg")
data = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1}
# data = {0:300,10:300,100:300,1000:300} # mesh element index : number of particles
#data = {1:3,10:3,100:3,1000:3} # mesh element index : number of particles

# data = {}
# for i in range(len(surface)):
#     if surface[i]==1:
#         data[i] = 1
# print(len(data))

# Generate positions per mesh element
nP = 0
x,y,z = [],[],[]
_a,_b,_c = [],[],[]
for i in data.keys():
    print(i,data[i])
    nP+=data[i]
    xr,yr,zr = genPoints(data[i],x1[i],x2[i],x3[i],y1[i],y2[i],y3[i], \
                            z1[i],z2[i],z3[i],a[i],b[i],c[i],d[i])

    print(xr,yr,zr)
    xr,yr,zr = offsetPoints(xr,yr,zr,a[i],b[i],c[i],indir[i])
    print(xr,yr,zr)
    x.extend(xr)
    y.extend(yr)
    z.extend(zr)

    surface_x, surface_y, surface_z = [x1[i],x2[i],x3[i],x1[i]], \
            [y1[i],y2[i],y3[i],y1[i]], [z1[i],z2[i],z3[i],z1[i]]
    zmin,zmax = min(zr), max(zr)

    plotPointsAndElement(surface_x,surface_y,surface_z,xr,yr,zr,a[i],b[i],c[i],indir[i])

    for j in range(0,data[i]):
        _a.append(a[i])
        _b.append(b[i])
        _c.append(c[i])

x,y,z = np.array(x),np.array(y),np.array(z)
print("total number of particles in source:",nP)



# Populate NC file with positions and velocities
ParticleFile='particleSource.nc'
p = ParticleDistribution()
p.SetAttr('Np', nP)
# p.ShowAvailablePdfs()

# Set positions of particles
p.SetAttr('x',x)
p.SetAttr('y',y)
p.SetAttr('z',z)

# Set velocities of particles
p.SetAttr(['vx','vy'],'Gaussian')
p.SetAttr(['vz'],LevyDistrib, x=np.linspace(0.001,10,1000), c=2, mu=0)

vpara = 1e4
vperp = 1e5
p.ScaleAttr(['vx','vz'],vperp)
p.ScaleAttr('vy',vpara)

# Write particle distribution in netcdf file
p.WriteParticleFile(ParticleFile)



# Rotate velocities parallel to mesh element normals
nP,x,y,z,vx,vy,vz = loadNC(ParticleFile)
v = np.vstack((vx,vy,vz)).T

norm = np.vstack((_a,_b,_c)).T  # length number of mesh elements
_norm = np.reshape(norm,(-1,3))  # length number of mesh elements

init = np.copy(norm) # length number of particles
init[:,0] = 0
init[:,1] = 0
init[:,2] = 1
# _init = np.reshape(init,(-1,3))

angles = angle_between(init,norm) #in radians, wrt surface norm
AxisRot = np.reshape(Cross(v, _norm),(-1,3))

vxr,vyr,vzr = [],[],[]
for i,angle in enumerate(angles):
    vx_rot, vy_rot, vz_rot = RotateVector(v[i], AxisRot[i], angle, Degree=False)
    vxr.append(vx_rot)
    vyr.append(vy_rot)
    vzr.append(vz_rot)

v_rotated = np.vstack((vxr,vyr,vzr)).T

# v = unit_vector(v)
# v_rotated = unit_vector(v_rotated)

_v = Norm(np.reshape(v,(-1,3)))
av = v[:,0]*0.01/_v
bv = v[:,1]*0.01/_v
cv = v[:,2]*0.01/_v
_v = np.vstack((av,bv,cv)).T

# Plot results, mesh element, positions, and velocities
for i in data.keys():
    # print(i,data[i])
    v_rot = fixNorm(v_rotated[:,0],v_rotated[:,1],v_rotated[:,2],indir[i]) # fixes up the size of the vector for plotting
    # _v = fixNorm(v[:,0],v[:,1],v[:,2],indir[i]) # fixes up the size of the vector for plotting
    _norm = fixNorm(a[i],b[i],c[i],indir[i])
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    surface_x, surface_y, surface_z = [x1[i],x2[i],x3[i],x1[i]], \
            [y1[i],y2[i],y3[i],y1[i]], [z1[i],z2[i],z3[i],z1[i]]
    ax.plot3D(surface_x,surface_y,surface_z)

    ax.scatter3D(x,y,z) # plotting all point, not just those on the mesh

    ax.quiver(x,y,z,v_rot[:,0],v_rot[:,1],v_rot[:,2],color='g') # rotated
    ax.quiver(x,y,z,_v[:,0],_v[:,1],_v[:,2],color='r') # not rotated
    ax.quiver(np.average(surface_x[0:3]),np.average(surface_y[0:3]),np.average(surface_z[0:3]),_norm[:,0],_norm[:,1],_norm[:,2],color='m') # rotated


    ax.set_xlim([min(surface_x)-0.005,min(surface_x)+0.01])
    ax.set_ylim([min(surface_y)-0.005,min(surface_y)+0.01])
    ax.set_zlim([min(surface_z)-0.005,min(surface_z)+0.01])

    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    plt.show()



# Save data in netcdf format
rootgrp = netCDF4.Dataset("particleSource_test.nc", "w", format="NETCDF4")
npp = rootgrp.createDimension("nP", nP)
xxx = rootgrp.createVariable("x","f8",("nP"))
yyy = rootgrp.createVariable("y","f8",("nP"))
zzz = rootgrp.createVariable("z","f8",("nP"))
vxx = rootgrp.createVariable("vx","f8",("nP"))
vyy = rootgrp.createVariable("vy","f8",("nP"))
vzz = rootgrp.createVariable("vz","f8",("nP"))
xxx[:] = x
yyy[:] = y
zzz[:] = z
vxx[:] = v_rotated[:,0]
vyy[:] = v_rotated[:,1]
vzz[:] = v_rotated[:,2]
rootgrp.close()





