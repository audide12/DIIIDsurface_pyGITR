"""
Generation of particles distribution for GITR.
@author: zack bergstrom

To Do:
points are close to surface, but not always within the surface
scale velocities: to energy, vpara, vperp
rotate to magnetic field reference frame?
need to more velocity rotations into the surface loop
"""
from pyGITR.particleSource_functions import *
import os

def makeParticleSource(data,geomFile,particleFile,particleWeightFile):
    x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,indir = loadCFG(geomFile=geomFile)
    
    
    particleFile1 = '/Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_temp.nc'

    # Generate positions per mesh element
    nP = 0
    x,y,z = [],[],[]
    _a,_b,_c = [],[],[]
    _indir = []
    skip=1
    for num,i in enumerate(data.keys()):
        if num%skip==0:
            # print(i,data[i])
            nP+=data[i]
            xr,yr,zr = genPoints(data[i],x1[i],x2[i],x3[i],y1[i],y2[i],y3[i], \
                                    z1[i],z2[i],z3[i],a[i],b[i],c[i],d[i])
            xr,yr,zr = offsetPoints(xr,yr,zr,a[i],b[i],c[i],indir[i])
            x.extend(xr)
            y.extend(yr)
            z.extend(zr)

            surface_x, surface_y, surface_z = [x1[i],x2[i],x3[i],x1[i]], \
                    [y1[i],y2[i],y3[i],y1[i]], [z1[i],z2[i],z3[i],z1[i]]
            #zmin,zmax = min(zr), max(zr)

            # plotPointsAndElement(surface_x,surface_y,surface_z,xr,yr,zr,a[i],b[i],c[i],indir[i])

            for j in range(0,data[i]):
                _a.append(a[i])
                _b.append(b[i])
                _c.append(c[i])
                _indir.append(indir[i])


    x,y,z = np.array(x),np.array(y),np.array(z)
    print("total number of particles in source:",nP)



    # Populate NC file with positions and velocities
    p = ParticleDistribution()
    p.SetAttr('Np', nP)

    # Set positions of particles
    p.SetAttr('x',x)
    p.SetAttr('y',y)
    p.SetAttr('z',z)

    # Set velocities of particles
    #p.SetAttr(['vx','vy'],'Gaussian')
    #p.SetAttr(['vz'],LevyDistrib, x=np.linspace(0.001,10,1000), c=2, mu=0)
    
    # p.SetAttr(['vz'],'Gaussian',sigma = 1.825e4,beta=3.16e14)
    # p.SetAttr(['vy','vx'],'Gaussian',sigma = 1.825e4,beta=3.16e14)

    # vpara = 1.0
    # vperp = 10.0 #5
    # #p.ScaleAttr(['vx','vz'],vperp)
    
    # p.ScaleAttr('vx',vperp)
    
    # p.ScaleAttr('vy',vpara)
    
    
    p.SetAttr(['vx','vy'],'Gaussian_Jerome')
    
    #p.SetAttr(['vz'],'Levy', x=np.linspace(0.001,10,1000), c=2, mu=0)   # without importance sampling

    p.SetAttr(['vz'],'Gaussian') # with the correct importance interval
    p.SetAttr(['vz'],'Levy', x=np.linspace(0.001,10,1000), c=2, mu=0)
    
    vpara = 1e4
    vperp = 1e3
    p.ScaleAttr(['vx','vz'],vperp)
    p.ScaleAttr('vy',vpara)
    
    w = np.zeros(nP)
    for i in range(len(w)):
        w[i] = p['weight'][i]/p['vz'][i]
    p.SetAttr(['weight'],w)


    # Write particle distribution in netcdf file
    p.WriteParticleFile(particleFile1)



    # Rotate velocities parallel to mesh element normals
    nP,x,y,z,vx,vy,vz = loadNC(particleFile1)
    v = np.vstack((vx,vy,vz)).T

    norm = np.vstack((_a,_b,_c)).T  # length number of particles

    init = np.copy(norm) # length number of particles
    init[:,0], init[:,1], init[:,2] = 0,0,1

    angles = angle_between(init,norm) # in radians, wrt surface norm
    print(angles)
    AxisRot = Cross(v,norm)

    vxr,vyr,vzr = [],[],[]
    for i,angle in enumerate(angles):
        vx_rot, vy_rot, vz_rot = RotateVector(v[i], AxisRot[i], angle, Degree=False)
        vxr.append(vx_rot)
        vyr.append(vy_rot)
        vzr.append(vz_rot)

    v_rotated = np.vstack((vxr,vyr,vzr)).T

    # fix v_rotated for the indir
    for i in range(len(v_rotated)):
        v_rotated[i] = v_rotated[i]*_indir[i]

    v_rot = fixNorm(v_rotated[:,0],v_rotated[:,1],v_rotated[:,2],1) # fixes up the size of the vector for plotting
    _v = fixNorm(v[:,0],v[:,1],v[:,2],1) # fixes up the size of the vector for plotting

    # Plot results, mesh element, positions, and velocities
    particle_idx = 0
    for num,i in enumerate(data.keys()):
        if num%skip==0:
            # print(i,data[i])

            _norm = fixNorm(a[i],b[i],c[i],indir[i])

            #fig = plt.figure()
            #ax = plt.axes(projection="3d")

            surface_x, surface_y, surface_z = [x1[i],x2[i],x3[i],x1[i]], \
                    [y1[i],y2[i],y3[i],y1[i]], [z1[i],z2[i],z3[i],z1[i]]
            #ax.plot3D(surface_x,surface_y,surface_z)

            # ax.scatter3D(x,y,z) # plotting all point, not just those on the mesh  -- from before
            
            # ax.scatter3D(x[particle_idx:particle_idx+data[i]],\
            #     y[particle_idx:particle_idx+data[i]],\
            #         z[particle_idx:particle_idx+data[i]]) # plotting all point, not just those on the mesh


            # ax.quiver(x,y,z,v_rot[:,0],v_rot[:,1],v_rot[:,2],color='g') # rotated
            # ax.quiver(x,y,z,_v[:,0],_v[:,1],_v[:,2],color='r') # not rotated
            # ax.quiver(np.average(surface_x[0:3]),np.average(surface_y[0:3]),np.average(surface_z[0:3]),_norm[:,0],_norm[:,1],_norm[:,2],color='m') # rotated

            # ax.quiver(x[particle_idx:particle_idx+data[i]],\
            #     y[particle_idx:particle_idx+data[i]],\
            #         z[particle_idx:particle_idx+data[i]],\
            #             v_rot[particle_idx:particle_idx+data[i],0],\
            #                 v_rot[particle_idx:particle_idx+data[i],1],\
            #                     v_rot[particle_idx:particle_idx+data[i],2],color='g') # rotated
            # ax.quiver(x[particle_idx:particle_idx+data[i]],\
            #     y[particle_idx:particle_idx+data[i]],\
            #         z[particle_idx:particle_idx+data[i]],\
            #             _v[particle_idx:particle_idx+data[i],0],\
            #                 _v[particle_idx:particle_idx+data[i],1],\
            #                     _v[particle_idx:particle_idx+data[i],2],color='r') # not rotated
            # ax.quiver(np.average(surface_x[0:3]),np.average(surface_y[0:3]),np.average(surface_z[0:3]),\
            #     _norm[:,0],_norm[:,1],_norm[:,2],color='m') # rotated


            particle_idx += data[i]

            # ax.set_xlim([min(surface_x)-0.005,max(surface_x)+0.01])
            # ax.set_ylim([min(surface_y)-0.005,max(surface_y)+0.01])
            # ax.set_zlim([min(surface_z)-0.005,max(surface_z)+0.01])

            # ax.set_xlabel('x [m]')
            # ax.set_ylabel('y [m]')
            # ax.set_zlabel('z [m]')
            # plt.show()

    # Save data in netcdf format
    
    ncFile = netCDF4.Dataset(particleWeightFile, "w", format="NETCDF4")
    particle_dim = ncFile.createDimension('particle_dim', nP) # number of particles
    p_number = ncFile.createVariable('particle_number',np.float64,('particle_dim'))
    p_weight = ncFile.createVariable('particle_wt',np.float64,('particle_dim'))
    p_number[:] = p_number
    p_weight[:] = p_weight
    ncFile.close()
        
    
    os.system("rm /Users/de/Research/DIIIDsurface_pyGITR/examples/DIMES/input/particleConf_temp.nc")
    
    rootgrp = netCDF4.Dataset(particleFile, "w", format="NETCDF4")
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

#if __name__ == "__main__":
    # x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,indir = loadCFG()
    # data = {}
    # for i in range(len(surface)):
    #     if surface[i]==1:
    #         data[i] = 1
    # print(len(data))

    #data = {3:3,30:3,300:3,3000:3,30000:3} # mesh element index : number of particles
    
    #makeParticleSource(data, "../input/gitrGeom.cfg", "particleSourceForGITRm.nc")

