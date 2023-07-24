import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from pyGITR.Particles import *
from pyGITR.Geom import *

import netCDF4
from netCDF4 import Dataset

def checkCoplanar(x,y,z,a,b,c,d):
    coplanar = False
    plane = a*x + b*y + c*z + d
    if round(plane,3) == 0:
        coplanar = True
    else: print(plane)
    return coplanar

# Load geomGITR
def loadCFG(geomFile = "gitrGeom.cfg"):
    with io.open(geomFile) as f:
        config = libconf.load(f)
        
    x1 = np.array(config['geom']['x1'])
    x2 = np.array(config['geom']['x2'])
    x3 = np.array(config['geom']['x3'])
    y1 = np.array(config['geom']['y1'])
    y2 = np.array(config['geom']['y2'])
    y3 = np.array(config['geom']['y3'])
    z1 = np.array(config['geom']['z1'])
    z2 = np.array(config['geom']['z2'])
    z3 = np.array(config['geom']['z3'])

    a = np.array(config['geom']['a'])
    b = np.array(config['geom']['b'])
    c = np.array(config['geom']['c'])
    d = np.array(config['geom']['d'])

    area = np.array(config['geom']['area'])
    plane_norm = np.array(config['geom']['plane_norm'])
    surface = np.array(config['geom']['surface'])
    inDir = np.array(config['geom']['inDir'])

    return x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,inDir

def offsetPoints(x,y,z,a,b,c,inDir):
    position = np.vstack((x,y,z)).T
    normal = np.array([a,b,c])
    norm = Norm(np.reshape(normal,(-1,3)))*inDir
    a_norm = a*0.0000001/norm
    b_norm = b*0.0000001/norm
    c_norm = c*0.0000001/norm
    norm = np.vstack((a_norm, b_norm, c_norm)).T
    offset_position = position+norm
    return offset_position[:,0], offset_position[:,1], offset_position[:,2]



def genPoints(numPoints,x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d):
    xr,yr,zr=[],[],[]
    genPoints = 0
    while genPoints < numPoints:
        x = (x1,x2,x3)
        y = (y1,y2,y3)
        z = (z1,z2,z3)

        xmin,xmax = min(x),max(x)
        ymin,ymax = min(y),max(y)
        zmin,zmax = min(z),max(z)

        switch = random.randint(0,3)
        if switch == 0:
            xrand = random.uniform(xmin,xmax)
            yrand = random.uniform(ymin,ymax)
            zrand = -(d + a*xrand + b*yrand)/c
            if checkCoplanar(xrand,yrand,zrand,a,b,c,d) and zrand >= zmin and zrand <= zmax:
                xr.append(xrand)
                yr.append(yrand)
                zr.append(zrand)
                genPoints+=1
        if switch == 1:
            xrand = random.uniform(xmin,xmax)
            zrand = random.uniform(zmin,zmax)
            yrand = -(d + a*xrand + c*zrand)/b
            if checkCoplanar(xrand,yrand,zrand,a,b,c,d) and yrand >= ymin and yrand <= ymax:
                xr.append(xrand)
                yr.append(yrand)
                zr.append(zrand)
                genPoints+=1
        if switch == 2:
            zrand = random.uniform(zmin,zmax)
            yrand = random.uniform(ymin,ymax)
            xrand = -(d + c*zrand + b*yrand)/a
            if checkCoplanar(xrand,yrand,zrand,a,b,c,d) and xrand >= xmin and xrand <= xmax:
                xr.append(xrand)
                yr.append(yrand)
                zr.append(zrand)
                genPoints+=1

    return xr,yr,zr

def fixNorm(a,b,c,inDir):
    normal = np.array([a,b,c])
    norm = Norm(np.reshape(normal,(-1,3)))*inDir
    a_norm = a*0.01/norm
    b_norm = b*0.01/norm
    c_norm = c*0.01/norm
    norm = np.vstack((a_norm, b_norm, c_norm)).T
    return norm

def plotPointsAndElement(surface_x,surface_y,surface_z,xr,yr,zr,a,b,c,inDir):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot3D(surface_x,surface_y,surface_z)
    ax.scatter3D(xr,yr,zr)
    norm = fixNorm(a,b,c,inDir)
    ax.quiver(np.average(surface_x[0:3]),np.average(surface_y[0:3]),np.average(surface_z[0:3]),norm[:,0],norm[:,1],norm[:,2],color='m')
    ax.set_xlim([min(surface_x),min(surface_x)+0.01])
    ax.set_ylim([min(surface_y),min(surface_y)+0.01])
    ax.set_zlim([min(surface_z),min(surface_z)+0.01])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    plt.show()

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    # return the value in radians
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(Dot(v1_u, v2_u), -1.0, 1.0))\

# Distribution can be also defined with a user-defined pdf. Must be formatted as
# f(x, **kwargs)
def LevyDistrib(x, c=1, mu=0):
    return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)

def loadNC(File):
    ParticleData = Dataset(File, "r", format="NETCDF4")
    # nP = len(ParticleData['nP'])
    nP = len(ParticleData.dimensions['nP']) # number of timesteps
    x = np.array(ParticleData['x'])
    y = np.array(ParticleData['y'])
    z = np.array(ParticleData['z'])
    vx = np.array(ParticleData['vx'])
    vy = np.array(ParticleData['vy'])
    vz = np.array(ParticleData['vz'])
    v = np.vstack((vx,vy,vz)).T

    return nP,x,y,z,vx,vy,vz


def checkRotation(surface_x,surface_y,surface_z,xr,yr,zr,a,b,c,inDir):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot3D(surface_x,surface_y,surface_z)
    ax.scatter3D(xr,yr,zr)
    ax.set_xlim([min(surface_x),min(surface_x)+0.01])
    ax.set_ylim([min(surface_y),min(surface_y)+0.01])
    ax.set_zlim([min(surface_z),min(surface_z)+0.01])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    plt.show()




# thetaB = -90 #degree by default, set Degree=False in RotateVector for radian
# AxisrotB = [1,0,0]  # [0,1,0]

# # Check results of the rotation
# B0 = [0,0,1]
# B = RotateVector(B0, AxisrotB, thetaB)
# # PlotVector(B)


# initial = [[0,0,1],[1,1,1]]
# final = [[1,1,1],[0,0,0]]
# angles = angle_between(initial,final)
# print(angles)


# initial_ = np.reshape(initial,(-1,3))
# final_ = np.reshape(final,(-1,3))
# cross = np.reshape(Cross(initial_, final_),(-1,))



# new = RotateVector(initial, cross, angles, Degree=False)
# PlotVector(new)
# plt.show()


# PlotVector(final)
# plt.show()












# def sample_areas(triangles, samples):
#   # compute and sum triangle areas
#   totalA = 0.0
#   areas = []
#   for t in triangles:
#     a = t.area()
#     areas.append(a)
#     totalA += a

#   # compute and sort random numbers from [0,1)
#   rands = sorted([random.random() for x in range(samples)])

#   # sample based on area
#   area_limit = 0.0
#   rand_index = 0
#   rand_value = rands[rand_index]
#   for i in range(len(areas)):
#     area_limit += areas[i]
#     while rand_value * totalA < area_limit:
#       # sample randomly over current triangle
#      triangles[i].add_random_sample()

#       # advance to next sorted random number
#       rand_index += 1;
#       if rand_index >= samples:
#         return
#       rand_value = rands[rand_index]