#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 20:10:56 2021

@author: guterl
"""

import numpy as np
import libconf
import io
import click
import os
import gmsh
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import collections as mc
from pyGITR.math_helper import *

class GeomGroup():
    '''Methods to handle groups of elements.'''

    def ShowGroups(self):
        print('-- Groups --')
        print('-- Dim :  ID :  Name')
        for p in self.Groups:
            print('-- {} : {} : "{}" '.format(2, p[0], p[1]))

    def GetGroups(self):
        self.PhysicalGroups = self.Model.getPhysicalGroups(2)
        self.Groups = [(g[1], self.Model.getPhysicalName(2, g[1]))
                       for g in self.PhysicalGroups]

    def GetGroupElements(self, GroupID: str or int or list):

        SelectedGroups = []
        if type(GroupID) != list:
            GroupID = [GroupID]
        for gID in GroupID:
            if type(gID) == str or type(gID) == int:
                Groups = [g[0]
                          for g in self.Groups if gID == g[0] or gID == g[1]]
                if len(Groups) == 0:
                    print('Cannot find {} in existing groups: {}'.format(
                        gID, self.Groups))
                else:
                    SelectedGroups.extend(Groups)
            else:
                print('GroupID must be a str or int or a list of str or int.')
                return np.array([])
        Idx = np.array([])

        if GroupID == []:
            ListEntities = self.Model.getEntities(2)
            if self.Verbose:
                print('Selected group:', [], ' Entities:', ListEntities)
            for L in [L[1] for L in ListEntities]:
                Elems = self.Mesh.getElements(2, L)[1][0]
                if self.Verbose:
                    print('Entities:', L, ' elems:', Elems)
                Idx = np.append(Idx, Elems)
        else:
            for G in SelectedGroups:
                ListEntities = self.Model.getEntitiesForPhysicalGroup(2, G)
                if self.Verbose:
                    print('Selected group:', G, ' Entities:', ListEntities)
                for L in ListEntities:
                    Elems = self.Mesh.getElements(2, L)[1][0]
                    if self.Verbose:
                        print('Entities:', L, ' elems:', Elems)
                    Idx = np.append(Idx, Elems)

        return Idx.astype('int')

    def GetGroupIdx(self, GroupID):
        if GroupID is not None and GroupID != []:
            Idx = self.GetGroupElements(GroupID)
            Idx = Idx-1
        else:
            Idx = np.arange(0, self.NElems, dtype=int)
        return Idx

class GeomInput(GeomGroup):

    def __init__(self):
        self.GeomInput = {}
        self.NElems = 0
        self.DefineDefaultAttributes()

    def DefineDefaultAttributes(self):
        self.DefaultsAttr = {}
        self.DefaultsAttr['3D'] = {}
        self.DefaultsAttr['3D']['Vector'] = {
            'x1': 0.0, 'x2': 0.0, 'x3': 0.0,
            'y1': 0.0, 'y2': 0.0, 'y3': 0.0,
            'z1': 0.0, 'z2': 0.0, 'z3': 0.0,
            'a': 0.0, 'b': 0.0, 'c': 0.0, 'd': 0.0,
            'inDir': 1, 'plane_norm': 0,
            'surface': 0,
            'potential': 0.0,
            'Z': 0,
            'area': 1.0
        }

        self.DefaultsAttr['3D']['Scalar'] = {'periodic': 0,
                                             'periodic_bc_x': 0,
                                             'periodic_bc_y': 0,
                                             'theta0': 0.0,
                                             'theta1': 0.0,
                                             'periodic_bc_x0': 0.0,
                                             'periodic_bc_x1': 0.0,
                                             'periodic_bc_y0': 0.0,
                                             'periodic_bc_y1': 0.0
                                             }

    def SetDefaultsAttr(self, GeoDim='3D'):
        assert GeoDim == '3D', 'only 3D implemented for the moment'
        assert self.NElems > 0, 'Number of elements must be >0 to set \
                                defaults input values'

        Defaults = self.DefaultsAttr['3D']
        for k, v in Defaults['Vector'].items():
            self.GeomInput[k] = np.full((self.NElems), v)
        for k, v in Defaults['Scalar'].items():
            self.GeomInput[k] = v

    def ShowDefaultsAttr(self, GeoDim='3D'):
        Defaults = self.DefaultsAttr[GeoDim]
        print('----- Default attributes -----')
        for k, v in Defaults.items():
            print('--- {} ---'.format(k))
            for k1, v1 in v.items():
                print(' - {} : {} ({})'.format(k1, v1, v1.__class__.__name__))


    def LoadElemsAttr(self):
        ''' Calculate properties of triagular elements.
        Triangle: N x 3 x 3 numpy array.
        '''

        self.SetDefaultsAttr()
        self.SetInputPoints()

        A = self.Triangles[:, 0, :].squeeze()
        B = self.Triangles[:, 1, :].squeeze()
        C = self.Triangles[:, 2, :].squeeze()

        AB = B-A
        AC = C-A
        BC = C-B
        BA = A-B
        CA = -AC
        CB = -BC

        dAB = Norm(AB)
        dBC = Norm(BC)
        dAC = Norm(AC)

        s = (dAB+dBC+dAC)/2
        self.GeomInput['area'] = np.sqrt(s*(s-dAB)*(s-dBC)*(s-dAC))
        self.normalVec = Cross(AB, AC)
        self.GeomInput['a'] = self.normalVec[:, 0]
        self.GeomInput['b'] = self.normalVec[:, 1]
        self.GeomInput['c'] = self.normalVec[:, 2]
        self.GeomInput['d'] = -Dot(self.normalVec, A)
        self.GeomInput['plane_norm'] = Norm(self.normalVec)

        self.Centroid = 1/3*(np.sum(self.Triangles, 1)).squeeze()
        # rr = sqrt(centroid(1).^2 + centroid(2).^2);
        # density(i) = interpn(y,z,dens,rr,centroid(3));
        # Set defaults parameters
        print('Attributes of triangular elements computed and stored in GeomInput ...')


    def SetInputPoints(self, GeoDim='3D'):
        for j, s in enumerate(['x', 'y', 'z']):
            for i in range(3):
                if self.Verbose: print( 'setting "{}{}"'.format(s, i+1))
                self.GeomInput['{}{}'.format(s, i+1)] = self.Triangles[:, i, j].squeeze()

    def SetAttr(self,Attribute,Value):
        if self.GeomInput.get(Attribute) is None:
                print('Scalar attribute "{}" not found in GeomInput ... Creating a new scalar')
                self.GeomInput[Attribute] = Value
        else:
            self.GeomInput[Attribute] = Value


    def SetElemAttr(self, GroupID, Attribute, Value):
        assert self.NElems > 0, 'Number of elements must be >0 to set \
                                defaults input values. Was a mesh imported?'
        if self.GeomInput.get(Attribute) is None:
            if GroupID == []:
                print('Element attribute "{}" not found in GeomInput ... Creating a new array for all elements')
                self.GeomInput[Attribute] = np.full((self.NElems), Value)
            else:
                raise ValueError('Cannot create a new attribute for specific groups. GroupID must be []')
        else:
            self.GeomInput[Attribute][self.GetGroupIdx(GroupID)] = Value




class GeomPlot(GeomGroup):
    '''Subclass for plotting methods.'''

    def PlotTriangles(self, GroupID=None, ax=None, cmap_name='viridis', ElemAttr=None, Alpha=0.1, EdgeColor='k', FaceColor='b'):
        if ax is None:
            ax = self.ax
        if ax is None:
            ax = plt.gca()
        cmap = mpl.cm.get_cmap(cmap_name)               # Get colormap by name
        # c = cmap(mpl.colors.Normalize(vmin, vmax)(v))   # Normalize value and get color
        Idx = self.GetGroupIdx(GroupID)
        # Create PolyCollection from coords
        pc = Poly3DCollection(
            self.Triangles[Idx, :, :], cmap=mpl.cm.jet, alpha=0.4)

        # Color set by values of element attribute if ElemAttr not None (e.g. ElemAttr='Z')
        if ElemAttr is not None and self.GeomInput.get(ElemAttr) is not None:
            pc.set_array(self.GeomInput.get(ElemAttr))
        else:

            pc.set_facecolor(FaceColor)
        pc.set_edgecolor(EdgeColor)
        pc.set_alpha(Alpha)

        # Add PolyCollection to axes
        ax.add_collection3d(pc)
        return pc

    def Plot(self, GroupID=None, ax=None, fig=None, **kwargs):
        if fig is None:
            self.fig = plt.figure()
        if ax is None:
            self.ax = self.fig.add_subplot(111, projection='3d')
        else:
            self.ax = ax
        self.PlotTriangles(GroupID, ax, **kwargs)
        xmin = np.min(self.Triangles[:, :, 0])
        xmax = np.max(self.Triangles[:, :, 0])
        ymin = np.min(self.Triangles[:, :, 1])
        ymax = np.max(self.Triangles[:, :, 1])
        zmin = np.min(self.Triangles[:, :, 2])
        zmax = np.max(self.Triangles[:, :, 2])
        mn = min([xmin, ymin, zmin])
        mx = max([xmax, ymax, zmax])
        self.SetAxisLim(mn, mx)

    def ShowCentroids(self, ax=None):
        if ax is None:
            ax = self.ax
        if ax is None:
            ax = plt.gca()
        ax.scatter(self.Centroid[:, 0], self.Centroid[:, 1],
                   self.Centroid[:, 2], marker='o', color='b')

    def ShowNormals(self, GroupID=None, ax=None, L=0.002, Color='b'):
        if ax is None:
            ax = self.ax
        if ax is None:
            ax = plt.gca()
        c = self.Centroid
        v = self.normalVec

        Idx = self.GetGroupIdx(GroupID)
        if self.Verbose:
            print('Normals Idx:', Idx)
        ax.quiver(c[Idx, 0], c[Idx, 1], c[Idx, 2], v[Idx, 0],
                  v[Idx, 1], v[Idx, 2], length=L, normalize=True, color=Color)
        # ax.add_collection(lc)

    def SetAxisLim(self, mn, mx):
        self.ax.set_xlim3d(mn, mx)
        self.ax.set_ylim3d(mn, mx)
        self.ax.set_zlim3d(mn, mx)

class GeomSetup(GeomInput, GeomPlot):

    def __init__(self, mesh=None, ax=None, Verbose=False):
        super().__init__()  # init the parent class as well to  set attributes
        self.ax = ax
        self.Verbose = Verbose
        if mesh is not None:
            self.ImportMesh(mesh)

    def ImportMesh(self, mesh):
        mesh = os.path.abspath(mesh)
        print('Loading {} ...'.format(mesh))
        assert os.path.exists(mesh), "Cannot read '{}'".format(mesh)

        try:
            gmsh.finalize()
        except:
            pass
        finally:
            gmsh.initialize()

        gmsh.open(mesh)

        mesh = gmsh.model.mesh
        self.Model = gmsh.model
        print('Model ' + gmsh.model.getCurrent() +
              ' (' + str(gmsh.model.getDimension()) + 'D)')
        self.GetGroups()
        self.ImportMeshElements(mesh)
        self.LoadElemsAttr()

    def ImportMeshElements(self, mesh) -> None:

        self.Mesh = mesh
        elem_types, elem_tags, node_tags = mesh.getElements(dim=2)
        if self.Verbose:
            print('elem_types:{}\n elem_tags:{} \n node_tags:{}'.format(
                elem_types, elem_tags, node_tags))
        idx, points, param = mesh.getNodes(2, -1, True)
        self.Vertices = np.asarray(points).reshape(-1, 3)
        self.VerticeTags = np.asarray(node_tags[0]).reshape(-1, 3)
        self.ElementsTags = np.asarray(elem_tags[0])
        self.idxVertices = np.zeros((int(np.max(idx))+1,), dtype=int)
        self.idxElements = np.zeros(
            (int(np.max(self.ElementsTags))+1,), dtype=int)
        for i in range(np.size(idx)):
            self.idxVertices[idx[i]] = int(i)
        for i in range(np.size(self.ElementsTags)):
            self.idxElements[self.ElementsTags[i]] = int(i)

        self.Triangles = self.Vertices[self.idxVertices[self.VerticeTags], :]
        self.NElems = self.Triangles.shape[0]
        assert len(self.Triangles.shape) > 1 and self.Triangles.shape[1] == 3 and self.Triangles.shape[
            2] == 3, "Points must be a Nx3x3 numpy arrays: N_triangle x 3_vertices x (x_vertice,y_vertice,z_vertice)"
        print('Triangular mesh elements imported ... Number of triangular elements: {}'.format(
            self.Triangles.shape[0]))




    def CheckGeomInput(self, GeoDim='3D'):

        Defaults = self.DefaultsAttr[GeoDim]
        Requirements = list(Defaults['Vector'].keys()) + list(Defaults['Scalar'].keys())
        for Field in Requirements:
            if self.GeomInput.get(Field) is None:
                raise ValueError('No field "{}" found to generate geometry cfg file. \n Requirements for "{}" cfg file: {}'.format(
                    Field, GeoDim, Requirements))

        Vec = list(Defaults['Vector'].keys())
        S = self.GeomInput[Vec[0]].shape
        for (k, v) in self.GeomInput.items():
            if k in Vec:
                if v.shape != S:
                    raise ValueError(
                        "Issue with shape of field '{}':{}.shape={} ; {}.shape={}".format(k, k, v.shape, Vec[0], S))

    def ConvertGeomInput(self):
        D = {}
        for (k,v) in self.GeomInput.items():
            if type(v) == np.ndarray:
                D[k] = v.tolist()
            else:
                D[k] = v
        return {'geom':D}


    def WriteGeomFile(self, FileName='gitrGeom.cfg', Folder='',GeoDim='3D', OverWrite=False):

        FileName = os.path.abspath(os.path.join(Folder,FileName))
        print('Writing geometry config into {} ...'.format(FileName))

        self.CheckGeomInput()

        if os.path.exists(FileName):
            if not OverWrite and not click.confirm('File {} exists.\n Do you want to overwrite it?'.format(FileName), default=True):
                return

        with io.open(FileName,'w') as f:
            libconf.dump(self.ConvertGeomInput(),f)




