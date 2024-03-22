#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 13 05:09:10 2023
Flow around obstacles. BR Riverdomain. 2d NS

@author: linojr
"""
from math import *
from proteus import (Domain,
                     Context,
                     Profiling)
from proteus.default_n import *                     

opts = Context.Options([
    ("T", 100.0, "Time interval [0, T]"),
    ("he",0.05, "maximum size of edges"),
    ("onlySaveFinalSolution",False,"Only save the final solution"),
    ("usePETSc",False,"use parallel solvers from PETSc"),
    ("cell_bottom",0.0,"the starting elevation of the flow cell"),
    ("cell_height",0.05,"the height of the flow cell"),
    ("stl", True, "use an stl for the structure")
], mutable=True)


usePETSc = opts.usePETSc


boundaryTags = {'upstream' : 1, 'wall' : 2, 'downstream' : 3, 'top': 4, 'bottom': 5, 'holes': 6}

nd = 3
if opts.stl:
    x = (-2.5, -2.5, opts.cell_bottom)
    L = (5.0, 5.0, opts.cell_height)
else:
    x = (0.0, 0.0, opts.cell_bottom)
    L = (1.0, 1.0, opts.cell_height)


vertices=[[x[0]     , x[1]     , x[2]],#0
          [x[0]+L[0], x[1]     , x[2]],#1
          [x[0]+L[0], x[1]+L[1], x[2]],#2
          [x[0]     , x[1]+L[1], x[2]],#3
          [x[0]     , x[1]     , x[2] + L[2]],#4
          [x[0]+L[0], x[1]     , x[2] + L[2]],#5
          [x[0]+L[0], x[1]+L[1], x[2]+L[2]],#6
          [x[0]     , x[1]+L[1], x[2]+L[2]]]#7
vertexFlags=[boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall'],
             boundaryTags['wall']]
facets=[[[0,1,2,3]],
        [[0,1,5,4]],
        [[1,2,6,5]],
        [[2,3,7,6]],
        [[3,0,4,7]],
        [[4,5,6,7]]]
facetFlags=[boundaryTags['bottom'],
            boundaryTags['wall'],
            boundaryTags['downstream'],
            boundaryTags['wall'],
            boundaryTags['upstream'],
            boundaryTags['top']]
volumes=[[[0,1,2,3,4,5]]]
regions=[[0.5*L[0],0.5*L[1],x[2]+0.5*L[2]]]
regionFlags=[0]

domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             regions = regions,
                                             regionFlags = regionFlags)
#                                             holes=holes)
domain.holes_ind=[]
domain.volumes=volumes
nLevels = 1
parallelPartitioningType = MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
domain.MeshOptions.setParallelPartitioningType('node')
#domain.MeshOptions.use_gmsh=True
domain.MeshOptions.he=opts.he
domain.MeshOptions.triangleFlag=0
domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((opts.he**3)/6.0,)
domain.boundaryTags = boundaryTags
domain.writePoly("cell")

#cutfem structure
from mangrove import sdf_vectorized
from mangrove import sdf_vectorized_stl

# Time stepping
T= opts.T
runCFL = 0.33
dt_fixed = 0.1
dt_init = 0.001
nDTout = int(T/dt_fixed)
dt_init = min(dt_init,0.5*dt_fixed)
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)] 
if opts.onlySaveFinalSolution == True:
    tnList = [0.0,dt_init,opts.T]


# Numerical parameters
ns_shockCapturingFactor  = 0.0
ns_lag_shockCapturing = True
ns_lag_subgridError = True

epsFact_density    = 1.5
epsFact_viscosity  = 1.5
epsFact_redistance = 0.33
epsFact_consrv_heaviside = 1.5
epsFact_consrv_dirac     = 1.5
epsFact_consrv_diffusion = 10.0

# Fluid
rho = 998.2
nu = 1.004e-6


# Gravity
g = [0.0,0.0,-9.81]

domain.MeshOptions.triangleOptions="VApq1.35q12feena{0:1.16f}".format(opts.he**3 / 6.0)
domain.MeshOptions.genMesh=True
