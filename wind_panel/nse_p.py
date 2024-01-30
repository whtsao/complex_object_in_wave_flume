"""
Incompressible Navier-Stokes flow around a cylinder in 2D.
"""
from __future__ import absolute_import
from __future__ import division
from builtins import object
from proteus import *
from proteus.default_p import *
import sys

import cell
from cell import *
from proteus.mprans import RANS2P

name="cell"

nd = cell.nd
domain = cell.domain

bcsTimeDependent = True
LevelModelType = RANS2P.LevelModel

coefficients = RANS2P.Coefficients(rho_0=cell.rho,
                                   nu_0=cell.nu,
                                   rho_1=cell.rho,
                                   nu_1=cell.nu,
                                   g=cell.g,
                                   nd=cell.nd,
                                   forceStrongDirichlet=False,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=100.0,
                                   useRBLES=0.0,
                                   useMetrics=1.0,
                                   particle_sdfList = [cell.sdf_vectorized_stl if cell.opts.stl else cell.sdf_vectorized],
                                   #particle_sdfList = [cell.sdf_vectorized_stl],
                                   particle_velocityList = [lambda t,x: (0.0,0.0,0.0)],
                                   use_ball_as_particle=0,
                                   nParticles=1,
                                   useExact=True)

coefficients.projection_direction=np.array([1.0,0.0,0.0])


def vel(x,t):
    return t*2.

def getDBC_p(x,flag):
    if flag == cell.boundaryTags['downstream']:
        return lambda x,t: (x[2] - cell.x[2] - cell.L[2])*cell.rho*cell.g[2]
    
def getDBC_u(x,flag):
    if flag == cell.boundaryTags['upstream']:
        return lambda x,t: vel(x,t)
#    elif flag == cell.boundaryTags['downstream']:
#        return lambda x,t: 0.0 #only for upwinding
def getDBC_v(x,flag):
    if flag == cell.boundaryTags['upstream']:
        return lambda x,t: 0.0
#    elif flag == cell.boundaryTags['downstream']:
#        return lambda x,t: 0.0
    
def getDBC_w(x,flag):
    if flag == cell.boundaryTags['upstream']:
        return lambda x,t: 0.0
#    elif flag == cell.boundaryTags['downstream']:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag == cell.boundaryTags['upstream']:
        return lambda x,t: -vel(x,t)
    elif flag == cell.boundaryTags['downstream']:
        return None
    else:
        return lambda x,t: 0.0
                          

def getAFBC_u(x,flag):
    if flag in [cell.boundaryTags['upstream'],cell.boundaryTags['downstream']]:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag in [cell.boundaryTags['upstream'],cell.boundaryTags['downstream']]:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag in [cell.boundaryTags['upstream'],cell.boundaryTags['downstream']]:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if flag == cell.boundaryTags['upstream']:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag == cell.boundaryTags['upstream']:#cell.boundaryTags['downstream']]:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    if flag == cell.boundaryTags['upstream']:#,cell.boundaryTags['downstream']]:
        return None
    else:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class Steady_p(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return (x[2] - cell.x[2] - cell.L[2])*cell.rho*cell.g[2]

class Steady_u(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return vel(x,t)

class Steady_v(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_w(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0



initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}
