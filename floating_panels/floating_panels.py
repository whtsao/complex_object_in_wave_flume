import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import proteus.TwoPhaseFlow.utils.Parameters as Parameters
from proteus import WaveTools as wt
import math

# dependencies for FSI
from proteus.mbd import CouplingFSI as fsi
import pychrono

#class spring_function(pychrono.ForceFunctorP):

#    def __call__(self, d_time, d_rest_length, d_length, d_vel, o_spring_object):
#        """
#        Overrides the spring calculation class with a custom calculation function. This current sets any value over the
#        spring maximum as the maximum to mnimize the the affect of the broken flexor and the simulation.

#        Parameters
#        ----------
#        d_time: float
#            Seconds
#        d_rest_length: float
#            meters
#        d_length: float
#            metters
#        d_vel: float
#            meeers per second
#        o_spring_object: objecti
#            Sping chrono object

#        Returns
#        -------
#        d_force: float
#            netwons

#        """

#        d_spring_constant = 4000.0 # [N/m] Spring force after initial tolerance

#        d_force = -d_spring_constant * (d_length-d_rest_length)

#        return d_force

opts= Context.Options([
    ("final_time",60.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.5,"Desired CFL restriction"),
    ("he",0.05,"he relative to Length of domain in x"),
    ("wave_height",0.05,"Wave height"),
    ("wave_period",2.,"Wave height"),
    ("mwl", 1., "still water line"),
    ("tank_length", 10., "tank length"),
    ("mooring",True,"True if the mooring lines are attached"),
    ("collision",False,"True if the mooring lines is collision body"),
    ("ic_angle",0.,"Initial pitch angle of the floating platform (deg)"),
    ("bodybool1",True,"add floating body if True"),
    ("bodybool2",True,"add floating body if True"),
    ("linkbool",True,"add floating body if True"),
    ])

# general options
# sim time
T = opts.final_time
# initial step
dt_init = 0.001
# CFL value
cfl = opts.cfl
# mesh size
he = opts.he
# rate at which values are recorded
sampleRate = 0.05


# for ALE formulation
movingDomain = True
# for added mass stabilization
addedMass = True

# physical options
# water density
rho_0 = 998.2
# water kinematic viscosity
nu_0 = 1.004e-6
# air density
rho_1 = 1.205
# air kinematic viscosity
nu_1 = 1.5e-5
# gravitational acceleration
g = np.array([0., -9.81, 0.])

# body options
fixed = False

# wave channel
water_level = opts.mwl
tank_length = opts.tank_length

# wave options
wave_period = opts.wave_period
wave_height = opts.wave_height
wave_direction = np.array([1., 0., 0.])
wave_type = 'Fenton'  #'Linear'
Nf = 8 # number of Fourier coefficients
wave = wt.MonochromaticWaves(period=wave_period,
                             waveHeight=wave_height,
                             mwl=water_level,
                             depth=water_level,
                             g=g,
                             waveDir=wave_direction,
                             waveType=wave_type,
                             Nf=8)
wavelength = wave.wavelength

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

# TANK
tank = st.Tank2D(domain, dim=(tank_length, 1.5*water_level))

# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
tank.setSponge(x_n=wavelength, x_p=wavelength)

# Panel dimensions and mass
body_w1 = 1.
body_h1 = 0.1

body_w2 = 1.
body_h2 = 0.1

gap = 0.1

#  --------w1--------       ---------w2------
#  |                |       |               |
#  |                h1=link=|               h2
#  --------w1--------       ---------w2------

thob = 400.
mb1 = thob*body_w1*body_h1
ib1 = mb1*(body_w1**2+body_h1**2)/12.
mb2 = thob*body_w2*body_h2
ib2 = mb2*(body_w2**2+body_h2**2)/12.

# initial displacement and angle
yst1 = mb1/body_w1/rho_0
yst2 = mb2/body_w2/rho_0
ic_angle = (opts.ic_angle/180.)*math.pi
d_front = 1.


# stiffness and damping coefficient of the link between panels
ki = 1000.
ci = 0.


if opts.bodybool1:
    # PANEL 1
    caisson1 = st.Rectangle(domain, dim=(body_w1, body_h1), coords=(0., 0.))
    # set barycenter in middle of caisson
    caisson1.setBarycenter([0., 0.])
    # caisson is considered a hole in the mesh
    caisson1.setHoles([[0., 0.]])
    # 2 following lines only for py2gmsh
    caisson1.holes_ind = np.array([0])
    tank.setChildShape(caisson1, 0)
    # translate caisson to middle of the tank
    caisson1.translate(np.array([d_front+0.5*body_w1, water_level-(yst1-0.5*body_h1)]))
    caisson1.rotate(rot = ic_angle)

if opts.bodybool2:
    # PANEL 2
    caisson2 = st.Rectangle(domain, dim=(body_w2, body_h2), coords=(0., 0.))
    # set barycenter in middle of caisson
    caisson2.setBarycenter([0., 0.])
    # caisson is considered a hole in the mesh
    caisson2.setHoles([[0., 0.]])
    # 2 following lines only for py2gmsh
    caisson2.holes_ind = np.array([0])
    tank.setChildShape(caisson2, 0)
    # translate caisson to middle of the tank
    caisson2.translate(np.array([d_front+body_w1+gap+0.5*body_w2, water_level-(yst2-0.5*body_h2)]))
    caisson2.rotate(rot = ic_angle)


#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

# SYSTEM

# create system
system = fsi.ProtChSystem()
# access chrono object
chsystem = system.getChronoObject()
# communicate gravity to system
# can also be set with:
# system.ChSystem.Set_G_acc(pychrono.ChVectorD(g[0], g[1], g[2]))
system.setGravitationalAcceleration(g)
# set maximum time step for system
system.setTimeStep(1e-4)

solver = pychrono.ChSolverMINRES()
chsystem.SetSolver(solver)

# FLOATING PANEL 1
if opts.bodybool1:
    # create floating body
    body = fsi.ProtChBody(system=system)
    # give it a name
    body.setName(b'my_body1')
    # attach shape: this automatically adds a body at the barycenter of the caisson shape
    body.attachShape(caisson1)
    # set 2D width (for force calculation)
    body.setWidth2D(1.)
    # access chrono object
    chbody = body.getChronoObject()
    # impose constraints
    chbody.SetBodyFixed(fixed)
    free_x = np.array([1., 1., 0.]) # translational
    free_r = np.array([0., 0., 1.]) # rotational
    body.setConstraints(free_x=free_x, free_r=free_r)
    # access pychrono ChBody
    # set mass
    # body.ChBody.SetMass(14.5)
    body.setMass(mb1)
    # set inertia
    # body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
    body.setInertiaXX(np.array([1., 1., ib1]))
    # record values
    body.setRecordValues(all_values=True)

# FLOATING PANEL 2
if opts.bodybool2:
    # create attached body
    body = fsi.ProtChBody(system=system)
    # give it a name
    body.setName(b'my_body2')
    # attach shape: this automatically adds a body at the barycenter of the caisson shape
    body.attachShape(caisson2)
    # set 2D width (for force calculation)
    body.setWidth2D(1.)
    # access chrono object
    chbody = body.getChronoObject()
    # impose constraints
    chbody.SetBodyFixed(fixed)
    free_x = np.array([1., 1., 0.]) # translational
    free_r = np.array([0., 0., 1.]) # rotational
    body.setConstraints(free_x=free_x, free_r=free_r)
    # access pychrono ChBody
    # set mass
    # body.ChBody.SetMass(14.5)
    body.setMass(mb2)
    # set inertia
    # body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
    body.setInertiaXX(np.array([1., 1., ib2]))
    # record values
    body.setRecordValues(all_values=True)

# LINK (connect two panels)
if opts.linkbool:
    TSDA1 = pychrono.ChLinkTSDA()
    body1_point = pychrono.ChVectorD(d_front+body_w1,water_level-(yst1-0.5*body_h1),0.0)
    body2_point = pychrono.ChVectorD(d_front+body_w1+gap,water_level-(yst2-0.5*body_h2),0.0)
    TSDA1.Initialize(system.subcomponents[0].ChBody,
                                    system.subcomponents[1].ChBody,
                                    False, body1_point, body2_point, auto_rest_length=True)
    #o_spring_force_functor = spring_function()
    #TSDA.RegisterForceFunctor(o_spring_force_functor)
    TSDA1.SetSpringCoefficient(ki)
    TSDA1.SetDampingCoefficient(ci)
    system.ChSystem.Add(TSDA1)


# MOORING
if opts.mooring:
    # length
    L = (water_level**2+d_front**2)**0.5 # m
    # submerged weight
    w = 0.0778  # kg/m
    # equivalent diameter (chain -> cylinder)
    d = 2.5e-3 # m
    # unstretched cross-sectional area
    A0 = (np.pi*d**2/4.)
    # density
    dens = w/A0+rho_0
    # number of elements for cable
    nb_elems = 50
    # Young's modulus
    E = (1.e10)/50**3/A0
    #E = (753.6e6)/50**3/A0
    #E = 1.e8 #5.44e10

    # fairleads coordinates
    fairlead = np.array([d_front, water_level-(yst1-0.5*body_h1), 0.])
    #fairlead1 = np.array([0.5*tank_length-0.5*body_w1, y0+body_h2, 0.])
    #fairlead2 = np.array([0.5*tank_length+0.5*body_w1, y0+body_h2, 0.])

    # anchors coordinates
    anchor = np.array([0.1, 0., 0.])
    #anchor1 = np.array([fairlead[0]-lx, 0., 0.])
    #anchor2 = np.array([fairlead[0]+lx, 0., 0.])
    #anchor1 = np.array([fairlead1[0]-lx, 0., 0.])
    #anchor2 = np.array([fairlead2[0]+lx, 0., 0.])

    # quasi-statics for finding shape of cable
    from pycatenary.cable import MooringLine
    # create lines
    EA = E*A0
    cat1 = MooringLine(L=L,
                    w=w*9.81,
                    EA=EA,
                    anchor=anchor,
                    fairlead=fairlead,
                    nd=2,
                    floor=True)
    
#    cat2 = MooringLine(L=L,
#                    w=w*9.81,
#                    EA=EA,
#                    anchor=anchor2,
#                    fairlead=fairlead2,
#                    nd=2,
#                    floor=True)

    cat1.computeSolution()
#    cat2.computeSolution()

    # ANCHOR
    # arbitrary body fixed in space
    body1 = fsi.ProtChBody(system)
    body1.barycenter0 = anchor #np.zeros(3)
    # fix anchor in space
    body1.ChBody.SetBodyFixed(True)

    # arbitrary body fixed in space
#    body2 = fsi.ProtChBody(system)
#    body2.barycenter0 = anchor2 #np.zeros(3)
    # fix anchor in space
#    body2.ChBody.SetBodyFixed(True)

    # MESH
    # initialize mesh that will be used for cables
    mesh = fsi.ProtChMesh(system)

    # FEM CABLES
    # moorings line 1
    m1 = fsi.ProtChMoorings(system=system,
                            mesh=mesh,
                            length=np.array([L]),
                            nb_elems=np.array([nb_elems]),
                            d=np.array([d]),
                            rho=np.array([dens]),
                            E=np.array([E]))
    m1.setName(b'mooring1')
    # send position functions from catenary to FEM cable
    m1.setNodesPositionFunction(cat1.s2xyz, cat1.ds2xyz)
    # sets node positions of the cable
    m1.setNodesPosition()
    # build cable
    m1.buildNodes()
    # apply external forces
    m1.setApplyDrag(True)
    m1.setApplyBuoyancy(True)
    m1.setApplyAddedMass(True)
    # set fluid density at cable nodes
    m1.setFluidDensityAtNodes(np.array([rho_0 for i in range(m1.nodes_nb)]))
    # sets drag coefficients
    m1.setDragCoefficients(tangential=1.15, normal=0.213, segment_nb=0)
    # sets added mass coefficients
    m1.setAddedMassCoefficients(tangential=0.269, normal=0.865, segment_nb=0)
    # small Iyy for bending
    m1.setIyy(0., 0)
    # attach back node of cable to body
    m1.attachBackNodeToBody(body)
    # attach front node to anchor
    m1.attachFrontNodeToBody(body1)

    # mooring line 2
#    m2 = fsi.ProtChMoorings(system=system,
#                            mesh=mesh,
#                            length=np.array([L]),
#                            nb_elems=np.array([nb_elems]),
#                            d=np.array([d]),
#                            rho=np.array([dens]),
#                            E=np.array([E]))
#    m2.setName(b'mooring2')
    # send position functions from catenary to FEM cable
#    m2.setNodesPositionFunction(cat2.s2xyz, cat2.ds2xyz)
    # sets node positions of the cable
#    m2.setNodesPosition()
    # build cable
#    m2.buildNodes()
    # apply external forces
#    m2.setApplyDrag(True)
#    m2.setApplyBuoyancy(True)
#    m2.setApplyAddedMass(True)
    # set fluid density at cable nodes
#    m2.setFluidDensityAtNodes(np.array([rho_0 for i in range(m2.nodes_nb)]))
    # sets drag coefficients
#    m2.setDragCoefficients(tangential=1.15, normal=0.213, segment_nb=0)
    # sets added mass coefficients
#    m2.setAddedMassCoefficients(tangential=0.269, normal=0.865, segment_nb=0)
    # small Iyy for bending
#    m2.setIyy(0., 0)
    # attach back node of cable to body
#    m2.attachBackNodeToBody(body)
    # attach front node to anchor
#    m2.attachFrontNodeToBody(body2)

    if opts.collision:
        # CONTACT MATERIAL
        # define contact material for collision detection
        material = pychrono.ChMaterialSurfaceSMC()
        material.SetKn(3e6)  # normal stiffness
        material.SetGn(1.)  # normal damping coefficient
        material.SetFriction(0.3)
        material.SetRestitution(0.2)
        material.SetAdhesion(0)

        # SEABED
        # create a box
        #seabed = pychrono.ChBodyEasyBox(100., 0.2, 1., 1000, True)
        seabed = pychrono.ChBodyEasyBox(100., 0.2, 1., 1000, True, True, material)
        # move box
        seabed.SetPos(pychrono.ChVectorD(0., -0.1-2.*d, 0.))
        # fix boxed in space
        seabed.SetBodyFixed(True)
        # add box to system
        system.ChSystem.Add(seabed)

        # add material to objects
        #seabed.SetMaterialSurface(material) # not valid for new Chrono
        m1.setContactMaterial(material)
        #m2.setContactMaterial(material)

#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# CAISSON

# set no-slip conditions on caisson
if opts.bodybool1:
    for tag, bc in caisson1.BC.items():
        bc.setNoSlip()

if opts.bodybool2:
    for tag, bc in caisson2.BC.items():
        bc.setNoSlip()

# TANK

# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# non material boundaries for sponge interface
tank.BC['sponge'].setNonMaterial()

# fix in space nodes on the boundaries of the tank
for tag, bc in tank.BC.items():
    bc.setFixedNodes()

# WAVE AND RELAXATION ZONES

smoothing = he*1.5
dragAlpha = 5*2*np.pi/wave_period/(1.004e-6)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.setGenerationZones(x_n=True,
                        waves=wave,
                        smoothing=smoothing,
                        dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_p=True,
                        dragAlpha=dragAlpha)



#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
smoothing = 1.5*he
nd = domain.nd

class P_IC:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = tank.dim[nd-1] - water_level
        phi = x[nd-1] - water_level
        p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
        return p
class zero:
    def uOfXT(self, x, t):
        return 0.0
class U_IC:
    def uOfXT(self, x, t):
        return 0.0
class V_IC:
    def uOfXT(self, x, t):
        return 0.0
class W_IC:
    def uOfXT(self, x, t):
        return 0.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(smoothing, x[nd-1]-water_level)
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.genMesh = True
domain.MeshOptions.he = he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)
#
st.assembleDomain(domain)




#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()

myTpFlowProblem.outputStepping.final_time = T
myTpFlowProblem.outputStepping.dt_output=sampleRate
myTpFlowProblem.outputStepping.dt_init=dt_init
myTpFlowProblem.domain = domain

myTpFlowProblem.SystemNumerics.cfl = cfl
myTpFlowProblem.SystemNumerics.useSuperlu=True #False


#outputStepping = TpFlow.OutputStepping(
#    final_time=T,
#    dt_init=dt_init,
#    # cfl=opts.cfl,
#    dt_output=sampleRate,
#    nDTout=None,
#    dt_fixed=None,
#)

#myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
#    ns_model=None,
#    ls_model=None,
#    nd=domain.nd,
#    cfl=cfl,
#    outputStepping=outputStepping,
#    structured=False,
#    he=he,
#    domain=domain,
#    initialConditions=initialConditions,
#    useSuperlu=True
#)


# Necessary for moving domains
#myTpFlowProblem.movingDomain = movingDomain
#params = myTpFlowProblem.Parameters

myTpFlowProblem.SystemPhysics.movingDomain = movingDomain
params = myTpFlowProblem.SystemPhysics

# PHYSICAL PARAMETERS
params['rho_0'] = rho_0  # water
params['rho_1'] = rho_1 # air
params['nu_0'] = nu_0  # water
params['nu_1'] = nu_1  # air
params['gravity'] = g
params['surf_tension_coeff'] = 0.0

#params.physical.densityA = rho_0  # water
#params.physical.densityB = rho_1  # air
#params.physical.kinematicViscosityA = nu_0  # water
#params.physical.kinematicViscosityB = nu_1  # air
#params.physical.gravity = g
#params.physical.surf_tension_coeff = 0.

params.addModel(Parameters.ParametersModelMoveMeshElastic,'move')
params.useDefaultModels()
# added mass estimation
if addedMass is True:
    params.addModel(Parameters.ParametersModelAddedMass,'addedMass')

m = myTpFlowProblem.SystemPhysics.modelDict

m['move'].p.initialConditions['hx'] = zero()
m['move'].p.initialConditions['hy'] = zero()
m['flow'].p.initialConditions['p'] = zero()
m['flow'].p.initialConditions['u'] = zero()
m['flow'].p.initialConditions['v'] = zero()
m['vof'].p.initialConditions['vof'] = VF_IC()
m['ncls'].p.initialConditions['phi'] = PHI_IC()
m['rdls'].p.initialConditions['phid'] = PHI_IC()
m['mcorr'].p.initialConditions['phiCorr'] = zero()
m['addedMass'].p.initialConditions['addedMass'] = zero()

# add for proteus/fct
m['flow'].p.coefficients.projection_direction = np.array([1.0, 0.0, 0.0])

# ADD RELAXATION ZONES TO AUXILIARY VARIABLES
m['flow'].auxiliaryVariables += domain.auxiliaryVariables['twp']
# ADD SYSTEM TO AUXILIARY VARIABLES
m['flow'].auxiliaryVariables += [system]
m['flow'].p.coefficients.NONCONSERVATIVE_FORM=0
if addedMass is True:
    # passed in added_mass_p.py coefficients
    m['addedMass'].auxiliaryVariables += [system.ProtChAddedMass]
    max_flag = 0
    max_flag = max(domain.vertexFlags)
    max_flag = max(domain.segmentFlags+[max_flag])
    max_flag = max(domain.facetFlags+[max_flag])
    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
    for s in system.subcomponents:
        if type(s) is fsi.ProtChBody:
            for i in s.boundaryFlags:
                flags_rigidbody[i] = 1
    m['addedMass'].p.coefficients.flags_rigidbody = flags_rigidbody


#m = myTpFlowProblem.Parameters.Models

# MODEL PARAMETERS
#ind = -1
# first model is mesh motion (if any)
#if movingDomain:
#    m.moveMeshElastic.index = ind+1
#    ind += 1
# navier-stokes
#m.rans2p.index = ind+1
#ind += 1
# volume of fluid
#m.vof.index = ind+1
#ind += 1
# level set
#m.ncls.index = ind+1
#ind += 1
# redistancing
#m.rdls.index = ind+1
#ind += 1
# mass correction
#m.mcorr.index = ind+1
#ind += 1
# added mass estimation
#if addedMass is True:
#    m.addedMass.index = ind+1
#    ind += 1

# ADD RELAXATION ZONES TO AUXILIARY VARIABLES
#m.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
# ADD SYSTEM TO AUXILIARY VARIABLES
#m.rans2p.auxiliaryVariables += [system]
#m.rans2p.p.coefficients.NONCONSERVATIVE_FORM=0
#if addedMass is True:
#    # passed in added_mass_p.py coefficients
#    m.addedMass.auxiliaryVariables += [system.ProtChAddedMass]
#    max_flag = 0
#    max_flag = max(domain.vertexFlags)
#    max_flag = max(domain.segmentFlags+[max_flag])
#    max_flag = max(domain.facetFlags+[max_flag])
#    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
#    for s in system.subcomponents:
#        if type(s) is fsi.ProtChBody:
#            for i in s.boundaryFlags:
#                flags_rigidbody[i] = 1
#    m.addedMass.p.coefficients.flags_rigidbody = flags_rigidbody
