from __future__ import absolute_import
from proteus import *
from proteus.default_n import *

import cell
from cell import *
from nse_p import *

systemStepExact=False

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
#stepController = HeuristicNL_dt_controller
#nonlinearIterationsFloor = 2
#nonlinearIterationsCeil=4
#nonlinearIterationsFloor = 3
#nonlinearIterationsCeil=4
#dtNLgrowFactor  = 1.5
#dtNLreduceFactor= 0.75

basis=C0_AffineLinearOnSimplexWithNodalBasis
elementQuadrature = SimplexGaussQuadrature(cell.nd,5)
elementBoundaryQuadrature = SimplexGaussQuadrature(cell.nd-1,5)

femSpaces = {0:basis,
             1:basis,
             2:basis,
             3:basis}

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,cell.nd,lag=True,hFactor=1.0)
shockCapturing = RANS2P.ShockCapturing(coefficients,cell.nd,0.0,lag=True)

massLumping = False

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None

matrix = SparseMatrix

if cell.usePETSc:    
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
    linearSmoother = SimpleNavierStokes3D
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

tolFac = 0.0
linTolFac = 0.01
nl_atol_res = 1.0e-6
l_atol_res = 1.0e-8

maxNonlinearIts = 100
maxLineSearches =0
#conservativeFlux = {0:'pwl-bdm-opt'}
