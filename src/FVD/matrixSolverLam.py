#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = "Ananda S. Kannan"
__copyright__ = "Copyright 2017, pyFy project"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ananda S. Kannan"
__email__ = "ansubru@gmail.com"

"""
################################################################################ LAMINAR FLOW SOLVER ################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

from Discretize2 import Discretize2
disc_obj2 = Discretize2()

from gaussSeidel2 import gaussSiedel2
gs_obj = gaussSiedel2()

from solFunc import solFunc
solFunc_obj = solFunc()

from IO import IO
IO_obj = IO("random")

from Corr2 import Corr2
corr_obj = Corr2()

from plotResults import plotResults
plt_obj = plotResults()

from rhieChow2 import rhieChow2
rc_obj = rhieChow2()

from calcResiduals import calcResiduals
res_obj = calcResiduals()

################################################################ INPUT TO SOLVER ###########################################################################################################################

grid = IO_obj.grid_nodes #generate a base grid node layout
NU = IO_obj.nu #viscosity nu (constant viscosity model)
alpha = IO_obj.alpha #Under relaxation
rho = IO_obj.rho #Density in Kg/m3
mu = NU*rho #constant for all nodes
#Boundary conditions (for u-velocity)
UA = IO_obj.UA
UB = IO_obj.UB
UC = IO_obj.UC
UD = IO_obj.UD

#Boundary conditions (for v-velocity)
VA = IO_obj.VA
VB = IO_obj.VB
VC = IO_obj.VC
VD = IO_obj.VD

#Initialize all variables
U = 0.0*grid
V = 0.0*grid
P = 0.0*grid
Pset = 0.0*grid
mdotw, mdote, mdotn, mdots  = 0.0*grid, 0.0*grid, 0.0*grid, 0.0*grid
k, omega, mut = 0.0*grid, 0.0*grid, 0.0*grid
resU, resV, resP, resB = 1.0 ,1.0 ,1.0 , 1.0
resplotU, resplotV, resplotP, resplotB = [1.0] ,[1.0],[1.0] ,[1.0]
outerIters = 0
eps = 1
#Apply BCs
U = solFunc_obj.applyBCs(U,UA,UB,UC,UD) #Apply U bcs
V = solFunc_obj.applyBCs(V,VA,VB,VC,VD) #Apply V bcs

########################################################################################################################################################################################################

caseType = "Laminar"

############################################################### Controls for iterations ###################################################################################################################

iters = 20 #number of gauss seidel sweeps
resInp = 1e-3 #residual for gauss seidel
residual = 1e-8 #residual for equations
interinp = 40 #number of outer iterations

################################################################################## SOLVER ##################################################################################################################

    #Step 1 : Solve for U,V using interpolated pressure using Discretize class obj
    # Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
    # --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
    # --> Implicit under-relaxation due to non-linearity in pde's
    # --> Returns face values of flux and velocities along with A and b matrices
while (outerIters < interinp):
    mdotwPrev, mdotePrev, mdotnPrev, mdotsPrev = mdotw, mdote, mdotn, mdots #mdot from previous iteration
    outerIters += 1
    print "Solving iteration %i"%(outerIters)
    aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp = disc_obj2.FOU_disc2( U,  V, mdotwPrev, mdotePrev, mdotnPrev, mdotsPrev , P)

    # #Step 1a : Solve for U using gauss seidel method
    # # --> Returns U* (newU) which will be corrected using rhie-chow interpolatio
    Ustar = gs_obj.gaussSeidel3u(U, aW, aE, aN, aS, aPmod,SUxmod, iters)
    #print Ustar

#     # #Step 1b : Solve for V using gauss seidel method
#     # # --> Returns V* (newV) which will be corrected using rhie-chow interpolation
#
    Vstar = gs_obj.gaussSeidel3u(V, aW, aE, aN, aS, aPmod,SUymod, iters)
#
#     # print ("APMOD AT %i iteration" %(outerIters))
#     # print Ustar
#     # print ("APMOD AT %i iteration" %(outerIters))
#     #Step 2 : RHIE-CHOW INTERPOLATION
#     # --> correct face velocities (EAST AND NORTH FACES ALONE)
#
    pcorre, pcorrn = rc_obj.rcInterp2(U, V, mdotw, mdote, mdotn, mdots, P, k, omega, mut, 'U', caseType)
#
#     # Step 2a : Correct face fluxes with rhie chow
    mstare, mstarn =  solFunc_obj.calcFaceMassFlux(Ustar,Vstar) # --> Get east and north face mass fluxes for Ustar and Vstar
    mdote,mdotn = solFunc_obj.correctFaceMassFlux(mstare,mstarn,pcorre,pcorrn) # --> Correct with rhie-chow to get corrected mdote and motn
    mdotw, mdots = solFunc_obj.getFaceMassFluxWS(mdote,mdotn) #--> Get corrected mdotw and mdots
#
#     #Step 3 : Solve P' equation using U velocities
#     # --> Calculates co-eff aP, aW, aW, aN, aS, at faces
#     # --> Calculates b matrix (Fw - Fe + Fs - Fn)
#     # --> Solve P' using gauss seidel
#
#     #Step 3a : Create B term (error in continuity)
    b = solFunc_obj.calcBterm(mdotw, mdote,mdotn ,mdots)
#
#     #Step 3b #Solve P' using co-eff for Pprime equation
#
    Pprime = gs_obj.gaussSeidel3u(Pset, aWpp, aEpp, aNpp, aSpp, aPpp, b, iters)
#
#     #Step 4 : Set pressure based on value at cell (2,2)
#     #Set P' level
    x = 1; y = 1; # Grid cell where P is fixed to zero
    Pset = solFunc_obj.setPress(Pprime,x,y)
#
#     #COPY Pset TO BOUNDARY
    Pset = solFunc_obj.setPsetbcs(Pset)
#
#     # #Step 5 : Pressure straddling
#
    mdoteNew, mdotnNew = corr_obj.massFluxcorr(Pset,mdote, mdotn, aEpp, aNpp)
#
    mdotwNew, mdotsNew = solFunc_obj.getFaceMassFluxWS(mdoteNew, mdotnNew)  # --> Get new mdotw and mdots
#
    uNew, vNew = corr_obj.velcorr(Ustar,Vstar,Pset, aPmod) # --> Correct u and v velcoities
#
#     # #Step 7 : Under relax P'
    Pnew = solFunc_obj.rlxP(P, Pset, alpha)
#
#     #Step 8 : Calculate residual
    # Step 10 : Calculate residual
    resU = res_obj.calcRes(U, aW, aE, aN, aS, SUxmod, aPmod)
    resV = res_obj.calcRes(V, aW, aE, aN, aS, SUymod, aPmod)
    resB = res_obj.calcResB(b)
    eps = max(resU, resV, resB)
    resplotU.append(resU)
    resplotV.append(resV)
    resplotB.append(resB)
#
#     #Replace U,V, mdotw, mdote, mdotn, mdots and P
    U = uNew
    V = vNew
    mdotw = mdotwNew; mdote = mdoteNew; mdotn = mdotnNew; mdots = mdotsNew
    P = Pnew
#
# # #Plot residuals
print U
xpt, ypt = plt_obj.genGrid(U)
plt_obj.plotdata2(U,V,ypt,resplotU, resplotV, resplotB,outerIters)
