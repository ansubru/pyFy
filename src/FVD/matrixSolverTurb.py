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
################################################################################ TURBULENT FLOW SOLVER ################################################################################################

import time
import math
import pylab

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
ypt = IO_obj.ypt

#Boundary conditions (for v-velocity)
VA = IO_obj.VA
VB = IO_obj.VB
VC = IO_obj.VC
VD = IO_obj.VD

#Initialize all variables
U = 0.0*grid
V = 0.0*grid
ubar = 0.0*grid
vbar = 0.0*grid
P = 0.0*grid
Pset = 0.0*grid
mut = 0.0*grid
muts = 0.0*grid
mdotw, mdote, mdotn, mdots  = 0.0*grid, 0.0*grid, 0.0*grid, 0.0*grid
resU, resV, resP, resB = 1.0 ,1.0 ,1.0 , 1.0
resplotU, resplotV, resplotk, resplotomg, resplotB, resplotPP = [1.0] ,[1.0],[1.0] ,[1.0], [1.0], [1.0]
outerIters = 0
eps = 1
#Apply BCs
U = solFunc_obj.applyBCs(U,UA,UB,UC,UD) #Apply U bcs
V = solFunc_obj.applyBCs(V,VA,VB,VC,VD) #Apply V bcs
kinit = 0.01 #initial guess
omegainit = 100 #initial guess
k = solFunc_obj.initializeK(kinit)
omegain = solFunc_obj.initializeOmega(omegainit)

# COPY omegaw TO BOUNDARY
omega = solFunc_obj.setPsetbcs(omegain)
#Initialize turbulent viscosity
#mut, muts = 0.0*grid, 0.0*grid
#mut, muts = solFunc_obj.calcMuts(k, omega)
#Timing simulations
start = time.clock()
#Plot output display control
n = 1 #plots displayed at multiples of n iterations
########################################################################################################################################################################################################

caseType = "Turbulent"

############################################################### Controls for iterations ###################################################################################################################

iters = 20 #number of gauss seidel sweeps
iters2 = 2
resInp = 1e-6 #residual for gauss seidel (internal residual for GS sweeps)
residual = 1e-5 #residual for all equations (convergence criteria)
interinp = 100 #number of outer iterations
################################################################################### SOLVER ##################################################################################################################

    #Step 1 : Solve for U,V using interpolated pressure using Discretize class obj
    # Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
    # --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
    # --> Implicit under-relaxation due to non-linearity in pde's
    # --> Returns face values of flux and velocities along with A and b matrices
#while (eps > residual):
while (outerIters < interinp):
    mdotwPrev, mdotePrev, mdotnPrev, mdotsPrev = mdotw, mdote, mdotn, mdots #mdot from previous iteration
    mutPrev = mut
    outerIters += 1

    print "Solving iteration %i"%(outerIters)
    aW, aE, aN, aS, aWp, aEp, aNp, aSp, aP, aPmod, SUmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp = disc_obj2.FOU_discTurb2( U, V, k, omega, mdotwPrev, mdotePrev, mdotnPrev, mdotsPrev , mutPrev, P, 'U')

    # #Step 1a : Solve for U using gauss seidel method
    # # --> Returns U* (newU) which will be corrected using rhie-chow interpolation
    Ustar1 = gs_obj.gaussSeidel6u(U, aW, aE, aN, aS, aPmod,SUxmod, iters2, "U")
    Ustar = gs_obj.gaussSeidel3u(U, aW, aE, aN, aS, aPmod, SUxmod, iters, "U")

    # #Step 1b : Solve for V using gauss seidel method
    # # --> Returns V* (newV) which will be corrected using rhie-chow interpolation
    Vstar1 = gs_obj.gaussSeidel6u(V, aW, aE, aN, aS, aPmod,SUymod, iters2, "V")
    Vstar = gs_obj.gaussSeidel3u(V, aW, aE, aN, aS, aPmod, SUymod, iters, "V")

    #Step 2 : RHIE-CHOW INTERPOLATION
    # --> correct face velocities (EAST AND NORTH FACES ALONE)

    pcorre, pcorrn = rc_obj.rcInterp2(U,  V, mdotw, mdote, mdotn, mdots ,P, k, omega, mut, 'U', caseType)

    ## Step 2a : Correct face fluxes with rhie chow
    mstare, mstarn =  solFunc_obj.calcFaceMassFluxIW(Ustar,Vstar) # --> Get east and north face mass fluxes for Ustar and Vstar (INTERPOLATION WEIGHTED)
    mdote,mdotn = solFunc_obj.correctFaceMassFlux(mstare,mstarn,pcorre,pcorrn) # --> Correct with rhie-chow to get corrected mdote and motn
    mdotw, mdots = solFunc_obj.getFaceMassFluxWS(mdote,mdotn) #--> Get corrected mdotw and mdots

    #Step 3 : Solve P' equation using U velocities
    # --> Calculates co-eff aP, aW, aW, aN, aS, at faces
    # --> Calculates b matrix (Fw - Fe + Fs - Fn)
    # --> Solve P' using gauss seidel

    #Step 3a : Create B term (error in continuity)
    b = solFunc_obj.calcBterm(mdotw, mdote,mdotn ,mdots)

    #Step 3b #Solve P' using co-eff for Pprime equation

    resInpPP = 1e-10
    itersSpcl = 3000
    Pprime1 = gs_obj.gaussSeidel6u(Pset, aWpp, aEpp, aNpp, aSpp, aPpp, b, iters2, "P")
    Pprime = gs_obj.gaussSeidel3u(Pset, aWpp, aEpp, aNpp, aSpp, aPpp, b, itersSpcl, "Pprime")

    #Step 4 : Set pressure based on value at cell (2,2)
    #Set P' level
    x = 1; y = 1; # Grid cell where P is fixed to zero
    Pset = solFunc_obj.setPress(Pprime,x,y)

    #COPY Pset TO BOUNDARY
    Pset = solFunc_obj.setPsetbcs(Pset)

    # #Step 5 : Pressure straddling

    mdoteNew, mdotnNew = corr_obj.massFluxcorr(Pset,mdote, mdotn, aEpp, aNpp)
    mdotwNew, mdotsNew = solFunc_obj.getFaceMassFluxWS(mdoteNew, mdotnNew)  # --> Get new mdotw and mdots
    uNew, vNew = corr_obj.velcorrTurb(Ustar,Vstar,Pset, aPmod) # --> Correct u and v velcoities


    # #Step 7 : Under relax P'
    Pnew = solFunc_obj.rlxP(P, Pset, alpha)

    ## Step 8 : Solve k equation
    ## Calculate turbulent viscosity
    aWk, aEk, aNk, aSk, aWpk, aEpk, aNpk, aSpk, aPk, aPmodk, SUmodk, SUxmodk, SUymodk, aWppk, aEppk, aNppk, aSppk, aPppk \
                                                    = disc_obj2.FOU_discTurb2(uNew, vNew, k, omega, mdotwNew, mdoteNew, mdotnNew, mdotsNew, muts, P, "k")

    # #Step 8a : Solve for k using gauss seidel method
    # # --> Returns newk which will be used to solve for omega
    kNew1 = gs_obj.gaussSeidel6u(k, aWk, aEk, aNk, aSk, aPmodk,SUmodk, iters2, "k")
    kNew = gs_obj.gaussSeidel3u(k, aWk, aEk, aNk, aSk, aPmodk, SUmodk, iters, "k")

    ## Step 9 : Solve omega equation
    ## update turbulent viscosity with kNew
    mutkNew, mutskNew = solFunc_obj.calcMuts(kNew, omega)
    aWomg, aEomg, aNomg, aSomg, aWpo, aEpo, aNpo, aSpo, aPomg, aPmodomg, SUmodomg, SUxmodo, SUymodo, aWppo, aEppo, aNppo, aSppo, aPppo \
        = disc_obj2.FOU_discTurb2(uNew, vNew, kNew, omega, mdotwNew, mdoteNew, mdotnNew, mdotsNew, mutskNew, P, "omega")

    #omegaNew = gs_obj.gaussSeidel4u(omega, aWomg, aEomg, aNomg, aSomg, aPmodomg,SUmodomg, resInp, "omega")
    omegaNew1 = gs_obj.gaussSeidel6u(omega, aWomg, aEomg, aNomg, aSomg, aPmodomg, SUmodomg, iters2, "omega")
    omegaNew = gs_obj.gaussSeidel3u(omega, aWomg, aEomg, aNomg, aSomg, aPmodomg, SUmodomg, iters, "omega")
#
#   #Step 10 : Calculate residual
    resU = res_obj.calcRes(U,  aW, aE, aN, aS, SUxmod, aPmod, b, "U")
    resV = res_obj.calcRes(V,  aW, aE, aN, aS, SUymod, aPmod, b, "V")
    resK = res_obj.calcRes(k,  aWk, aEk, aNk, aSk, SUmodk, aPmodk, b, "K")
    resPP = res_obj.calcRes(Pset, aWpp, aEpp, aNpp, aSpp,  b, aPpp, b, "PP")
    resomg = res_obj.calcRes(omega, aWomg, aEomg, aNomg, aSomg, SUmodomg, aPmodomg, b, "omega")
    resB = res_obj.calcRes(omega, aWomg, aEomg, aNomg, aSomg, SUmodomg, aPmodomg, b, "B")

    eps = max(resU, resV, resK, resomg, resB, resPP) # max residual for the run

    if eps > 1e5:
        print ("*#!%&*!%/*%/*!¤%&    F***K!! WE DIVERGED :( ...     *#!%&*!%/*%/*!¤%&")
        break

    resplotU.append(resU)
    resplotV.append(resV)
    resplotk.append(resK)
    resplotPP.append(resPP)
    resplotomg.append(resomg)
    resplotB.append(resB)
    elapsed = (time.clock() - start)
    print ("Ures: %e, Vres: %e, PPres: %e, kres: %e, omegares: %e, Continuity: %e, maxRes(eps): %e"%(resU, resV, resPP, resK, resomg, resB, eps))
    print ("%d seconds elapsed (Execution time)"%(elapsed))
#
    #Replace U,V, mdotw, mdote, mdotn, mdots,P, mut, k, omega
    U = uNew
    V = vNew
    mdotw = mdotwNew; mdote = mdoteNew; mdotn = mdotnNew; mdots = mdotsNew
    P = Pnew
    k = kNew
    omega = omegaNew
    mut, muts = solFunc_obj.calcMuts(k, omega)

# # #Plots at simulation end
print ("***!!***!!***!!***!!***!!***    HOUSTON!! WE HAVE CONVERGED :D ...     ***!!***!!***!!***!!***!!***")
plt_obj.plotdataTurb(U, V, k, omega, ypt, mut, resplotU, resplotV, resplotk, resplotomg, resplotB, resplotPP, outerIters)


