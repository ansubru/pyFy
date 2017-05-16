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
################################################################################ MAIN SOLVER ################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

################################################################ INPUT TO SOLVER ###########################################################################################################################


from IO import IO
IO_obj = IO("random")
grid = IO_obj.grid_nodes #generate a base grid node layout
NU = IO_obj.nu #viscosity nu (constant viscosity model)
rho = IO_obj.rho #Density in Kg/m3
mu = NU*rho #constant for all nodes
delXW = IO_obj.delXW #distance to Western neighbour
delXE = IO_obj.delXE #distance to Eastern neighbour
delYN = IO_obj.delYN #distance to Northern neighbour
delYS = IO_obj.delYS #distance to Northern neighbour
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
mdotw, mdote, mdotn, mdots,  = 0.0*grid, 0.0*grid, 0.0*grid, 0.0*grid,

# Apply bcs to U an V
i = np.size(U,0) #get indices for rows
j = np.size(U,1) #get indices for columns

for m in range(0, j):
    for n in range(0, i):
        if (m >= 0 and n == 0): # A Boundary
            U[m][n] = UA  # pad bc velocities
            V[m][n] = VA  # pad bc velocities
        if (m == 0 and n != j-1 and n != 0): # B Boundary
            U[m][n] = UB  # pad bc velocities
            V[m][n] = VB  # pad bc velocities
        if (m >= 0 and n == j-1): # C Boundary
            U[m][n] = UC  # pad bc velocities
            V[m][n] = VC  # pad bc velocities
        if (m == i-1  and n != j-1 and n != 0): # D Boundary
            U[m][n] = UD  # pad bc velocities
            V[m][n] = VD  # pad bc velocities

from Discretize2 import Discretize2
disc_obj2 = Discretize2()

from gaussSeidel2 import gaussSiedel2
gs_obj = gaussSiedel2()

from solFunc import solFunc
solFunc_obj = solFunc()


def gaussSeidel3u(U, aW, aE, aN, aS, aP, SUx, iterinp):
    """Solves gauss seidel for a fixed number of iterations."""

    u = 1.0*U

    # Initialize all relevant variables:u, v, p etc.
    i = np.size(u, 0)
    j = np.size(u, 1)

    iter = 0  # Mock up error parameter only to run the While
    while iter < iterinp:
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes:

                    uN = u[m - 1][n]
                    uS = u[m + 1][n]
                    uE = u[m][n + 1]
                    uW = u[m][n - 1]

                    u[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                               uN * aN[m][n] + SUx[m][n]) / aP[m][n]

        iter += 1
    return u

############################################################### Controls for iterations ###################################################################################################################

iters = 5 #number of gauss seidel sweeps
alpha = 0.1 #Under relaxation
outerIters = 0
################################################################################### SOLVER ##################################################################################################################

    #Step 1 : Solve for U,V using interpolated pressure using Discretize class obj
    # Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
    # --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
    # --> Implicit under-relaxation due to non-linearity in pde's
    # --> Returns face values of flux and velocities along with A and b matrices
while (outerIters < 1):
    outerIters += 1
    print "Solving iteration %i"%(outerIters)
    aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp = disc_obj2.FOU_disc2( U,  V, mdotw, mdote, mdotn, mdots , P)

    # #Step 1a : Solve for U using gauss seidel method
    # # --> Returns U* (newU) which will be corrected using rhie-chow interpolatio
    Ustar = gaussSeidel3u(U, aW, aE, aN, aS, aPmod,SUxmod, iters)
    print Ustar
    # #Step 1b : Solve for V using gauss seidel method
    # # --> Returns V* (newV) which will be corrected using rhie-chow interpolation

    Vstar = gaussSeidel3u(V, aW, aE, aN, aS, aPmod,SUymod, iters)

    # print ("APMOD AT %i iteration" %(outerIters))
    # print Ustar
    # print ("APMOD AT %i iteration" %(outerIters))
    #Step 2 : RHIE-CHOW INTERPOLATION
    # --> correct face velocities (EAST AND NORTH FACES ALONE)

    from rhieChow2 import rhieChow2
    rc_obj = rhieChow2()
    pcorre, pcorrn = rc_obj.rcInterp2(U,  V, mdotw, mdote, mdotn, mdots , P)

    # Step 2a : Correct face fluxes with rhie chow
    mstare, mstarn =  solFunc_obj.calcFaceMassFlux(Ustar,Vstar) # --> Get east and north face mass fluxes for Ustar and Vstar
    mdote,mdotn = solFunc_obj.correctFaceMassFlux(mstare,mstarn,pcorre,pcorrn) # --> Correct with rhie-chow to get corrected mdote and motn
    mdotw, mdots = solFunc_obj.getFaceMassFluxWS(mdote,mdotn) #--> Get corrected mdotw and mdots


    #Step 3 : Solve P' equation using U velocities
    # --> Calculates co-eff aP, aW, aW, aN, aS, at faces
    # --> Calculates b matrix (Fw - Fe + Fs - Fn)
    # --> Solve P' using gauss seidel

    #Step 3a : Create B term (error in continuity)
    b = solFunc_obj.calcBterm(mdotw, mdote,mdotn ,mdots)


    #Step 3b #Solve P' using co-eff for Pprime equation

    Pprime = gaussSeidel3u(Pset, aWpp, aEpp, aNpp, aSpp, aPpp, b, iters)

    #Step 4 : Set pressure based on value at cell (2,2)
    #Set P' level
    x = 1; y = 1; # Grid cell where P is fixed to zero
    Pset = solFunc_obj.setPress(Pprime,x,y)

    #COPY Pset TO BOUNDARY
    Pset = solFunc_obj.setPsetbcs(Pset)



    # #Step 5 : Pressure straddling
    from Corr2 import Corr2
    corr_obj = Corr2()
    mdoteNew, mdotnNew = corr_obj.massFluxcorr(Pset,mdote, mdotn, aEpp, aNpp)

    mdotwNew, mdotsNew = solFunc_obj.getFaceMassFluxWS(mdoteNew, mdotnNew)  # --> Get new mdotw and mdots

    uNew, vNew = corr_obj.velcorr(Ustar,Vstar,Pset, aPmod) # --> Correct u and v velcoities

    # #Step 7 : Under relax P'
    Pnew = solFunc_obj.rlxP(P, Pset, alpha)

    #Replace U,V, mdotw, mdote, mdotn, mdots and P
    U = uNew
    V = vNew
    mdotw = mdotwNew; mdote = mdoteNew; mdotn = mdotnNew; mdots = mdotsNew
    P = Pnew

# # #Plot residuals
# from plotResults import plotResults
# plt_obj = plotResults()
# xpt, ypt = plt_obj.genGrid(U)
# plt_obj.plotContours(xpt,ypt,U,V,P,Pprime,outerIters)
