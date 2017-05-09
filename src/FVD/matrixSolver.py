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
################################################################################MODULE FOR DISCRETIZATION################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

########################################################################FUNCTION DEF######################################################################################################################

def calcBterm (rho,ufw, ufe, ufs, ufn,dx,dy):
    b = rho*ufw*dy - rho*ufe*dy + rho*ufs*dx - rho*ufn*dx
    return b

def fixCoeffsP (coeff,side):
    i = np.size(coeff, 0)
    j = np.size(coeff, 1)

    if side in 'W':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  != (i-1) and n == 0 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == 0 :
                    coeff[m][n] = 0.0

    if side in 'E':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  != (i-1) and n == j-1 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == j-1 :
                    coeff[m][n] = 0.0

    if side in 'S':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  == (i-1) and n != j-1 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == j-1 :
                    coeff[m][n] = 0.0
    return coeff

def setPress(P):
    i = np.size(P, 0)
    j = np.size(P, 1)
    for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  == 1 and n == 1 :
                    P[m][n] = 0.0
                else:
                    P[m][n] = P[m][n] - P[1][1]
    return P

def rlxP(P, PP, coeff):
    i = np.size(P, 0)
    j = np.size(P, 1)
    pPnew = 0.0*P
    for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                pPnew[m][n] = P[m][n] + coeff*PP[m][n]
    return pPnew
#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
###############------------CREATE IO OBJECT--------------################

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

from Discretize2 import Discretize2
disc_obj2 = Discretize2()

from gaussSeidel2 import gaussSiedel2
gs_obj = gaussSiedel2()

from gaussSiedel import gaussSiedel
gs_obj1 = gaussSiedel()

###################################Controls for iterations####################################################################

iters = 5 #number of gauss seidel sweeps
alpha = 0.5 #Under relaxation

###################################Controls for iterations####################################################################
#Step 1a : Solve for U using interpolated pressure using Discretize class obj
# Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
# --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# --> Implicit under-relaxation due to non-linearity in pde's
# --> Returns face values of flux and velocities along with A and b matrices

Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, A, Bx, By = disc_obj2.FOU_disc2(U,P, UA, UB, UC, UD)
print Fw[1, :]
print Fe[1, :]
print Fn[1, :]
print Fs[1, :]


# i = np.size(U,0) #get indices for rows
# j = np.size(U,1) #get indices for columns
#
#
# #Step 2 : Solve for U using gauss seidel method
# # --> Returns U* (newU) which will be corrected using rhie-chow interpolation
# print "Solving for Ustar using gauss seidel"
#
# Ustar = gs_obj.gaussSeidel2u(U, aW, aE, aN, aS, aPmod,SUxmod, iters)
#
# print "Done with %d sweeps" %(iters,)
# #Step 3 : Discretize newU using FOU for rhie-chow interpolation using Discretize class obj
# # --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# # --> Returns face values of flux and velocities along with A and b matrices
#
# from rhieChow import rhieChow
# rc_obj = rhieChow()
# ufwrc, uferc, ufnrc, ufsrc, pcorrw, pcorre, pcorrn, pcorrs = rc_obj.rcInterp(Ustar,P)
#
# # Get velocities from face to grid nodes
# from interpToGrid import interpToGrid
# intrgrd_obj = interpToGrid()
# Uprc = intrgrd_obj.gridInterpU(ufwrc, uferc)
#
# #Step 1a : Solve for V using interpolated pressure using Discretize class obj
# # Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
# # --> Calculates co-ef aP, aW, aW, aN, aS, and Sources
# # --> Implicit under-relaxation due to non-linearity in pde's
# # --> Returns face values of flux and velocities along with A and b matrices
#
# Fev, Fwv, Fnv, Fsv, ufev, ufwv, ufnv, ufsv, aWv, aEv, aNv, aSv,aWpv, aEpv, aNpv, aSpv, aPv, aPmodv, SUxmodv, SUymodv, Av, Bxv, Byv = disc_obj.FOU_disc(V,P, VA, VB, VC, VD)
#
# i = np.size(V,0) #get indices for rows
# j = np.size(V,1) #get indices for columns
#
#
# #Step 2 : Solve for V using gauss seidel method
# # --> Returns V* (newV) which will be corrected using rhie-chow interpolation
#
# print "Solving for Vstar using gauss seidel"
#
# Vstar = gs_obj.gaussSeidel2u(V, aWv, aEv, aNv, aSv, aPmodv,SUymodv, iters)
#
# print "Done with %d sweeps" %(iters,)
#
# #Step 3 : Discretize newV using FOU for rhie-chow interpolation using Discretize class obj
# # --> Calculates co-eff aPv, aWv, aWv, aNv, aSv, and Sources
# # --> Returns face values of flux and velocities along with A and b matrices
#
# from rhieChow import rhieChow
# rc_obj = rhieChow()
# vfwrc, vferc, vfnrc, vfsrc, pcorrwv, pcorrev, pcorrnv, pcorrsv = rc_obj.rcInterp(Vstar,P)
#
# # Get velocities from face to grid nodes
# from interpToGrid import interpToGrid
# intrgrd_obj = interpToGrid()
# Vprc = intrgrd_obj.gridInterpU(vfwrc, vferc)
#
#
# #Step 4 : Solve P' equation using U velocities
# # --> Calculates co-eff aP, aW, aW, aN, aS, at faces
# # --> Calculates b matrix (Fw - Fe + Fs - Fn)
# # --> Solve P' using gauss seidel
#
# #Create the b term (continuity)
#
# b = 0.0*Ustar
# for m in range(i): #loop through rows
#     for n in range(j): #loop through columns
#         b[m][n] = calcBterm(rho,ufwrc[m][n],uferc[m][n],vfsrc[m][n],vfnrc[m][n], dx, dy)
#
# #Get face interpolated final coeffs for pressure correction equation
#
# from pressCorr import pressCorr
# pressCorr_obj = pressCorr()
# aWp, aEp, aSp, aNp, aPp= pressCorr_obj.pcorr(Ustar,P)
#
# #Fix boundaries for naumann bc's (zero coeffs at cells close to the bcs)
#
# aWp = fixCoeffsP(aWp,side="W")
# aEp = fixCoeffsP(aEp,side="E")
# aSp = fixCoeffsP(aSp,side="S")
#
# #Solve P'
# print "Solving P prime equation using gauss seidel"
# Pprime = gs_obj.gaussSeidel2u(P, aWp, aEp, aNp, aSp, aPp, b , iters)
# print "Done with %d sweeps" %(iters,)
#
# print Ustar
#
# #Step 5 : Set pressure based on value at cell (2,2)
#
# #Set P' level
# Pset = setPress(Pprime)
#
# #Step 6 : Fix face velocities - Pressure straddling
# from veloCorr import veloCorr
# velcorr_obj = veloCorr()
# uNew = velcorr_obj.velcorrU(Uprc,Pset, aPmod)
# vNew = velcorr_obj.velcorrV(Vprc,Pset,aPmodv)
#
# #Step 7 : Under relax P'
# Pnew = rlxP(P, Pset, alpha )
#
# #Replace U,V and P
# U = uNew
# V = vNew
# P = Pnew
#
# # Plot residuals
#
# from plotResults import plotResults
# plt_obj = plotResults()
# xpt, ypt = plt_obj.genGrid(U)

#plt_obj.plotContours(xpt,ypt*0.1,U,V,P,Pprime,1)

