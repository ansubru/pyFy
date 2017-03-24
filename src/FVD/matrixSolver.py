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
def DDChk(i,j,A):
    sum = np.zeros((i, 1), dtype=float)
    flag = []
    for m in range(i):
        for n in range(j):
            if m != n:
                sum[m] = sum[m] + abs(A[m][n])

    for m in range(i):
        for n in range(j):
            if m == n:
                if A[m][n] > sum[m]:
                    flag.append("DD")
                else:
                    flag.append("NXD")

    for m in range(flag.__len__()):
        if flag[m] in "NXD":
            print ("Critical error! Matrix is not diagonally dominant at row ",m)
            state = "fail"
        else:
            state = "pass"
    print state
    return state

def calcBterm (rho,ufw, ufe, ufs, ufn):
    b = rho*ufw - rho*ufe + rho*ufs - rho*ufn
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
#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
###############------------CREATE IO OBJECT--------------################

from IO import IO
IO_obj = IO("random")
grid = IO_obj.grid_nodes #generate a base grid node layout
NU = IO_obj.nu #viscosity nu (constant viscosity model)
rho = IO_obj.rho #Density in Kg/m3
mu = NU*rho #constant for all nodes
dx = IO_obj.x_dis # x-grid spacing
dy = IO_obj.y_dis # y-grid spacing

#Boundary conditions
UA = IO_obj.UA
UB = IO_obj.UB
UC = IO_obj.UC
UD = IO_obj.UD
VA = IO_obj.VA
VB = IO_obj.VB
VC = IO_obj.VC
VD = IO_obj.VD

#Initialize all variables

U = 0.0*grid
V = 0.0*grid
P = 0.0*grid

from Discretize import Discretize
disc_obj = Discretize()

from gaussSeidel2 import gaussSiedel2
gs_obj = gaussSiedel2()


#Step 1 : Solve for U using interpolated pressure using Discretize class obj
# Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
# --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# --> Implicit under-relaxation due to non-linearity in pde's
# --> Returns face values of flux and velocities along with A and b matrices

Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, A, Bx, By = disc_obj.FOU_disc(U,P, UA, UB, UC, UD)

## A matrix checks!
#Check if matrix A is diagonally dominant
i = np.size(U,0) #get indices for rows
j = np.size(U,1) #get indices for columns
testVal = DDChk(i, j, A)

#Step 2 : Solve for U using gauss seidel method
# --> Returns U* (newU) which will be corrected using rhie-chow interpolation
iters = 5
Ustar = gs_obj.gaussSeidel2u(U, aW, aE, aN, aS, aPmod,SUxmod, iters)

#Step 3 : Discretize newU using FOU for rhie-chow interpolation using Discretize class obj
# --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# --> Returns face values of flux and velocities along with A and b matrices

from rhieChow import rhieChow
rc_obj = rhieChow()
ufwrc, uferc, ufnrc, ufsrc, pcorrw, pcorre, pcorrn, pcorrs = rc_obj.rcInterp(Ustar,P)

#Step 4 : Solve P' equation
# --> Calculates co-eff aP, aW, aW, aN, aS, at faces
# --> Calculates b matrix (Fw - Fe + Fs - Fn)
# --> Solve P' using gauss seidel

#Create the b term (continuity)

b = 0.0*Ustar
for m in range(i): #loop through rows
    for n in range(j): #loop through columns
        b[m][n] = calcBterm(rho,ufwrc[m][n],uferc[m][n],ufsrc[m][n],ufnrc[m][n])

#Get face interpolated final coeffs for pressure correction equation

from pressCorr import pressCorr
pressCorr_obj = pressCorr()
aWp, aEp, aSp, aNp, aPp= pressCorr_obj.pcorr(Ustar,P)

#Fix boundaries for naumann bc's (zero pressure at cells close to the bcs)

aWp = fixCoeffsP(aWp,side="W")
aEp = fixCoeffsP(aEp,side="E")
aSp = fixCoeffsP(aSp,side="S")

#Solve P'
Pprime = gs_obj.gaussSeidel2u(P, aWp, aEp, aNp, aSp, aPp, b , iters)

#Step 5 : Set pressure based on value at cell (2,2)

#Set P' level
Pset = setPress(Pprime)
print Pset

#Step 6 : Fix face velocities - Pressure straddling
#Unew =



