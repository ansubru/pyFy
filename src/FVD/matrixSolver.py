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

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
###############------------CREATE IO OBJECT--------------################

from IO import IO
IO_obj = IO("random")
grid = IO_obj.grid_nodes #generate a base grid node layout

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

from gaussSiedel import gaussSiedel
gs_obj = gaussSiedel()


#Step 1 : Solve for U using interpolated pressure using Discretize class obj
# Discretize U velocity using FOU --> P is checkerboarded (no rhie-chow interpolation)
# --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# --> Implicit under-relaxation due to non-linearity in pde's
# --> Returns face values of flux and velocities along with A and b matrices

Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS, aP, aPmod, A, Bx, By = disc_obj.FOU_disc(U,P, UA, UB, UC, UD)

## A matrix checks!
#Check if matrix A is diagonally dominant
i = np.size(A,0) #get indices for rows
j = np.size(A,1) #get indices for columns
testVal = DDChk(i, j, A)

#Step 2 : Solve for U using gauss seidel method
# --> Returns U* (newU) which will be corrected using rhie-chow interpolation
if testVal in "pass":
    reltol = 0.001
    toltype = "r"
    solve = gs_obj.gauss_seidel(np.array(A), np.array(Bx), toltype, reltol)
np.array(solve)
newU = solve.reshape(np.size(grid,0),np.size(grid,1))

#Step 3 : Discretize newU using FOU for rhie-chow interpolation using Discretize class obj
# --> Calculates co-eff aP, aW, aW, aN, aS, and Sources
# --> Returns face values of flux and velocities along with A and b matrices

from rhieChow import rhieChow
rc_obj = rhieChow()

ufwrc, uferc, ufnrc, ufsrc, pcorrw, pcorre, pcorrn, pcorrs = rc_obj.rcInterp(newU,P)

from interpToFace import interpToFace

interpFc_obj = interpToFace()

aWw, aEe, aSs, aNn, aPe = interpFc_obj.faceInterp(newU,P)

pp = pprint.PrettyPrinter(indent=6)
pp.pprint(aPe)





