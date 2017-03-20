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
# Discretize U velocity using FOU
Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, A, Bx, By = disc_obj.FOU_disc(U,P, UA, UB, UC, UD)
print Bx
## A matrix checks!

#Check if matrix A is diagonally dominant
i = np.size(A,0) #get indices for rows
j = np.size(A,1) #get indices for columns
testA = np.ones(shape=(i, j))
testVal = DDChk(i, j, A)

if testVal in "pass":
    reltol = 0.001
    toltype = "r"
    solve = gs_obj.gauss_seidel(np.array(A), np.array(Bx), toltype, reltol)
np.array(solve)
newU = solve.reshape(np.size(grid,0),np.size(grid,1))
print newU








