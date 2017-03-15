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
i = np.size(U,0) #get indices for rows
j = np.size(U,1) #get indices for columns

from Discretize import Discretize
disc_obj = Discretize()

# Discretize U velocity using FOU
Px, dPx, dPy, Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aP, aW, aE, aN, aS, SUx, SPx, SUy, SPy = disc_obj.FOU_disc(U, UA, UB, UC, UD)

#Generate A, b and x matrices for Gauss Seidel method
AU = []
amat = [0 for y in range(i)]
for x in range(i):
    for m in range(i): #loop through rows
        for n in range(j): #loop through columns
            print m,n
            amat[m] = [ aP[m][n], (-aW[m][n]), (-aE[m][n]), (-aN[m][n]), (-aS[m][n]) ]
    AU.append(amat[m])
    print AU

AU = np.asarray(AU)
print (np.size(AU,0))

# Discretize V velocity using FOU
Px, dPx, dPy, Fe, Fw, Fn, Fs, vfe, vfw, vfn, vfs, aP, aW, aE, aN, aS, SVx, SPx, SVy, SPy = disc_obj.FOU_disc(V, VA, VB, VC, VD)



'''___MAIN___'''

A = np.array([[4.0, -2.0, 1.0], [1.0, -3.0, 2.0], [-1.0, 2.0, 6.0]])
b = [1.0, 2.0, 3.0]
x = [1, 1, 1]

n = 20

#print gauss(A, b, x, n)
#print solve(A, b)




