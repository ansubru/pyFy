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

from Discretize import Discretize
disc_obj = Discretize()

from gaussSiedel import gaussSiedel
gs_obj = gaussSiedel()
# Discretize U velocity using FOU
Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, A, Bx, By = disc_obj.FOU_disc(U,P, UA, UB, UC, UD)
#Check if matrix A is diagonally dominant
i = np.size(A,0) #get indices for rows
j = np.size(A,1) #get indices for columns
sum = np.array([])
for m in range(i):
    for n in range(j):
        if m != n:
            sum[i] = A[]

#Ugs = gs_obj.gauss_seidel(AU, x0=None, eps=1e-5, max_iteration=100)



# Discretize V velocity using FOU
#Px, dPx, dPy, Fe, Fw, Fn, Fs, vfe, vfw, vfn, vfs, aP, aW, aE, aN, aS, SVx, SPx, SVy, SPy = disc_obj.FOU_disc(V, VA, VB, VC, VD)
pp = pprint.PrettyPrinter(indent=1)
print 'A:'
pp.pprint(A)
print 'B:'
pp.pprint(Bx)







