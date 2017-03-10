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

# Discretize U velocity using FOU
Px, dPx, dPy, Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aP, aW, aE, aN, aS, SUx, SPx, SUy, SPy = disc_obj.FOU_disc(U, UA, UB, UC, UD)
print 'aW:'
pprint.pprint(aW)
print 'aE:'
pprint.pprint(aE)
print 'aS:'
pprint.pprint(aS)
print 'aN:'
pprint.pprint(aN)
print 'aP:'
pprint.pprint(aP)






