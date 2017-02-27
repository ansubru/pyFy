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
#############################################################################################################################################################################################################
###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###

###############------------CREATE IO OBJECT--------------################

from IO import IO
IO_obj = IO("random")
xNodes = IO_obj.x #Number of x xnodes
yNodes = IO_obj.y #Number of ynodes
grid = IO_obj.grid_zeros
grid_xpos = IO_obj.grid_x #x co-ordinates
grid_ypos = IO_obj.grid_y #y co-ordinates
aE = aW = IO_obj.delta_y #area of east and west faces
aN = aS = IO_obj.delta_x #area of north and south faces
NU = IO_obj.nu #viscosity nu (constant viscosity model)
#Create array with x,y co-ordinates from generated grid

###############------------------------------------------################

x = 1
y = 2
z = 3

###############-------CREATE INTERPOLATE OBJECT----------################

from Interpolate import Interpolate
interp = Interpolate(x, y, z)
print interp.lin_interp()

###############------------------------------------------################
