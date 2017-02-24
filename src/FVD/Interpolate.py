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
################################################################################MODULE CLASS FOR INTERPOLATION SCHEMES################################################################################################

import os
import re

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class Interpolate(object):

    """A class for interpolating relevant nodes

    Attributes:
       xnode_val: Position of xNode
       ynode_val: Position of yNode
           z_val: Distance between nodal points
    """

    def __init__(self, x, y, z):
        """Return object"""
        self.x = x
        self.y = y
        self.z = z

    def lin_interp(self):
        """Function performs a linear interpolation"""
        lininterp = (self.x+self.y)/2
        return lininterp

    def CD_interp(self):
        """Function that performs a Central differencing"""
        CDinterp = (self.x-self.y)/self.z
        return CDinterp

##############################################################################################################################################################################################################

