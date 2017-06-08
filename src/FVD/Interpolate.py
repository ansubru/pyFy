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

    def __init__(self):
        """Return object"""


    def lin_interp(self,x,y):
        """Function performs a linear interpolation"""
        lininterp = (x + y)/2
        return lininterp

    def CD_interp(self,x,y,z):
        """Function that performs a Central differencing"""
        CDinterp = (x-y)/z
        return CDinterp

    def weighted_interp(self,NP,NB,f):
        """Function that performs interpolation based on provided interpolation weights"""
        WTinterp = (f*NB) + ((1.0-f)*NP)
        return WTinterp

##############################################################################################################################################################################################################


