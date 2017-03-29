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
class interpToGrid(object):

    """A class for interpolating u and v face velo to grid

    Attributes:
       U --> Feed in U face velocities (ufe and ufw) (x-direction)
       v --> Feed in V velocity(vfn and vfs) (y-direction)
    """

    def __init__(self):
        """Return object"""

    def gridInterpU(self, ufw, ufe):
        """Class to interpolate to grid"""

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        uw = 1.0*ufw
        ue = 1.0*ufe

        i = np.size(uw, 0)
        j = np.size(uw, 1)

        uPfinal = 0.0*uw


        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                uPfinal[m][n] = Interp_obj.lin_interp(uw[m][n],ue[m][n])

        return uPfinal

    def gridInterpV(self, vfs, vfn):
        """Class to interpolate to grid"""

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        vs = 1.0*vfs
        vn = 1.0*vfn

        i = np.size(vs, 0)
        j = np.size(vs, 1)

        vPfinal = 0.0*vs

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                vPfinal[m][n] = Interp_obj.lin_interp(vs[m][n],vn[m][n])

        return vPfinal
