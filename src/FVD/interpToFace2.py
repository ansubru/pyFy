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
class interpToFace2(object):

    """A class to interpolate aP to faces

    Attributes:

    """

    def __init__(self):
        """Return object"""

    def faceInterp2(self, matU, matV, mdotw, mdote, mdotn, mdots ,matP):
        """Function to interpolate to faces"""

        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        grid = IO_obj.grid_nodes #generate a base grid node layout
        u = np.array(1.0 * matU)
        v = np.array(1.0 * matU)
        p = np.array(1.0 * matP)



        from Discretize2 import Discretize2
        disc_obj = Discretize2()

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        # Get co-efficients

        Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS,anotW, anotE, anotN, anotS, aP, aPmod, SUxmod, SUymod, A, Bx, By = disc_obj.FOU_disc2( u, UA, UB, UC, UD, v, VA, VB, VC, VD , p)

        i = np.size(u, 0)
        j = np.size(u, 1)
        aPe = 0.0*u  # interpolated co-eff to faces
        aPw, aPs, aPn = 0.0*u, 0.0*u, 0.0*u

        c = 0
        d = 0
        q = 0
        z = 0
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m > 1 and m <= i - 3 and n > 1 and n <= j - 3):  # Internal nodes (except boundary + first grid nodes)

                    # West faces
                    aPw[m][n] = Interp_obj.lin_interp(aPmod[m][n-1], aPmod[m][n])

                    # East faces
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n+1])

                    # South faces
                    aPs[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m+1][n])

                    # North faces
                    aPn[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m-1][n])

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    c+=1
                    aPe[m][n] = 0.0  # Neumann bc at all boundaries
                    aPw[m][n] = Interp_obj.lin_interp(aPmod[m][n - 1], aPmod[m][n])

                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    d += 1
                    aPw[m][n] = 0.0 # Naumann bc at all boundaries
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    q += 1
                    aPn[m][n] = 0.0 # Naumann bc at all boundaries
                    aPs[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m + 1][n])

                if (m == i - 2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    z += 1
                    aPs[m][n] = 0.0 # Naumann bc at all boundaries
                    aPn[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m - 1][n])
        return aPw, aPe, aPs, aPn
