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

    def gridInterpU2(self, ufe):
        """Function to interpolate u velocities to grid nodes"""
        from IO import IO
        IO_obj = IO("random")
        uA = IO_obj.UA
        uB = IO_obj.UB
        uC = IO_obj.UC
        uD = IO_obj.UD
        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        ue = 1.0*ufe

        i = np.size(ue, 0)
        j = np.size(ue, 1)

        uPfinal = 0.0*ue

        # Apply bcs U
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0):  # A Boundary
                    uPfinal[m][n] = uA  # pad bc velocities
                if (m == 0 and n != j - 1 and n != 0):  # B Boundary
                    uPfinal[m][n] = uB  # pad bc velocities
                if (m >= 0 and n == j - 1):  # C Boundary
                    uPfinal[m][n] = uC  # pad bc velocities
                if (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                    uPfinal[m][n] = uD  # pad bc velocities

        # Interpolate U face velocities to nodes
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m > 1 and m <= i - 3 and n > 1 and n <= j - 3):  # Internal nodes (except boundary + first grid nodes)
                    uPfinal[m][n] = Interp_obj.lin_interp(ue[m][n],ue[m][n-1])
                if(m > 0 and m < i-1 and n == j-2): # Boundary face (EAST):  # first grid nodes
                    uPfinal[m][n] = Interp_obj.lin_interp(ue[m][n], 0.0) #wall face velocity zero
                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    uPfinal[m][n] = Interp_obj.lin_interp(ue[m][n], 0.0)  # wall face velocity zero
                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    uPfinal[m][n] = Interp_obj.lin_interp(ue[m][n], 0.0)  # wall face velocity zero
                if (m == i-2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    uPfinal[m][n] = Interp_obj.lin_interp(ue[m][n], 0.0)  # wall face velocity zero
        return uPfinal

    def gridInterpV2(self, vfn):
        """Function to interpolate u velocities to grid nodes"""
        from IO import IO
        IO_obj = IO("random")
        vA = IO_obj.VA
        vB = IO_obj.VB
        vC = IO_obj.VC
        vD = IO_obj.VD
        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        vn = 1.0*vfn

        i = np.size(vn, 0)
        j = np.size(vn, 1)

        vPfinal = 0.0*vn

        # Apply bcs V
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0):  # A Boundary
                    vPfinal[m][n] = vA  # pad bc velocities
                elif (m == 0 and n != j - 1 and n != 0):  # B Boundary
                    vPfinal[m][n] = vB  # pad bc velocities
                elif (m >= 0 and n == j - 1):  # C Boundary
                    vPfinal[m][n] = vC  # pad bc velocities
                elif (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                    vPfinal[m][n] = vD  # pad bc velocities

        # Interpolate V face velocities to nodes
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m > 1 and m <= i - 3 and n > 1 and n <= j - 3):  # Internal nodes (except boundary + first grid nodes)
                    vPfinal[m][n] = Interp_obj.lin_interp(vn[m][n],vn[m-1][n])
                if(m > 0 and m < i-1 and n == j-2): # Boundary face (EAST):  # first grid nodes
                    vPfinal[m][n] = Interp_obj.lin_interp(vn[m][n], 0.0) #wall face velocity zero
                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    vPfinal[m][n] = Interp_obj.lin_interp(vn[m][n], 0.0)  # wall face velocity zero
                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    vPfinal[m][n] = Interp_obj.lin_interp(vn[m][n], 0.0)  # wall face velocity zero
                if (m == i-2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    vPfinal[m][n] = Interp_obj.lin_interp(vn[m][n], 0.0)  # wall face velocity zero
        return vPfinal
