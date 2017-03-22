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
class interpToFace(object):

    """A class for applying rhie-chow interpolation

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       p --> Feed in pressure
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""

    def faceInterp(self, matU, matP):
        """Class to interpolate to faces"""



        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        grid = IO_obj.grid_nodes #generate a base grid node layout
        u = np.array(1.0 * matU)
        p = np.array(1.0 * matP)
        #Boundary conditions
        UA = IO_obj.UA
        UB = IO_obj.UB
        UC = IO_obj.UC
        UD = IO_obj.UD
        VA = IO_obj.VA
        VB = IO_obj.VB
        VC = IO_obj.VC
        VD = IO_obj.VD


        from Discretize import Discretize
        disc_obj = Discretize()

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        # Get velocities

        Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS, aP, aPmod, A, Bx, By = disc_obj.FOU_disc(u,p, UA, UB, UC, UD)


        i = np.size(matU, 0)
        j = np.size(matU, 1)

        aEe, aWw, aNn, aSs, aPe = 0.0*u, 0.0*u, 0.0*u, 0.0*u ,0.0*u  # interpolated co-eff to faces
        print aPmod

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m == 0 and n == 0):  # cells in the top left edge

                    # West faces
                    aWw[m][n] = aW[m][n]

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m][n], aS[m + 1][n])

                    # North faces
                    aNn[m][n] = aN[m][n]

                    # aP term
                    aPe = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                elif (m == 0 and n != (j - 1)):  # First row bordering the BC --> B (sans edges)

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m][n], aS[m + 1][n])

                    # North faces
                    aNn[m][n] = aN[m][n]

                    # aP term
                    aPe = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                elif (m == 0 and n == (j - 1)):  # Top right edge

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = aE[m][n]

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m][n], aS[m + 1][n])

                    # North faces
                    aNn[m][n] = aN[m][n]

                    # aP term
                    print m,n
                    aPe[m][n] = aPmod[m][n]

                elif (m == (i - 1) and n == 0):  # bottom left edge

                    # West faces
                    aWw[m][n] = aW[m][n]

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = aS[m][n]

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aN[m][n], aN[m - 1][n])

                    # aP term
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                elif (m == (i - 1) and n != (j - 1)):  # Bottom row bordering the BC --> D  (sans edges)

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = aS[m][n]

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aN[m][n], aN[m - 1][n])

                    # aP term
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                elif (m != (i - 1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                    # West faces
                    aWw[m][n] = aW[m][n]

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m + 1][n], aS[m][n])

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aN[m][n], aN[m - 1][n])

                    # aP term
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

                elif (m == (i - 1) and n == (j - 1)):  # Bottom right edge

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = aE[m][n]

                    # South faces
                    aSs[m][n] = aS[m][n]

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aNn[m][n], aNn[m - 1][n])

                    # aP term
                    aPe[m][n] = aPmod[m][n]


                elif (m != (i - 1) and n == (j - 1)):  # Right (last) column bordering the BC --> C  (sans edges)

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = aE[m][n]

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m][n], aS[m + 1][n])

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aN[m][n], aN[m - 1][n])

                    # aP term
                    aPe[m][n] = aPmod[m][n]

                elif (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # All other elements

                    # West faces
                    aWw[m][n] = Interp_obj.lin_interp(aW[m][n - 1], aW[m][n])

                    # East faces
                    aEe[m][n] = Interp_obj.lin_interp(aE[m][n], aE[m][n + 1])

                    # South faces
                    aSs[m][n] = Interp_obj.lin_interp(aS[m][n], aS[m + 1][n])

                    # North faces
                    aNn[m][n] = Interp_obj.lin_interp(aN[m][n], aN[m - 1][n])

                    # aP term
                    aPe[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n + 1])

            return aWw, aEe, aSs, aNn, aPe