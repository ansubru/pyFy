#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = "Ananda S. Kannan"
__copyright__ = "Copyright 2017, pyFy project"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ananda S. Kannan"
__email__ = "ansubru@gmail.com"
__credits to__ https://github.com/stangmechanic/NE155_Homework_3/blob/master/GaussSeidel.py

"""
################################################################################MODULE FOR GAUSS SEIDEL METHOD################################################################################################

import os
import numpy as np
import re
import sys
import pprint
import scipy.sparse as sps
import scipy.sparse.linalg as spla

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###

class gaussSiedel2(object):

    def __init__(self):
        """Return object"""

    def gaussSeidel2u(self, matU, aWinp, aEinp, aNinp, aSinp, aPinp, SUxin, iterinp):
        """Solves gauss seidel for a fixed number of iterations."""
        print("Trying to solve equations for U using the Gauss-Seidel2 method")

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")
        #Boundary conditions
        UA = IO_obj.UA
        UB = IO_obj.UB
        UC = IO_obj.UC
        UD = IO_obj.UD

        #Initialize all relevant variables:u, v, p etc.
        u = 1.0*matU
        usolve, uP, uW, uE, uN, uS = 0.0*matU, 0.0*matU, 0.0*matU, 0.0*matU,0.0*matU, 0.0*matU
        aW = np.array(1.0*aWinp)
        aE = np.array(1.0*aEinp)
        aS = np.array(1.0*aSinp)
        aN = np.array(1.0*aNinp)
        aP = np.array(1.0*aPinp)
        SUx = np.array(1.0*SUxin)
        i = np.size(matU, 0)
        j = np.size(matU, 1)

        iter = 0 # Mock up error parameter only to run the While

        while iter < iterinp:
            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == 0 and n == 0):  # cells in the top left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = UA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = UB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = UB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == 0 and n == (j-1):  # Top right edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = UC
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = UB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == (i-1) and n == 0:  # bottom left edge
                        uP[m][n] = u[m][n]
                        uW[m][n] = UA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = UD
                        uN[m][n] = u[m-1][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)
                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = UD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = UA
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = UC
                       uS[m][n] = UD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]


                    elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = UC
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n]+ SUx[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUx[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

            iter += 1

            return usolve

    def gaussSeidel2v(self, matV, aWinp, aEinp, aNinp, aSinp, aPinp, SUyin , iterinp):
        """Solves gauss seidel for a fixed number of iterations."""
        print("Trying to solve equations for V using the Gauss-Seidel2 method")

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")
        #Boundary conditions
        VA = IO_obj.VA
        VB = IO_obj.VB
        VC = IO_obj.VC
        VD = IO_obj.VD

        #Initialize all relevant variables:u, v, p etc.
        u = 1.0*matV
        usolve, uP, uW, uE, uN, uS = 0.0*matV, 0.0*matV, 0.0*matV, 0.0*matV,0.0*matV, 0.0*matV
        aW = np.array(aWinp)
        aE = np.array(aEinp)
        aS = np.array(aSinp)
        aN = np.array(aNinp)
        aP = np.array(aPinp)
        SUy = np.array(SUyin)
        i = np.size(matV, 0)
        j = np.size(matV, 1)

        iter = 0 # Mock up error parameter only to run the While

        while iter > iterinp:
            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == 0 and n == 0):  # cells in the top left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = VA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]

                    elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]

                    elif m == 0 and n == (j-1):  # Top right edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = VC
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == (i-1) and n == 0:  # bottom left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = VA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = VD
                        uN[m][n] = u[m-1][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = VD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = VA
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = VC
                       uS[m][n] = VD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]


                    elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = VC
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

            iter =+ 1

            return usolve

    def gaussSeidelP(self, matP, aWinp, aEinp, aNinp, aSinp, aPinp, Bin , iterinp):
        """Solves gauss seidel for a fixed number of iterations."""
        print("Trying to solve equations for P using the Gauss-Seidel2 method")

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")

        #Initialize all relevant variables:u, v, p etc.
        u = 1.0*matP
        usolve, uP, uW, uE, uN, uS = 0.0*matP, 0.0*matP, 0.0*matP, 0.0*matP,0.0*matP, 0.0*matP
        aW = np.array(aWinp)
        aE = np.array(aEinp)
        aS = np.array(aSinp)
        aN = np.array(aNinp)
        aP = np.array(aPinp)
        b = np.array(Bin)
        i = np.size(matP, 0)
        j = np.size(matP, 1)

        iter = 0 # Mock up error parameter only to run the While

        while iter > iterinp:
            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == 0 and n == 0):  # cells in the top left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                        print usolve[m][n]

                    elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]

                    elif m == 0 and n == (j-1):  # Top right edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == (i-1) and n == 0:  # bottom left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m][n]
                        uN[m][n] = u[m-1][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n]
                       uS[m][n] = u[m][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]


                    elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

        iter =+ 1
        print usolve
        return usolve

