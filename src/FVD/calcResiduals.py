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
################################################################################MODULE FOR CALCULATING RESIDUALS################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint



#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class calcResiduals(object):

    """A class to calculate residuals

    Attributes:

    """

    def __init__(self):
        """Return object"""

    def calcRes(self, U, V, mdotw, mdote, mdotn, mdots, P,b):
        """Function to calculate residuals"""

        ###############------------CREATE IO OBJECT--------------###############


        from Discretize2 import Discretize2
        disc_obj = Discretize2()

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        from IO import IO
        IO_obj = IO("random")

        # Boundary conditions (for u-velocity)
        UA = IO_obj.UA
        UB = IO_obj.UB
        UC = IO_obj.UC
        UD = IO_obj.UD

        # Boundary conditions (for v-velocity)
        VA = IO_obj.VA
        VB = IO_obj.VB
        VC = IO_obj.VC
        VD = IO_obj.VD

        # Get co-efficients

        aW, aE, aN, aS, aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp = disc_obj.FOU_disc2(U, V, mdotw, mdote, mdotn, mdots, P)

        i = np.size(P, 0)
        j = np.size(P, 1)

        RU = 0.0
        RV = 0.0
        RP = 0.0
        RB = 0.0

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    sumU = aE[m][n]*U[m][n+1] + aW[m][n]*U[m][n-1] + aN[m][n]*U[m-1][n] + aS[m][n]*U[m+1][n] + SUxmod[m][n] -aPmod[m][n]*U[m][n]
                    sumV = aE[m][n] * V[m][n + 1] + aW[m][n] * V[m][n - 1] + aN[m][n] * V[m - 1][n] + aS[m][n] * V[m + 1][n] + SUymod[m][n] - \
                           aPmod[m][n] * V[m][n]
                    sumP = aEpp[m][n] * P[m][n + 1] + aWpp[m][n] * P[m][n - 1] + aNpp[m][n] * P[m - 1][n] + aSpp[m][n] * P[m + 1][n] + b[m][n] - \
                           aPpp[m][n] * P[m][n]
                    sumB = RB + b[m][n]



        RU = abs(RU + sumU)
        RV = abs(RV + sumV)
        RP = abs(RP + sumP)
        RB = abs(RB + sumB)

        return RU, RV, RP, RB
