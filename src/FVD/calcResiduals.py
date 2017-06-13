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

    def calcRes(self, var ,aW, aE, aN, aS, SU, aP, B, flag):
        """Function to calculate residuals given a variable"""

        sumvar = 0.0 * var
        if flag in ['omega']:
            sumvar[2:-2, 2:-2] = abs((var[2:-2, 1:-3] * aW[2:-2, 2:-2] + var[2:-2, 3:-1] * aE[2:-2, 2:-2] + var[3:-1, 2:-2] * aS[2:-2, 2:-2] + \
                                  var[1:-3, 2:-2] * aN[2:-2, 2:-2] + SU[2:-2, 2:-2]) - aP[2:-2, 2:-2]*var[2:-2, 2:-2])
        elif flag in ['B', 'b']:
            sumvar = B
        else:
            sumvar[1:-1, 1:-1] = abs((var[1:-1, 0:-2] * aW[1:-1, 1:-1] + var[1:-1, 2:] * aE[1:-1, 1:-1] + var[2:, 1:-1] * aS[1:-1,1:-1] + \
                                              var[0:-2, 1:-1] * aN[1:-1, 1:-1] + SU[1:-1, 1:-1]) - aP[1:-1, 1:-1]*var[1:-1, 1:-1])
        Rvar = sumvar.sum()
        return Rvar

    def calcResomg(self, var ,aW, aE, aN, aS, SU, aP):
        """Function to calculate residuals given a variable"""

        i = np.size(var, 0)
        j = np.size(var, 1)

        Rvar = 0

        for m in range(2,i-2):  # loop through rows
            for n in range(2,j-2):  # loop through columns

                    sumvar = aE[m][n]*var[m][n+1] + aW[m][n]*var[m][n-1] + aN[m][n]*var[m-1][n] + aS[m][n]*var[m+1][n] + SU[m][n] - aP[m][n]*var[m][n]
                    Rvar = Rvar + abs(sumvar)

        return (Rvar)

    def calcResB(self, B):
        """Function to calculate residuals given B"""

        i = np.size(B, 0)
        j = np.size(B, 1)

        Rvar = 0

        for m in range(1,i-1):  # loop through rows
            for n in range(1,j-1):  # loop through columns
                    sumvar = B[m][n]
                    Rvar = Rvar + abs(sumvar)
        return (Rvar)
