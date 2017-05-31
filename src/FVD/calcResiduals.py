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

    def calcRes(self, var ,aW, aE, aN, aS, SU, aP):
        """Function to calculate residuals given a variable"""

        i = np.size(var, 0)
        j = np.size(var, 1)

        Rvar = 0

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    sumvar = aE[m][n]*var[m][n+1] + aW[m][n]*var[m][n-1] + aN[m][n]*var[m-1][n] + aS[m][n]*var[m+1][n] + SU[m][n] -aP[m][n]*var[m][n]
                    Rvar = Rvar + abs(sumvar)

        return (Rvar)

    def calcResomg(self, var ,aW, aE, aN, aS, SU, aP):
        """Function to calculate residuals given a variable"""

        i = np.size(var, 0)
        j = np.size(var, 1)

        Rvar = 0

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != 1 and n != 1 and m != (i - 1) and n != (j - 1) and m != (i - 2) and n != (j - 2)):  # special only for omega:
                    sumvar = aE[m][n]*var[m][n+1] + aW[m][n]*var[m][n-1] + aN[m][n]*var[m-1][n] + aS[m][n]*var[m+1][n] + SU[m][n] -aP[m][n]*var[m][n]
                    Rvar = Rvar + abs(sumvar)

        return (Rvar)
