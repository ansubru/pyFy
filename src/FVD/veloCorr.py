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
################################################################################MODULE FOR PRESSURE CORRECTION################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class veloCorr(object):

    """A class for discretizing relevant nodes and returning cell and face coefficients

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def velcorrU(self,matU,matPP, aPu):
        """A function that corrects velocity by pressure straddling
           Note! the face pressures are replaced by grid pressures (W, E, N, S): as pw - pe = PW/2 + PP /2 - PP/2 + PE/2
        """
        u = 1.0*matU
        PP = 1.0*matPP
        aP = 1.0 * aPu
        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        dx = IO_obj.x_dis # x-grid spacing
        dy = IO_obj.y_dis # y-grid spacing

        i = np.size(matU, 0)
        j = np.size(matU, 1)

        ppW, ppE, uPfinal = 0.0*matPP, 0.0*matPP, 0.0*matU

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if n == 0:  # First column
                    ppW[m][n] = PP[m][n]
                    ppE[m][n] = PP[m][n+1]
                    uPfinal[m][n] = u[m][n] + (dy * (ppW[m][n] - ppE[m][n]) / (dx*aP[m][n]))

                if n == j-1:  # Last column
                    ppW[m][n] = PP[m][n-1]
                    ppE[m][n] = PP[m][n]
                    uPfinal[m][n] = u[m][n] + (dy * (ppW[m][n] - ppE[m][n]) / (dx*aP[m][n]))

                else: # all other elements
                    ppW[m][n] = PP[m][n-1]
                    ppE[m][n] = PP[m][n+1]
                    uPfinal[m][n] = u[m][n] + (dy * (ppW[m][n] - ppE[m][n]) / (dx*aP[m][n]))

        return uPfinal

    def velcorrV(self, matV, matPP, aPv):
        """A function that corrects velocity by pressure straddling
           Note! the face pressures are replaced by grid pressures (W, E, N, S): as pw - pe = PW/2 + PP /2 - PP/2 + PE/2
        """
        u = 1.0 * matV
        PP = 1.0 * matPP
        aP = 1.0 * aPv
        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        dx = IO_obj.y  # y-grid spacing
        dy = IO_obj.x  # x-grid spacing

        i = np.size(matV, 0)
        j = np.size(matV, 1)

        ppS, ppN, uPfinal = 0.0 * matPP, 0.0 * matPP, 0.0 * matV

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if m == 0:  # First row
                    ppS[m][n] = PP[m+1][n]
                    ppN[m][n] = PP[m][n]
                    uPfinal[m][n] = u[m][n] + (dy * (ppS[m][n] - ppN[m][n]) / (dx*aP[m][n]))

                if m == i-1:  # Last row
                    ppS[m][n] = PP[m][n]
                    ppN[m][n] = PP[m-1][n]
                    uPfinal[m][n] = u[m][n] + (dy * (ppS[m][n] - ppN[m][n]) / (dx*aP[m][n]))

                else:  # all other elements
                    ppS[m][n] = PP[m+1][n]
                    ppN[m][n] = PP[m-1][n]
                    uPfinal[m][n] = u[m][n] + (dy * (ppS[m][n] - ppN[m][n]) / (dx*aP[m][n]))

        return uPfinal











