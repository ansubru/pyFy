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
    def velcorr(self,matU,matP,matPP):
        """A function that returns co-efficients of discretization using first order upwind"""
        u = 1.0*matU
        Px = 1.0*matP
        PP = 1.0*matPP
        unw, une, uns, unn = 0.0*matU, 0.0*matU, 0.0*matU, 0.0*matU

        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.x_dis # x-grid spacing
        dy = IO_obj.y_dis # y-grid spacing
        alpha = IO_obj.alpha #under-relaxation

        #Obtain face interpolated coeffs
        from interpToFace import interpToFace
        interpFc_obj = interpToFace()

        #Obtain corrected face velocities
        from rhieChow import rhieChow
        rc_obj = rhieChow()
        ufwrc, uferc, ufnrc, ufsrc, pcorrw, pcorre, pcorrn, pcorrs = rc_obj.rcInterp(u,Px)

        aW, aE, aS, aN = interpFc_obj.faceInterp(u,Px) # coe-effs for P prime equation

        i = np.size(matU, 0)
        j = np.size(matU, 1)

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m == 0 and n == 0):  # cells in the top left edge

                    if(m==n==0): #cells in the top left edge

                        #East faces
                        unw[m][n] = ufwrc[m][n] + aw[m][n]*(PP[m][n] - PP[m][n])

                        #West faces
                        une[m][n] = uferc[m][n] + aE[m][n]*(PP[m][n] - PP[m][n+1])

