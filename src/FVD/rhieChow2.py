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
################################################################################MODULE FOR RHIE-CHOW INTERPOLATION################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class rhieChow2(object):

    """A class for applying rhie-chow interpolation

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       p --> Feed in pressure
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""



 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def rcInterp2(self, matU, matV, mdotw, mdote, mdotn, mdots, matP):
        """A function that provides necessary corrections for face mass fluxes"""

        def rhicE(PEE, PE, PP, PW):
            coeff1e = (PEE - 3.0 * PE + 3.0 * PP - PW)
            return coeff1e

        def rhicN(PNN, PN, PP, PS):
            coeff1n = (PNN - 3.0 * PN + 3.0 * PP - PS)
            return coeff1n

        u = 1.0*matU

        v = 1.0*matV

        Px = 1*matP

        #Get BC's
        from IO import IO
        IO_obj = IO("random")
        NU = IO_obj.nu  # viscosity nu (constant viscosity model)
        rho = IO_obj.rho  # Density in Kg/m3
        dx = IO_obj.dx  # x-grid spacing
        dy = IO_obj.dy  # y-grid spacing

        # NOTE : RHIE CHOW uses aPmod interpolated to the faces
        # Obtain face interpolated coeffs
        from Discretize2 import Discretize2
        disc_obj2 = Discretize2()
        aW, aE, aN, aS, aWw, aEe, aNn, aSs, aP, aPmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp = disc_obj2.FOU_disc2(u, v, mdotw, mdote, mdotn, mdots, Px)

        i = np.size(matU, 0)
        j = np.size(matU, 1)

        pcorre, pcorrn, = 0.0 * matU, 0.0 * matU
        coeff1eMat = 0.0 * matU
        coeff1nMat = 0.0 * matU

        for m in range(i): #loop through rows
            for n in range(j): #loop through columns

                if (m >= 1 and m < i-2 and n > 0 and n < j-2):  #East faces
                    PP = Px[m][n]
                    PW = Px[m][n - 1]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]

                    coeff1e = (PEE - 3.0 * PE + 3.0 * PP - PW)
                    coeff1eMat[m][n] = coeff1e
                    coeff2e = (rho*dy[m][n]*dy[m][n])/(4.0*aEe[m][n])
                    pcorre[m][n] = coeff1e*coeff2e

                if (m > 1 and m <= i-2 and n > 0 and n <= j-2):  #North faces
                    PP = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]
                    PS = Px[m + 1][n]

                    coeff1n = (PNN - 3.0 * PN + 3.0 * PP - PS)   # Pn = Pp (zero gradient bc)
                    coeff1nMat[m][n] = coeff1n
                    coeff2n = (rho*dx[m][n]**2)/(4.0*aNn[m][n])
                    pcorrn[m][n] = coeff1n*coeff2n

        return pcorre, pcorrn










