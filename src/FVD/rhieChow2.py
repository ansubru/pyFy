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
    def rcInterp2(self, matU, matV, matP):
        """A function that returns corrected face velocities and fluxes"""

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
        mu = NU * rho  # constant for all nodes
        dx = IO_obj.dx  # x-grid spacing
        dy = IO_obj.dy  # y-grid spacing

        # Boundary conditions
        UA = IO_obj.UA
        UB = IO_obj.UB
        UC = IO_obj.UC
        UD = IO_obj.UD
        VA = IO_obj.VA
        VB = IO_obj.VB
        VC = IO_obj.VC
        VD = IO_obj.VD

        #obtain face velocities using the Discretize class obj
        from Discretize import Discretize
        disc_obj = Discretize()
        Fe, Fw, Fn, Fs, ufe, ufw, vfn, vfs, aW, aE, aN, aS, aWp, aEp, aNp, aSp, anotmodP, aP, SUxmod, SUymod, A, Bx, By = disc_obj.FOU_discdisc_obj2.FOU_disc2( v, UA, UB, UC, UD, v, VA, VB, VC, VD , Px)

        i = np.size(matU, 0)
        j = np.size(matU, 1)

        pcorre, pcorrn, = 0.0 * matU, 0.0 * matU

        for m in range(i): #loop through rows
            for n in range(j): #loop through columns

                if (m >= 1 or m <= i-3 and n >= 1 or n <= j-3):  # Internal nodes (except boundary + first grid nodes)

                    PP = Px[m][n]
                    PW = Px[m][n - 1]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]
                    PS = Px[m + 1][n]

#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]

#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    vfs[m][n] = vfs[m][n] + pcorrn[m][n]

                elif(m > 0 and n == j-2): # Boundary face (EAST)
                    PP = Px[m][n]
                    PW = Px[m][n - 1]
                    PE = Px[m][n + 1]
                    PEE = Px[m][n + 1]

                    # East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy / (dx * 4.0 * aP[m][n])
                    pcorre[m][n] = coeff1e / coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]

                elif (m == 1 and n > 0):  # Boundary face (North)
                    PP = Px[m][n]
                    PN = Px[m - 1][n]
                    PNN = Px[m - 1][n]
                    PS = Px[m + 1][n]

                    # North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)  # Pn = Pp (zero gradient bc)
                    coeff2n = dx / (dy * 4.0 * aP[m][n])
                    pcorrn[m][n] = coeff1n / coeff2n
                    vfs[m][n] = vfs[m][n] + pcorrn[m][n]

        return ufe, vfn, pcorre, pcorrn










