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
class Corr2(object):

    """A class for correcting face mass fluxes and velocity using pressure straddling

    Attributes:
       U --> Feed in U grid velocities (x-direction)
       v --> Feed in V grid velocity (y-direction)
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def velcorr(self,matU,matV, matPP, aPu):
        """A function that corrects velocity by pressure straddling (aPu is coeff of Pprime equation)
           Note! the face pressures are replaced by grid pressures (W, E, N, S): as pw - pe = PW/2 + PP /2 - PP/2 + PE/2
        """
        u = 1.0*matU
        v= 1.0*matV
        PP = 1.0*matPP
        aP = 1.0 * aPu
        ###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        dx = IO_obj.dx # x-grid spacing
        dy = IO_obj.dy # y-grid spacing

        # Boundary conditions
        uA = IO_obj.UA
        uB = IO_obj.UB
        uC = IO_obj.UC
        uD = IO_obj.UD
        vA = IO_obj.VA
        vB = IO_obj.VB
        vC = IO_obj.VC
        vD = IO_obj.VD

        i = np.size(u, 0)
        j = np.size(v, 1)

        ppW, ppE,ppN, ppS, uPfinal,vPfinal  = 0.0*matPP, 0.0*matPP, 0.0*matU, 0.0*matU, 0.0*matU, 0.0*matU

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes

                    ppW[m][n] = PP[m][n-1]
                    ppE[m][n] = PP[m][n+1]
                    ppS[m][n] = PP[m+1][n]
                    ppN[m][n] = PP[m-1][n]
                    uPfinal[m][n] = u[m][n] + (dy[m][n] * (ppW[m][n] - ppE[m][n]) / (2.0*aP[m][n]))
                    vPfinal[m][n] = v[m][n] + (dx[m][n] * (ppS[m][n] - ppN[m][n]) / (2.0* aP[m][n]))

            # PAD BCS
        # Apply bcs U and V
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0):  # A Boundary
                    uPfinal[m][n] = uA  # pad bc velocities
                    vPfinal[m][n] = vA  # pad bc velocities
                if (m >= 0 and n == j - 1):  # C Boundary
                    uPfinal[m][n] = uC  # pad bc velocities
                    vPfinal[m][n] = vC  # pad bc velocities
                if (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                    uPfinal[m][n] = uD  # pad bc velocities
                    vPfinal[m][n] = vD  # pad bc velocities
                if (m == 0 and n != j - 1 and n != 0):  # B Boundary
                    uPfinal[m][n] = uB  # pad bc velocities
                    vPfinal[m][n] = vB  # pad bc velocities

        return uPfinal, vPfinal

    def massFluxcorr(self,p,mdote, mdotn, aEpp, aNpp):
        """A function that corrects face velocity fluxes by pressure straddling
                """
        mdoteNew = 0.0 * mdote; mdotnNew = 0.0 * mdote;
        PPE = 0.0 * mdote;  PPN = 0.0 * mdote;
        i = np.size(mdote, 0)
        j = np.size(mdote, 1)

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes

                    PPE= p[m][n+1] ; PPN = p[m-1][n] ; PPP = p[m][n]

                    mdoteNew[m][n] = mdote[m][n] + aEpp[m][n] * (PPP - PPE)  # East face mass flux
                    mdotnNew[m][n] = mdotn[m][n] + aNpp[m][n] * (PPP - PPN)  # North face mass flux

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    mdoteNew[m][n] = 0.0  # East face

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    mdotnNew[m][n] = 0.0   # North face

        return mdoteNew, mdotnNew













