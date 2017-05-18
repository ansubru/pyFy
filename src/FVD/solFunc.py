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
################################################################################MODULE FOR SOME ELEMENTARY SOLVER FUNCTIONS################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class solFunc(object):

    """A class for performing basic solver operations

    Attributes:
      calculating coninuity B term
      under relaxing P
      setting pressure at particular point
      Applying bcs for U and V
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###

    def applyBCs(self,U,UA,UB,UC,UD):
        """A function to Apply BC's"""
        i = np.size(U, 0)
        j = np.size(U, 1)
        u = 0.0*U

        # Apply bcs U and V
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0):  # A Boundary
                    u[m][n] = UA  # pad bc velocities
                if (m >= 0 and n == j - 1):  # C Boundary
                    u[m][n] = UC  # pad bc velocities
                if (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                    u[m][n] = UD  # pad bc velocities
                if (m == 0 and n != j - 1 and n != 0):  # B Boundary
                    u[m][n] = UB  # pad bc velocities
        return u

    def applyBCmus(self,kinit, omegainit, mut,muk, muomega):
        """A function to Apply turbulent viscosity BC's"""
        from IO import IO
        IO_obj = IO("random")
        sigmak = IO_obj.sigmak
        sigmaomega = IO_obj.sigmaomega
        cw2 = IO_obj.cw2
        i = np.size(mut, 0)
        j = np.size(mut, 1)
        kloc = 0.0*kinit
        omegaloc = 0.0*omegainit
        mutloc = 0.0*mut
        mukloc = 0.0 * mut
        muomegaloc = 0.0*mut

        # Apply bcs U and V
        for m in range(0, j):
            for n in range(0, i):
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    kloc[m][n] = kinit
                    omegaloc[m][n] = omegainit
                    mutloc[m][n] = kloc[m][n]/omegaloc[m][n]
                    mukloc[m][n] = mutloc[m][n]/sigmak
                    muomegaloc[m][n] = mutloc[m][n] / sigmaomega
        return u

    def setPress(self,P,x,y):
        """A function to set a fixed pressure at a grid point x,y"""
        i = np.size(P, 0)
        j = np.size(P, 1)
        refP = P[x][y]
        for m in range(i):  # loop through rows
            for n in range(j):

                P[m][n] = P[m][n] - refP

        return P

    def setPsetbcs(self,P):
        """"""
        i = np.size(P, 0)
        j = np.size(P, 1)
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m == 0:
                    P[m][n] = P[m+1][n]
                if n == 0:
                    P[m][n] = P[m][n+1]
                if m == i-1:
                    P[m][n] = P[m-1][n]
                if n == j-1:
                    P[m][n] = P[m][n-1]
        return P

    def rlxP(self, P, PP, coeff):
        """A function to under relax P"""
        i = np.size(P, 0)
        j = np.size(P, 1)
        pPnew = 0.0 * P
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                pPnew[m][n] = P[m][n] + coeff * PP[m][n]
        return pPnew

    def calcFaceMassFlux(self, u, v):
        """A function to calc face mass fluxes"""
        from IO import IO
        IO_obj = IO("random")
        rho = IO_obj.rho  # Density in Kg/m3
        dx = IO_obj.dx  # x-grid spacing
        dy = IO_obj.dy  # y-grid spacing
        from Interpolate import Interpolate
        Interp_obj = Interpolate()
        i = np.size(u, 0)
        j = np.size(u, 1)
        mfe, mfn = 0.0 * u, 0.0 * u

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    mfe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n + 1]) * dy[m][n] * rho  # East face
                    mfn[m][n] = Interp_obj.lin_interp(v[m][n], v[m - 1][n]) * dx[m][n] * rho  # North face

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    mfe[m][n] = 0.0  # East face


                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    mfn[m][n] = 0.0   # North face
        return  mfe, mfn

    def calcFaceMassFluxIW(self, u, v):
        """A function to calc face mass fluxes using interpolation weights"""
        from IO import IO
        IO_obj = IO("random")
        rho = IO_obj.rho  # Density in Kg/m3
        dx = IO_obj.dx  # x-grid spacing
        dy = IO_obj.dy  # y-grid spacing
        # Interpolation weights
        fxe = IO_obj.fxe
        fyn = IO_obj.fyn

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
        i = np.size(u, 0)
        j = np.size(u, 1)
        mfe, mfn = 0.0 * u, 0.0 * u

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    mfe[m][n] = Interp_obj.weighted_interp(u[m][n], u[m][n + 1], fxe[m][n]) * dy[m][n] * rho  # East face
                    mfn[m][n] = Interp_obj.weighted_interp(v[m][n], v[m - 1][n], fyn[m][n]) * dx[m][n] * rho  # North face

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    mfe[m][n] = 0.0  # East face

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    mfn[m][n] = 0.0   # North face
        return  mfe, mfn

    def correctFaceMassFlux(self, mfe, mfn, pcorre, pcorrn):
        """A function to correct face velocities"""
        i = np.size(mfe, 0)
        j = np.size(mfe, 1)

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    mfe[m][n] = mfe[m][n] + pcorre[m][n]
                    mfn[m][n] = mfn[m][n] + pcorrn[m][n]

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    mfe[m][n] = 0.0  # East face

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    mfn[m][n] = 0.0   # North face

        return  mfe, mfn

    def getFaceMassFluxWS(self, mfe, mfn):
        """A function to return west and south mass fluxes (given east and north fluxes)"""
        from Interpolate import Interpolate
        Interp_obj = Interpolate()
        i = np.size(mfe, 0)
        j = np.size(mfn, 1)
        mfw, mfs = 0.0 * mfe, 0.0 * mfe
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes

                    mfw[m][n] = mfe[m][n-1]
                    mfs[m][n] = mfn[m+1][n]

                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    mfw[m][n] = 0.0  # West face

                if (m == i - 2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    mfs[m][n] = 0.0  # South face

        return mfw, mfs

    def calcBterm(self, mdotw, mdote, mdotn, mdots):
        """A function to calculate the continuity error source needed to solve the pressure correction equation"""
        # Calculate b
        b = 0.0*mdotw
        i = np.size(b, 0)
        j = np.size(b, 1)
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    b[m][n] = mdotw[m][n] - mdote[m][n] + mdots[m][n] - mdotn[m][n]
        return b
