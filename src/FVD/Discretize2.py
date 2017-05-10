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
################################################################################MODULE FOR DISCRETIZATION################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class Discretize2(object):

    """A class for discretizing relevant nodes and returning cell and face coefficients

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def FOU_disc2(self, matU, uA, uB, uC, uD, matV, vA, vB, vC, vD , matP):
        """A function that returns co-efficients of discretization using first order upwind"""

        def underlxaP(ap,coeff):
            """A function to under-relax ap"""
            aPchg = ap/coeff
            return aPchg

        def underlxSU(sux, coeff1,coeff2, vel):
            """A function to under-relax SU
              coeff1: aPmod
              coeff2: alpha
              """
            SUchg = sux + (coeff1*(1-coeff2)*vel)
            return SUchg

        def cloneRowMat(F, G, l, p):
            """A function to clone pth row from a matrix G to into the lth row of a matrix F (p > l)"""
            j = np.size(F, 1)
            for m in range(p):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == p - 1):
                        F[m][n] = G[l - 1][n]
            return F

        def cloneColMat(F, G, l, p):
            """A function to clone lth column into the pth one (p > l)"""
            Fnew = F
            j = np.size(Fnew, 0)
            for m in range(j):  # loop through rows
                for n in range(p):  # loop through columns
                    if (n == p - 1):
                        Fnew[m][n] = G[m][l - 1]
            return Fnew



###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.dx # x-grid spacing
        dy = IO_obj.dy # y-grid spacing
        alpha = IO_obj.alpha #under-relaxation
        delXW = IO_obj.delXW  # distance to Western neighbour
        delXE = IO_obj.delXE  # distance to Eastern neighbour
        delYN = IO_obj.delYN  # distance to Northern neighbour
        delYS = IO_obj.delYS  # distance to Northern neighbour

##############-------CREATE INTERPOLATE OBJECT----------################

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
###############------------------------------------------################
#Initialize all relevant variables:u, v, p etc.
        u = 1.0*matU
        v = 1.0*matV
        Px = 1.0*matP
        dPx, dPy, Fe, Fw, Fn, Fs =   0.0*matU,  0.0*matU,  0.0*matU, \
                                          0.0*matU,  0.0*matU, 0.0*matU

#Calculate all relevant co-efficients
        aP, aPmod, aW, aWp, aE, aEp, aN, aNp, aS, aSp = 0.0*matU, 0.0*matU, 0.0*matU,  0.0*matU,  0.0*matU,0.0*matU,\
                                                        0.0*matU, 0.0*matU,  0.0*matU,  0.0*matU
        SUx, SUxmod, SPx, SUy, SUymod, SPy, DxE, DxW, DyN , DyS = 0.0*matU, 0.0*matU, 0.0*matU,  0.0*matU,  0.0*matU,0.0*matU ,0.0*matU, 0.0*matU , 0.0*matU, 0.0*matU

        i = np.size(matU,0)
        j = np.size(matU,1)

#Diffusion conductance
        for m in range(0, j):
            for n in range(0, i):
                DxE[m][n] = mu / delXE[m][n]
                DxW[m][n] = mu / delXW[m][n]
                DyN[m][n] = mu / delYN[m][n]
                DyS[m][n] = mu / delYS[m][n]

    # Apply bcs U
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0): # A Boundary
                    u[m][n] = uA  # pad bc velocities
                elif (m == 0 and n != j-1 and n != 0): # B Boundary
                    u[m][n] = uB  # pad bc velocities
                elif (m >= 0 and n == j-1): # C Boundary
                    u[m][n] = uC  # pad bc velocities
                elif (m == i-1  and n != j-1 and n != 0): # D Boundary
                    u[m][n] = uD  # pad bc velocities

    # Apply bcs V
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0):  # A Boundary
                    v[m][n] = vA  # pad bc velocities
                elif (m == 0 and n != j - 1 and n != 0):  # B Boundary
                    v[m][n] = vB  # pad bc velocities
                elif (m >= 0 and n == j - 1):  # C Boundary
                    v[m][n] = vC  # pad bc velocities
                elif (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                    v[m][n] = vD  # pad bc velocities
        A = []
    #Initialize face velocities
        ufe, ufw, vfn, vfs = 0.0 * u, 0.0 * u, 0.0 * u, 0.0 * u  # interpolated u face velocities at each grid node

    #Calculate face velocities
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns
                if (m > 0 and n >= 2 ):  # west faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n - 1], u[m][n]) #West face
                elif (m > 0 and n <= j-3):
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n + 1]) #East face
                elif (m >= 2 and n>0):
                    vfn[m][n] = Interp_obj.lin_interp(v[m][n], v[m - 1][n])  # North face
                elif (m <= i-3 and n > 0):
                    vfs[m][n] = Interp_obj.lin_interp(v[m][n], v[m + 1][n]) #South face

    # Discretize convection diffusion equation (CD - FOU)
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns

                if(m != 0 and n != 0 and m != (i-1) and n != (j-1)): #Internal nodes
#source terms
                    dPx[m][n] = (Interp_obj.CD_interp(Px[m][n-1], Px[m][n+1], dx[m][n]))
                    dPy[m][n] = (Interp_obj.CD_interp(Px[m+1][n], Px[m-1][n], dy[m][n]))
                    SUx[m][n] = (0.0 + dPx[m][n])*dy[m][n]
                    SUy[m][n] = (0.0+ dPy[m][n])*dx[m][n]
                    SPx[m][n], SPy[m][n] = 0.0, 0.0
#West faces
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (DxW[m][n] + max(0.0, Fw[m][n]))*dy[m][n]
#East faces
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (DxE[m][n] + max(0.0, -Fe[m][n]))*dy[m][n]
#South faces
                    Fs[m][n] = rho*vfs[m][n]
                    aS[m][n] = (DyS[m][n] + max(0.0, Fs[m][n]))*dx[m][n]
#North faces
                    Fn[m][n] = rho*vfn[m][n]
                    aN[m][n] = (DyN[m][n] + max(0.0, -Fn[m][n]))*dx[m][n]
#aP term
                    aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                + Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                - SPx[m][n])

#Under-relaxation (due to non-linearity in PDE's

                    aPmod[m][n] = underlxaP(aP[m][n], alpha)
                    SUxmod[m][n] = underlxSU(SUx[m][n], aPmod[m][n], alpha, u[m][n])
                    SUymod[m][n] = underlxSU(SUy[m][n], aPmod[m][n], alpha, u[m][n])

#'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
#Assigning coeff to grid nodes
                    tempU[m][n] = aPmod[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

#coeff- Terms for p' equation (HACK!)
                    aWp[m][n] = (DxW[m][n] + max(0.0, Fw[m][n]))*dy[m][n]
                    aEp[m][n] = (DxE[m][n] + max(0.0, -Fe[m][n]))*dy[m][n]
                    aSp[m][n] = (DyS[m][n] + max(0.0, Fs[m][n]))*dx[m][n]
                    aNp[m][n] = (DyN[m][n] + max(0.0, -Fn[m][n]))*dx[m][n]

                else: #Boundary nodes __> These values should not be really used as they dont mean anything!!!
                    ufw[m][n] = 0.0; Fw[m][n] = rho*ufw[m][n]
                    ufe[m][n] = 0.0; Fe[m][n] = rho*ufe[m][n]
                    vfn[m][n] = 0.0; Fn[m][n] = rho*vfn[m][n]
                    vfs[m][n] = 0.0; Fs[m][n] = rho*vfs[m][n]

        #B matrix generation
        Bx = []
        By = []
        for l in range(np.size(matU,0)):
            for m in range(np.size(matU,1)):
                Bx.append(SUxmod[l][m])
                By.append(SUymod[l][m])
        return  Fe, Fw, Fn, Fs, ufe, ufw, vfn, vfs, aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, A, Bx, By













