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
    def FOU_disc2(self, matU, matV, mdotw, mdote, mdotn, mdots, matP):
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

        def assignCoeffs(aWw,aEe,aSs,aNn,rho, dx, dy):
            # West faces
            aW = (rho*dy**2)/aWw

            # East faces
            aE = (rho*dy**2)/aEe

            # South faces
            aS = (rho*dx**2)/aSs

            # North faces
            aN = (rho*dx**2)/aNn

            return aW, aE, aS, aN


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

        # Boundary conditions
        uA = IO_obj.UA
        uB = IO_obj.UB
        uC = IO_obj.UC
        uD = IO_obj.UD
        vA = IO_obj.VA
        vB = IO_obj.VB
        vC = IO_obj.VC
        vD = IO_obj.VD

##############-------CREATE INTERPOLATE OBJECT----------################

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
###############------------------------------------------################
#Initialize all relevant variables:u, v, p etc.
        u = 1.0*matU
        v = 1.0*matV
        Px = 1.0*matP
        dPx = 0.0*matU; dPy = 0.0*matU; Fe = 0.0*matU; Fw= 0.0*matU;  Fn= 0.0*matU; Fs = 0.0*matU;


#Calculate all relevant co-efficients
        aP = 0.0*matU; aPmod = 0.0*matU; aW = 0.0*matU; aWp = 0.0*matU; aE = 0.0*matU; aEp= 0.0*matU; aN= 0.0*matU; aNp= 0.0*matU; aS= 0.0*matU; aSp = 0.0*matU;
        aWpp = 0.0 * matU; aEpp = 0.0 * matU; aSpp = 0.0 * matU; aNpp = 0.0 * matU; aPpp = 0.0 * matU;

        SUx = 0.0*matU; SUxmod = 0.0*matU; SPx = 0.0*matU; SUy = 0.0*matU; SUymod = 0.0*matU; SPy = 0.0*matU;

        DxE = 0.0*matU; DxW = 0.0*matU; DyN = 0.0*matU; DyS = 0.0*matU;

        i = np.size(matU,0)
        j = np.size(matU,1)


#Diffusion conductance
        for m in range(0, j):
            for n in range(0, i):
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    DxE[m][n] = mu / delXE[m][n]
                    DxW[m][n] = mu / delXW[m][n]
                    DyN[m][n] = mu / delYN[m][n]
                    DyS[m][n] = mu / delYS[m][n]

    # Apply bcs U and V
        for m in range(0, j):
            for n in range(0, i):
                if (m >= 0 and n == 0): # A Boundary
                    u[m][n] = uA  # pad bc velocities
                    v[m][n] = vA  # pad bc velocities
                elif (m == 0 and n != j-1 and n != 0): # B Boundary
                    u[m][n] = uB  # pad bc velocities
                    v[m][n] = vB  # pad bc velocities
                elif (m >= 0 and n == j-1): # C Boundary
                    u[m][n] = uC  # pad bc velocities
                    v[m][n] = vC  # pad bc velocities
                elif (m == i-1  and n != j-1 and n != 0): # D Boundary
                    u[m][n] = uD  # pad bc velocities
                    v[m][n] = vD  # pad bc velocities

        A = []
    #Initialize face velocities
        ufe, ufw, vfn, vfs = 0.0 * u, 0.0 * u, 0.0 * u, 0.0 * u  # u face velocities at each grid node

    #Calculate face velocities from mass fluxes (mdot)
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  #Internal nodes
                    ufw[m][n] =  mdotw[m][n]/(rho*dy[m][n]) # West face
                    ufe[m][n] =  mdote[m][n]/(rho*dy[m][n]) # East face
                    vfn[m][n] =  mdotn[m][n]/(rho*dx[m][n]) # North face
                    vfs[m][n] =  mdots[m][n]/(rho*dx[m][n]) # South face

    # Discretize convection diffusion equation (CD - FOU) -- !!!!!!!!!!!!!!!COEFFICIENTS FOR U AND V EQUATIONS!!!!!!!!!!!!
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns

                if(m != 0 and n != 0 and m != (i-1) and n != (j-1)): #Internal nodes
#source terms
                    dPx[m][n] = -0.5*(Interp_obj.CD_interp(Px[m][n-1], Px[m][n+1], dx[m][n]))
                    dPy[m][n] = -0.5*(Interp_obj.CD_interp(Px[m+1][n], Px[m-1][n], dy[m][n]))
                    SUx[m][n] = - dPx[m][n]*dx[m][n]*dy[m][n]
                    SUy[m][n] = - dPy[m][n]*dy[m][n]*dx[m][n]
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
                    aP[m][n] = (DxW[m][n]*dy[m][n] + DxE[m][n]*dy[m][n] + DyS[m][n]*dx[m][n] + DyN[m][n]*dx[m][n]) + max(0.0, -Fw[m][n])*dy[m][n] \
                                + max(0.0, Fe[m][n])*dy[m][n] + max(0.0, -Fs[m][n])*dx[m][n] + max(0.0, Fn[m][n])*dx[m][n] - SPx[m][n]

#Under-relaxation (due to non-linearity in PDE's

                    aPmod[m][n] = aP[m][n]/alpha
                    SUxmod[m][n] = SUx[m][n] + (aPmod[m][n]*(1-alpha)*u[m][n])
                    SUymod[m][n] = SUy[m][n] + (aPmod[m][n]*(1-alpha)*v[m][n])

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

        #Interpolate apmod to faces for pressure correction equation and rhiechow
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes

                    # West faces
                    aWp[m][n] = Interp_obj.lin_interp(aPmod[m][n-1], aPmod[m][n])

                    # East faces
                    aEp[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m][n+1])

                    # South faces
                    aSp[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m+1][n])

                    # North faces
                    aNp[m][n] = Interp_obj.lin_interp(aPmod[m][n], aPmod[m-1][n])

                #HACK FOR HAVING NON ZERO VALUES AT BC NODES

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    aEp[m][n] = Interp_obj.lin_interp(aPmod[m][n], 0)  # Neumann bc at all boundaries

                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    aWp[m][n] = Interp_obj.lin_interp(aPmod[m][n], 0)  # Neumann bc at all boundaries

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    aNp[m][n] = Interp_obj.lin_interp(aPmod[m][n], 0)  # Neumann bc at all boundaries

                if (m == i - 2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    aSp[m][n] = Interp_obj.lin_interp(aPmod[m][n], 0)  # Neumann bc at all boundaries

        # !!!!!!!!!!!!!!!COEFFICIENTS FOR Pprime EQUATION# !!!!!!!!!!!!

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns

                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    aWpp[m][n] = (rho * dy[m][n] ** 2) / aWp[m][n]
                    aEpp[m][n] = (rho * dy[m][n] ** 2) / aEp[m][n]
                    aSpp[m][n] = (rho * dx[m][n] ** 2) / aSp[m][n]
                    aNpp[m][n] = (rho * dx[m][n] ** 2) / aNp[m][n]

                if (m > 0 and m < i - 1 and n == j - 2):  # Boundary face (EAST):  # first grid nodes
                    aEpp[m][n] = 0.0

                if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
                    aWpp[m][n] = 0.0

                if (m == 1 and n > 0 and n < j - 1):  # Boundary face (NORTH):  # first grid nodes
                    aNpp[m][n] = 0.0

                if (m == i - 2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
                    aSpp[m][n] = 0.0

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                    aPpp[m][n] = aWpp[m][n] + aEpp[m][n] + aSpp[m][n] + aNpp[m][n]


        return  aW, aE, aN, aS,aWp, aEp, aNp, aSp, aP, aPmod, SUxmod, SUymod, aWpp, aEpp, aNpp, aSpp, aPpp













