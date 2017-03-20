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
class Discretize(object):

    """A class for discretizing relevant nodes and returning cell and face coefficients

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def FOU_disc(self, matU, matP, uA, uB, uC, uD):
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

###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.x_dis # x-grid spacing
        dy = IO_obj.y_dis # y-grid spacing
        alpha = IO_obj.alpha #under-relaxation

#Boundary conditions: Face fluxes
        FA = 1.0*IO_obj.FA
        FB = 1.0*IO_obj.FB
        FC = 1.0*IO_obj.FC
        FD = 1.0*IO_obj.FD
#Boundary conditions: Diffusion conductance
        DA = 1.0*IO_obj.DA
        DB = 1.0*IO_obj.DB
        DC = 1.0*IO_obj.DC
        DD = 1.0*IO_obj.DD
#Boundary conditions: Boundary velocities
        VA, UA = uA, uA #West boundary
        VB, UB = uB, uB #North boundary
        VC, UC = uC, uC #East boundary
        VD, UD = uD, uD #South boundary

##############-------CREATE INTERPOLATE OBJECT----------################

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
###############------------------------------------------################
#Initialize all relevant variables:u, v, p etc.
        u = 1.0*matU
        Px = 1.0*matP
        dPx, dPy, Fe, Fw, Fn, Fs =   0.0*matU,  0.0*matU,  0.0*matU, \
                                          0.0*matU,  0.0*matU, 0.0*matU
        ufe, ufw, ufn, ufs = 0.0*u, 0.0*u, 0.0*u, 0.0*u #interpolated face velocities at each grid node
#Calculate all relevant co-efficients
        aP, aPmod, aW, aE, aN, aS, SUx, SUxmod, SPx, SUy, SUymod, SPy = 0.0*matU, 0.0*matU, 0.0*matU,  0.0*matU,  0.0*matU,\
                                             0.0*matU,  0.0*matU,  0.0*matU,  0.0*matU, 0.0*matU, 0.0*matU, 0.0*matU
        i = np.size(matU,0)
        j = np.size(matU,1)
#Diffusion conductance for an equidistant grid
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.x_dis # x-grid spacing


        Dx = (NU*rho)/IO.x_dis
        Dy = (NU*rho)/IO.y_dis
        coeff = 0.5 #Co-efficient for equidistant grid
        A = []

        for m in range(i): #loop through rows
            for n in range(j): #loop through columns
                if(m==n==0): #cells in the top left edge
#source term
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n+1], dx))
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m][n], dy))
                    SUx[m][n] = (((2.0*DA + FA)*UA) + ((2.0*DB + FB)*UB) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DA + FA)*VA) + ((2.0*DB + FB)*VB) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DA + FA) + (2.0*DB + FB))*dy, -((2.0*DA + FA) + (2.0*DB + FB))*dx
#West faces
                    ufw[m][n] = UA
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = 0.0

#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx
#North faces
                    ufn[m][n] = UB
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = 0.0
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
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif(m == 0 and n != (j-1)): #First row bordering the BC --> B (sans edges)
#source terms
                    dPx[m][n] = (Interp_obj.CD_interp(Px[m][n-1], Px[m][n+1], dx))*coeff
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m][n], dy))
                    SUx[m][n] = (((2.0*DB + FB)*UB) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DB + FB)*VB) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DB + FB))*dy, -((2.0*DB + FB))*dx
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx
#North faces
                    ufn[m][n] = UB
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = 0.0
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
                    tempU[m][n+1] = aE[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)
                elif(m == 0 and n == (j-1)): #Top right edge
#source terms
                    dPx[m][n] = (Interp_obj.CD_interp(Px[m][n-1], Px[m][n], dx))*coeff
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m][n], dy))
                    SUx[m][n] = (((2.0*DB + FB)*UB) + ((2.0*DC + FC)*UC) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DB + FB)*VB) + ((2.0*DC + FC)*UC) +  dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DB + FB) + (2.0*DC + FC))*dy, -((2.0*DB + FB) + (2.0*DC + FC))*dx

#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = UC
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = 0.0
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx
#North faces
                    ufn[m][n] = UB
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = 0.0
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
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)
                elif(m == (i-1) and n == 0): #bottom left edge
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n+1], dx))

                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m-1][n], dy))
                    SUx[m][n] = (((2.0*DA + FA)*UA) + ((2.0*DD + FD)*UD) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DA + FA)*VA) + ((2.0*DD + FD)*VD) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DA + FA) + (2.0*DD + FD))*dy, -((2.0*DA + FA) + (2.0*DD + FD))*dx
#West faces
                    ufw[m][n] = UA
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = 0.0
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = UD
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = 0.0
#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
                    tempU[m][n+1] = aE[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)
                elif(m  == (i-1) and n != (j-1)): #Bottom row bordering the BC --> D  (sans edges)
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n-1], Px[m][n+1], dx))
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m-1][n], dy))
                    SUx[m][n] = (((2.0*DD + FD)*UD) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DD + FD)*VD) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DD + FD))*dy, -((2.0*DD + FD))*dx
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = UD
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = 0.0
#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
                    tempU[m][n+1] = aE[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)
                elif(m  != (i-1) and n == 0): #First column bordering the BC --> A  (sans edges)
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n+1], dx))
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m-1][n], dy))
                    SUx[m][n] = (((2.0*DA + FA)*UA) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DA + FA)*VA) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DA + FA))*dy, -((2.0*DA + FA))*dx
#West faces
                    ufw[m][n] = UA
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = 0.0
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m+1][n], u[m][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx

#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif(m == (i-1) and n == (j-1)): #Bottom right edge
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n-1], Px[m][n], dx))
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m-1][n], dy))
                    SUx[m][n] = (((2.0*DD + FD)*UD) + ((2.0*DC + FC)*UC) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DD + FD)*VD) + ((2.0*DC + FC)*VC) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DD + FD) + (2.0*DC + FC))*dy, -((2.0*DD + FD) + (2.0*DC + FC))*dx
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = UC
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = 0.0
#South faces
                    ufs[m][n] = UD
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = 0.0
#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif(m != (i-1) and n == (j-1)): #Right (last) column bordering the BC --> C  (sans edges)
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n-1], Px[m][n], dx))
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m-1][n], dy))
                    SUx[m][n] = (((2.0*DC + FC)*UC) + dPx[m][n])*dy
                    SUy[m][n] = (((2.0*DC + FC)*VC) + dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = -((2.0*DC + FC))*dy, -((2.0*DC + FC))*dx
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = UC
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = 0.0
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx
#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)
                elif(m != 0 and n != 0 and m != (i-1) and n != (j-1)): #All other elements
#source terms
                    dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n-1], Px[m][n+1], dx))
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m+1][n], Px[m-1][n], dy))
                    SUx[m][n] = (0.0 + dPx[m][n])*dy
                    SUy[m][n] = (0.0+ dPy[m][n])*dx
                    SPx[m][n], SPy[m][n] = 0.0, 0.0
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m][n-1], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = (Dx + max(0.0, Fw[m][n]))*dy
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = (Dx - max(0.0, -Fe[m][n]))*dy
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = (Dy + max(0.0, Fs[m][n]))*dx
#North faces
                    ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m-1][n])
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = (Dy - max(0.0, -Fn[m][n]))*dx
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
#B matrix generation
        Bx = []
        By = []
        for l in range(np.size(matU,0)):
            for m in range(np.size(matU,1)):
                Bx.append(SUxmod[l][m])
                By.append(SUymod[l][m])


        return  Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, A, Bx, By










