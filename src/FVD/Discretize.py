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
    def FOU_disc(self, mat, uA, uB, uC, uD):
        """A function that returns co-efficients of discretization using first order upwind"""

###############------------CREATE IO OBJECT--------------################

        from IO import IO
        IO_obj = IO("random")
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.x_dis # x-grid spacing
        dy = IO_obj.y_dis # y-grid spacing

#Boundary conditions: Face fluxes
        FA = IO_obj.FA
        FB = IO_obj.FB
        FC = IO_obj.FC
        FD = IO_obj.FD
#Boundary conditions: Diffusion conductance
        DA = IO_obj.DA
        DB = IO_obj.DB
        DC = IO_obj.DC
        DD = IO_obj.DD
#Boundary conditions: Boundary velocities
        VA = UA = 1.0*uA #West boundary
        VB = UB = 1.0*uB #North boundary
        VC = UC = 1.0*uC #East boundary
        VD = UD = 1.0*uD #South boundary

##############-------CREATE INTERPOLATE OBJECT----------################

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
###############------------------------------------------################
#Initialize all relevant variables:u, v, p etc.
        u, Px, dPx, dPy, Fe, Fw, Fn, Fs = 1.0*mat, 1.0*mat,  1.0*mat,  1.0*mat,  1.0*mat, \
                                          1.0*mat,  1.0*mat, 1.0*mat
        ufe, ufw, ufn, ufs = u, u, u, u #interpolated face velocities at each grid node
#Calculate all relevant co-efficients
        aP, aW, aE, aN, aS, SUx, SPx, SUy, SPy = 1.0*mat,  1.0*mat,  1.0*mat,  1.0*mat,\
                                             1.0*mat,  1.0*mat,  1.0*mat,  1.0*mat, 1.0*mat,
        i = np.size(mat,0)
        j = np.size(mat,1)
        print i,j
#Diffusion conductance for an equidistant grid
        NU = IO_obj.nu #viscosity nu (constant viscosity model)
        rho = IO_obj.rho #Density in Kg/m3
        mu = NU*rho #constant for all nodes
        dx = IO_obj.x_dis # x-grid spacing

        Dx = (NU*rho)/IO.x_dis
        Dy = (NU*rho)/IO.y_dis
        coeff = 0.5 #Co-efficient for equidistant grid

        for n in range(i): #loop through rows
            for m in range(j): #loop through columns
                if(m==n==0): #cells in the top left edge
#source term
                    dPx[m][n] = (Interp_obj.CD_interp(Px[m][n], Px[m+1][n], dx))
                    dPy[m][n] = (Interp_obj.CD_interp(Px[m][n+1], Px[m][n], dy))
                    SUx[m][n] = ((2.0*DA + FA)*UA) + ((2.0*DB + FB)*UB) + dPx[m][n]
                    SUy[m][n] = ((2.0*DA + FA)*VA) + ((2.0*DB + FB)*VB) + dPy[m][n]
                    SPx[m][n] = SPy[m][n] = -((2.0*DA + FA) + (2.0*DB + FB))
#West faces
                    ufw[m][n] = UA
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = 0.0
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    ufn[m][n] = UB
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = 0.0
#aP term
                    aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                               + Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                               - SPx[m][n])
                elif(n == 0 and m != (j-1)): #First row bordering the BC --> B (sans edges)
#source terms
                    print m
                    dPx[m][n] = (Interp_obj.CD_interp(Px[m-1][n], Px[m+1][n], dx))*coeff
                    dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n+1], Px[m][n], dy))
                    SUx[m][n] = ((2.0*DB + FB)*UB) + dPx[m][n]
                    SUy[m][n] = ((2.0*DB + FB)*VB) + dPy[m][n]
                    SPx[m][n] = SPy[m][n] = -((2.0*DB + FB))
#West faces
                    ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m][n])
                    Fw[m][n] = rho*ufw[m][n]
                    aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    Fe[m][n] = rho*ufe[m][n]
                    aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    Fs[m][n] = rho*ufs[m][n]
                    aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    ufn[m][n] = UB
                    Fn[m][n] = rho*ufn[m][n]
                    aN[m][n] = 0.0
#aP term
                    aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                + Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                - SPx[m][n])
                #elif(m == 0 and n == (j-1)): #top right edge
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m-1][n], Px[m][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n+1], Px[m][n], dy))
                    #SUx[m][n] = ((2.0*DC + FC)*UC) + ((2.0*DB + FB)*UB) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DC + FC)*VC) + ((2.0*DB + FB)*VB) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DC + FC) + (2.0*DB + FB))
#West faces
                    #ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m],[n])
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    #ufe[m][n] = UC
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = 0.0
#South faces
                    #ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    #ufn[m][n] = UB
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = 0.0
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(n < (i-1) and m == 0): #First coulumn bordering the BC --> A (sans edges)
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m+1][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n+1], Px[m][n-1], dy))
                    #SUx[m][n] = ((2.0*DA + FA)*UA) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DA + FA)*VA) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DA + FA))
#West faces
                    #ufw[m][n] = UA
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = 0.0
#East faces
                    #ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    #ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(n == i and m == 0): #bottom left edge
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m+1][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n-1], dy))
                    #SUx[m][n] = ((2.0*DA + FA)*UA) + ((2.0*DD + FD)*UD) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DA + FA)*VA) + ((2.0*DD + FD)*VD) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DA + FA) + (2.0*DD + FD))
#West faces
                    #ufw[m][n] = UA
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = 0.0
#East faces
                    #ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    #ufs[m][n] = UD
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = 0.0
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(n == i and m < (j-1)): #Bottom row bordering the BC --> D  (sans edges)
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m-1][n], Px[m+1][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n-1], dy))
                    #SUx[m][n] = ((2.0*DD + FD)*UD) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DD + FD)*VD) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DD + FD))
#West faces
                    #ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m],[n])
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    #ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    #ufs[m][n] = UD
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = 0.0
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(n == (i-1) and m == (j-1)): #Bottom right edge
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m-1][n], Px[m][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n], Px[m][n-1], dy))
                    #SUx[m][n] = ((2.0*DD + FD)*UD) + ((2.0*DC + FC)*UC) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DD + FD)*VD) + ((2.0*DC + FC)*VC) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DD + FD) + (2.0*DC + FC))
#West faces
                    #ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m],[n])
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    #ufe[m][n] = UC
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = 0.0
#South faces
                    #ufs[m][n] = UD
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = 0.0
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(n < (i-1) and m == (j-1)): #Right column bordering the BC --> C  (sans edges)
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m-1][n], Px[m][n], dx))
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n+1], Px[m][n-1], dy))
                    #SUx[m][n] = ((2.0*DC + FC)*UC) + dPx[m][n]
                    #SUy[m][n] = ((2.0*DC + FC)*VC) + dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = -((2.0*DC + FC))
#West faces
                    #ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m],[n])
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    #ufe[m][n] = UC
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = 0.0
#South faces
                    #ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
                #elif(m != 0 and n != 0 and m != (j-1) and n != (i-1)): #All other elements
#source terms
                    #print m,n
                    #dPx[m][n] = coeff*(Interp_obj.CD_interp(Px[m-1][n], Px[m+1][n], dx))
                    #dPy[m][n] = coeff*(Interp_obj.CD_interp(Px[m][n+1], Px[m][n-1], dy))
                    #SUx[m][n] = 0.0 + dPx[m][n]
                    #SUy[m][n] = 0.0+ dPy[m][n]
                    #SPx[m][n] = SPy[m][n] = 0.0
#West faces
                    #ufw[m][n] = Interp_obj.lin_interp(u[m-1][n], u[m][n])
                    #Fw[m][n] = rho*ufw[m][n]
                    #aW[m][n] = Dx + max(0.0, Fw[m][n])
#East faces
                    #ufe[m][n] = Interp_obj.lin_interp(u[m][n], u[m+1][n])
                    #Fe[m][n] = rho*ufe[m][n]
                    #aE[m][n] = Dx - max(0.0, -Fe[m][n])
#South faces
                    #ufs[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n+1])
                    #Fs[m][n] = rho*ufs[m][n]
                    #aS[m][n] = Dy + max(0.0, Fs[m][n])
#North faces
                    #ufn[m][n] = Interp_obj.lin_interp(u[m][n], u[m][n-1])
                    #Fn[m][n] = rho*ufn[m][n]
                    #aN[m][n] = Dy - max(0.0, -Fn[m][n])
#aP term
                    #aP[m][n] = (aW[m][n] + aE[m][n] + aN[m][n] + aS[m][n]
                                #+ Fe[m][n] - Fw[m][n] + Fn[m][n] - Fs[m][n]
                                #- SPx[m][n])
        return  Px, dPx, dPy, Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aP, aW, aE, aN, aS, SUx, SPx, SUy, SPy










