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
class rhieChow(object):

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
    def rcInterp(self, matU, matP):
        """A function that returns corrected face velocities and fluxes"""

        def rhicW(PE, PP, PW, PWW):
            coeff1w = (PE - 3.0 * PP + 3.0 * PW - PWW)
            return coeff1w

        def rhicE(PEE, PE, PP, PW):
            coeff1e = (PEE - 3.0 * PE + 3.0 * PP - PW)
            return coeff1e

        def rhicS(PN, PP, PS, PSS):
            coeff1s = (PN - 3.0 * PP + 3.0 * PS - PSS)
            return coeff1s

        def rhicN(PNN, PN, PP, PS):
            coeff1n = (PNN - 3.0 * PN + 3.0 * PP - PS)
            return coeff1n

        u = 1.0*matU
        Px = 1.0*matP
        #Get BC's
        from IO import IO
        IO_obj = IO("random")
        grid = IO_obj.grid_nodes  # generate a base grid node layout
        NU = IO_obj.nu  # viscosity nu (constant viscosity model)
        rho = IO_obj.rho  # Density in Kg/m3
        mu = NU * rho  # constant for all nodes
        dx = IO_obj.x_dis  # x-grid spacing
        dy = IO_obj.y_dis  # y-grid spacing

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
        Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aW, aE, aN, aS, aP1notmod, aP, A, Bx, By = disc_obj.FOU_disc(u, Px, UA, UB, UC, UD)

        i = np.size(matU, 0)
        j = np.size(matU, 1)

        pcorrw, pcorre, pcorrn, pcorrs = 0.0 * matU, 0.0 * matU, 0.0 * matU, 0.0 * matU
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns
                if(m==n==0): #First cell in the top left edge

                    PW = PWW = Px[m][n]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == 1): #Second cell bordering left edge

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n == 1): #Second row second cell

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == 1 and ((i-1)-m) > 1): #Second column

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == 1 and ((i-1)-m) == 1): #element before bottom element of second column

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m+1][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n != 1 and ((j-1)-n) > 1): # Cells on bordering bc 'B' sans 2 edge cells

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n != 1 and ((j-1)-n) == 1): # cell before right edge

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+1]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == (j-1)): # Right edge cell

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n == 0): # second cell in the first column

                    PW = PWW = Px[m][n]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n == 1): # second cell in second row

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n != 1 and ((j-1)-n) > 1): # second row without edges

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n != 1 and ((j-1)-n) == 1): # cell before second row edge cell

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n+1]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == 0 and ((i-1)-m) > 1): # Cells on bordering bc 'A' sans 2 edge cells

                    PW = PWW = Px[m][n]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == 0 and ((i-1)-m) == 1): # Cell before bottom left edge

                    PW = PWW = Px[m][n]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m+1][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-2) and n != 1 and ((j-1)-n) > 1): # last but one row without edges

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m+1][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-2) and n != 1 and ((j-1)-n) == 1): # last but one row, last but one column

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+1]
                    PS = PSS = Px[m+1][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-1) and n == 0): #  Bottom left edge cell

                    PW = PWW = Px[m][n]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-1) and n == 1): #  Cell to the right of bottom left edge cell -- BC 'D'

                    PW = PWW = Px[m][n-1]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-1) and n != 1 and ((j-1)-n) > 1): #  Cells bordering BC 'D' without sans two edge cells

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = PSS = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-1) and n != 1 and ((j-1)-n) == 1): #  Cell before bottom right edge

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n+1]
                    PS = PSS = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == (j-2) and ((i-1)-m) > 1): #  last but one column

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n+1]
                    PS = Px[m-1][n]
                    PSS = Px[m-2][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-1) and n == (j-1)): # Right edge cell

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = PSS = Px[m][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == (i-2) and n == (j-1)): # Cell above right edge cell

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = PSS = Px[m+1][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == (j-1) and ((i-1)-m) > 1): # Cells bordering BC 'C' sans two edges

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 1 and n == (j-1) and ((i-1)-m) == 1): # Cell before top right edge

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n == (j-1) and ((i-1)-m) == 1): # Cell before top right edge

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = PEE = Px[m][n]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = PNN = Px[m-1][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m != 1 and n != (j-1) and ((i-1)-m) > 1 and ((j-1)-m) > 1): # cells in the interior (two cell gap)

                    PW = Px[m][n-1]
                    PWW = Px[m][n-2]
                    PP = Px[m][n]
                    PE = Px[m][n+1]
                    PEE = Px[m][n+2]
                    PS = Px[m+1][n]
                    PSS = Px[m+2][n]
                    PN = Px[m-1][n]
                    PNN = Px[m-2][n]

#West faces
                    coeff1w = rhicW(PE, PP, PW, PWW)   # Pw = Pp (zero gradient bc)
                    coeff2w = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w/coeff2w
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e = rhicE(PEE, PE, PP, PW)
                    coeff2e = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e/coeff2e
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s = rhicS(PN, PP, PS, PSS)  # Pn = Pp (zero gradient bc)
                    coeff2s = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s/coeff2s
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n = rhicN(PNN, PN, PP, PS)   # Pn = Pp (zero gradient bc)
                    coeff2n = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n/coeff2n
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

        return ufw, ufe, ufn, ufs, pcorrw, pcorre, pcorrn, pcorrs










