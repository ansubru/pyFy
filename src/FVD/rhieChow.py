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
        Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs,aP, aPmod, A, Bx, By = disc_obj.FOU_disc(u, Px, UA, UB, UC, UD)

        i = np.size(matU, 0)
        j = np.size(matU, 1)
        coeff1w, coeff1e, coeff1s, coeff1n = 0.0 * matU, 0.0 * matU, 0.0 * matU, 0.0 * matU
        coeff2w, coeff2e, coeff2s, coeff2n = 0.0 * matU, 0.0 * matU, 0.0 * matU, 0.0 * matU
        pcorrw, pcorre, pcorrn, pcorrs = 0.0 * matU, 0.0 * matU, 0.0 * matU, 0.0 * matU
        for m in range(i): #loop through rows
            for n in range(j): #loop through columns
                if(m==n==0): #First cell in the top left edge
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m][n]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n+2] - (3.0*Px[m][n+1]) + (3.0*Px[m][n]) - Px[m][n-1])
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == 1): #Second cell bordering left edge (sans edges)
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n-1]) - Px[m][n-1]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n+2] - (3.0*Px[m][n+1]) + (3.0*Px[m][n]) - Px[m][n-1]) # Pw = Pp (zero gradient bc)
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n != 1 and n < (j-2) ): #First row bordering the BC --> B (sans 2 edge cells on either sides)
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n-1]) - Px[m][n-2]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n+2] - (3.0*Px[m][n+1]) + (3.0*Px[m][n]) - Px[m][n]) # Pw = Pp (zero gradient bc)
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == (j-2) ): #Cell before right edge
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n-1]) - Px[m][n-2]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n+1] - (3.0*Px[m][n+1]) + (3.0*Px[m][n]) - Px[m][n-1]) # Pw = Pp (zero gradient bc)
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == (j-1) ): #Cell at right edge
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n-1]) - Px[m][n-2]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m][n-1]) # Pw = Pp (zero gradient bc)
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]

                elif(m == 0 and n == 0 ): #Second cell in the first column
#West faces
                    coeff1w[m][n] = (Px[m][n+1] - (3.0*Px[m][n]) + (3.0*Px[m][n-1]) - Px[m][n-2]) # Pw = Pp (zero gradient bc)
                    coeff2w[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorrw[m][n] = coeff1w[m][n]/coeff2w[m][n]
                    ufw[m][n] = ufw[m][n] + pcorrw[m][n]
#East faces
                    coeff1e[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m][n-1]) # Pw = Pp (zero gradient bc)
                    coeff2e[m][n] = dy/(dx*4.0*aP[m][n])
                    pcorre[m][n] = coeff1e[m][n]/coeff2e[m][n]
                    ufe[m][n] = ufe[m][n] + pcorre[m][n]
#South faces
                    coeff1s[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m+1][n]) - Px[m+2][n]) # Pn = Pp (zero gradient bc)
                    coeff2s[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrs[m][n] = coeff1s[m][n]/coeff2s[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrs[m][n]
#North faces
                    coeff1n[m][n] = (Px[m][n] - (3.0*Px[m][n]) + (3.0*Px[m][n]) - Px[m+1][n]) # Pn = Pp (zero gradient bc)
                    coeff2n[m][n] = dx/(dy*4.0*aP[m][n])
                    pcorrn[m][n] = coeff1n[m][n]/coeff2n[m][n]
                    ufs[m][n] = ufs[m][n] + pcorrn[m][n]














