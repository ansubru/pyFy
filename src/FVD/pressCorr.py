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
class pressCorr(object):

    """A class for discretizing relevant nodes and returning cell and face coefficients

    Attributes:
       U --> Feed in U velocity (x-direction)
       v --> Feed in V velocity (y-direction)
       Specify U or V boundary values
    """

    def __init__(self):
        """Return object"""


 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###
    def pcorr(self,matU, matP):
        """A function that returns co-efficients for P prime equation"""

        def assignCoeffs(aWw,aEe,aSs,aNn):
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
        dx = IO_obj.x_dis # x-grid spacing
        dy = IO_obj.y_dis # y-grid spacing
        alpha = IO_obj.alpha #under-relaxation

##############-------CREATE INTERPOLATE OBJECT----------################

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
###############------------------------------------------################
#Initialize all relevant variables:u, v, p etc.
        u = 1.0*matU
        Px = 1.0*matP
        aE, aW, aN, aS, aP = 0.0*u, 0.0*u, 0.0*u, 0.0*u , 0.0*u  # initialize co-effs
        A = []

        #Obtain face interpolated coeffs
        from interpToFace import interpToFace
        interpFc_obj = interpToFace()

        aWw, aEe, aSs, aNn = interpFc_obj.faceInterp(u,Px) # All these are aP's interpolated to faces


        i = np.size(matU, 0)
        j = np.size(matU, 1)

        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if (m == 0 and n == 0):  # cells in the top left edge

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif m == 0 and n == (j-1):  # Top right edge

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]
                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m+1][n] = aS[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif m == (i-1) and n == 0:  # bottom left edge

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)

                   aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                   aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                   #'A' matrix generation
                   amat = np.array([])
                   tempU = np.array(0.0*matU)
                #Assigning coeff to grid nodes
                   tempU[m][n] = aP[m][n]
                   tempU[m][n+1] = aE[m][n]
                   tempU[m][n-1] = aW[m][n]
                   tempU[m-1][n] = aN[m][n]
                   for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                   A.append(amat)

                elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
#Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)


                elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
                    #Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

                elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                    aW[m][n], aE[m][n], aS[m][n], aN[m][n] = assignCoeffs(aWw[m][n], aEe[m][n], aSs[m][n], aNn[m][n])
                    aP[m][n] = aW[m][n] + aE[m][n] + aS[m][n] + aN[m][n]

                    #'A' matrix generation
                    amat = np.array([])
                    tempU = np.array(0.0*matU)
#Assigning coeff to grid nodes
                    tempU[m][n] = aP[m][n]
                    tempU[m][n-1] = aW[m][n]
                    tempU[m][n+1] = aE[m][n]
                    tempU[m+1][n] = aS[m][n]
                    tempU[m-1][n] = aN[m][n]
                    for a in range(i):
                        for b in range(j):
                            amat = np.append(amat,tempU[a][b])
                    A.append(amat)

        return aW, aE, aN, aS, aP
