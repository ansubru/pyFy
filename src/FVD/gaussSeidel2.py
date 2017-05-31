#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = "Ananda S. Kannan"
__copyright__ = "Copyright 2017, pyFy project"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ananda S. Kannan"
__email__ = "ansubru@gmail.com"
__credits to__ https://github.com/stangmechanic/NE155_Homework_3/blob/master/GaussSeidel.py

"""
################################################################################MODULE FOR GAUSS SEIDEL METHOD################################################################################################

import os
import numpy as np
import re
import sys
import pprint
import scipy.sparse as sps
import scipy.sparse.linalg as spla

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###

class gaussSiedel2(object):

    def __init__(self):
        """Return object"""

    def gaussSeidel3u(self, u, aW, aE, aN, aS, aP, SUx, iters):
        """Solves gauss seidel for a fixed number of iterations."""

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")

        #Initialize all relevant variables:u, v, p etc.
        i = np.size(u, 0)
        j = np.size(u, 1)
        ugs = 0.0*u
        uold = 0.0*u

        iter = 0 # Mock up error parameter only to run the While
        residual = 1
        sumRes = 0.0

        while iters > iter:

            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  #Internal nodes:

                        uold[m][n] = u[m][n]
                        uN = u[m - 1][n]
                        uS = u[m + 1][n]
                        uE = u[m][n + 1]
                        uW = u[m][n - 1]

                        ugs[m][n]=(uW* aW[m][n]+ uE*aE[m][n]+uS*aS[m][n]+\
                                           uN*aN[m][n]+SUx[m][n])/aP[m][n]
            iter += 1
            u = ugs #update u

        print "Solved G-S with %i iterations"%(iter)
        return u

    def gaussSeidel4u(self, u, aW, aE, aN, aS, aP, SUx, resinp, flag):
        """Solves gauss seidel for a fixed number of iterations."""

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")

        #Initialize all relevant variables:u, v, p etc.
        i = np.size(u, 0)
        j = np.size(u, 1)
        ugs = 1.0*u
        uold = 0.0*u

        iter = 0 # Mock up error parameter only to run the While
        diffs = 1
        sumRes = 0.0
        resHist  = [diffs]

        while diffs > resinp:

            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if flag in ['u', 'U', 'v', 'V', 'p', 'P', 'k', 'K']:
                        if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  #Internal nodes (for u, v and k):

                            uold[m][n] = ugs[m][n]
                            uN = u[m - 1][n]
                            uS = u[m + 1][n]
                            uE = u[m][n + 1]
                            uW = u[m][n - 1]

                            ugs[m][n]=(uW* aW[m][n]+ uE*aE[m][n]+uS*aS[m][n]+\
                                               uN*aN[m][n]+SUx[m][n])/aP[m][n]
                            sumRes = sumRes + (-(uold[m][n] - ugs[m][n]))

                    else:
                        if (m != 0 and n != 0 and m != 1 and n != 1 and m != (i - 1) and n != (j - 1) and m != (i - 2) and n != (j - 2)):  # special only for omega:

                            uold[m][n] = u[m][n]
                            uN = u[m - 1][n]
                            uS = u[m + 1][n]
                            uE = u[m][n + 1]
                            uW = u[m][n - 1]

                            ugs[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                                         uN * aN[m][n] + SUx[m][n]) / aP[m][n]
                            sumRes = sumRes + (abs(uold[m][n] - ugs[m][n])) / max(uold[m][n], 1e-10)
            iter += 1
            u = ugs #update u
            residual = sumRes/np.size(u) #update residual
            resHist.append((residual))
            diffs = abs(resHist[-1] - resHist[-2])

            if iter > 300 :
                break

        print "Solved G-S with %i iterations"%(iter)
        return u

