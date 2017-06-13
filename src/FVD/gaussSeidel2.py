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
from multiprocessing import Pool


#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
def computeError(u, uold):
    """Computes absolute error using an L2 norm for the solution.
    This requires that self.u and self.old_u must be appropriately
    setup."""
    v = abs(u - uold)/max(uold.any(), 1e-12)
    sumRes = v.sum()
    return sumRes

class gaussSiedel2(object):

    def __init__(self):
        """Return object"""


    def gaussSeidel3u(self, u, aW, aE, aN, aS, aP, SUx, iters, flag):
        """Solves gauss seidel for a fixed number of iterations."""

        iter = 0 # Mock up error parameter only to run the While

        while iters > iter:

            if flag in ['omega']:

                u[2:-2, 2:-2] = (u[2:-2, 1:-3] * aW[2:-2, 2:-2] + u[2:-2, 3:-1] * aE[2:-2, 2:-2] + u[3:-1, 2:-2] * aS[2:-2, 2:-2] + \
                             u[1:-3, 2:-2] * aN[2:-2, 2:-2] + SUx[2:-2, 2:-2]) / aP[2:-2, 2:-2]
            else:
                u[1:-1, 1:-1] = (u[1:-1, 0:-2] * aW[1:-1, 1:-1] + u[1:-1, 2:] * aE[1:-1, 1:-1] + u[2:, 1:-1] * aS[1:-1,1:-1] + \
                                 u[0:-2, 1:-1] * aN[1:-1, 1:-1] + SUx[1:-1, 1:-1]) / aP[1:-1, 1:-1]
            iter += 1

        print "Solved G-S for %s with %i iterations"%(flag,iter)
        return u

    def gaussSeidel4u(self, u, aW, aE, aN, aS, aP, SUx, resinp, flag):
        """Solves gauss seidel for a fixed number of iterations."""

        iter = 0 # Mock up error parameter only to run the While
        diffs = 1
        residual = 1
        resHist  = [diffs]

        while residual > resinp:
            uold = u
            if flag in ['omega']:

                u[2:-2, 2:-2] = (u[2:-2, 1:-3] * aW[2:-2, 2:-2] + u[2:-2, 3:-1] * aE[2:-2, 2:-2] + u[3:-1, 2:-2] * aS[2:-2, 2:-2] + \
                             u[1:-3, 2:-2] * aN[2:-2, 2:-2] + SUx[2:-2, 2:-2]) / aP[2:-2, 2:-2]
            else:
                u[1:-1, 1:-1] = (u[1:-1, 0:-2] * aW[1:-1, 1:-1] + u[1:-1, 2:] * aE[1:-1, 1:-1] + u[2:, 1:-1] * aS[1:-1,1:-1] + \
                                 u[0:-2, 1:-1] * aN[1:-1, 1:-1] + SUx[1:-1, 1:-1]) / aP[1:-1, 1:-1]

            sumRes = computeError(u, uold)
            residual = sumRes #/ np.size(u)  # update residual
            print residual
            # resHist.append((residual))
            # diffs = abs(resHist[-1] - resHist[-2])
            if iter > 5000 :
                break

            iter += 1
        print "Solved G-S for %s with %i iterations and a final residual of %e"%(flag, iter, diffs)

        return u

    def gaussSeidel5u(self, u, aW, aE, aN, aS, aP, SUx, resinp, flag):
        """Solves gauss seidel for a fixed number of iterations."""

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")

        # Initialize all relevant variables:u, v, p etc.
        i = np.size(u, 0)
        j = np.size(u, 1)
        ugs = 1.0 * u
        uold = 0.0 * u

        iter = 0  # Mock up error parameter only to run the While
        diffs = 1
        sumRes = 1
        resHist = [diffs]

        while diffs > resinp:
            # while 50 > iter:
            if flag in ['u', 'U', 'v', 'V', 'p', 'P', 'k', 'K']:

                for m in range(i):  # loop through rows
                    for n in range(j):  # loop through columns

                        if (m != 0 and n != 0 and m != (i - 1) and n != (
                            j - 1)):  # Internal nodes (for u, v and k):

                            uold[m][n] = ugs[m][n]
                            uN = u[m - 1][n]
                            uS = u[m + 1][n]
                            uE = u[m][n + 1]
                            uW = u[m][n - 1]

                            ugs[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                                         uN * aN[m][n] + SUx[m][n]) / aP[m][n]
                            sumRes = sumRes + (abs((uold[m][n] - ugs[m][n])) / max(1e-16, abs(uold[m][n])))

            else:
                for m in range(2, i - 2):  # loop through rows
                    for n in range(2, j - 2):  # loop through columns

                        uold[m][n] = u[m][n]
                        uN = u[m - 1][n]
                        uS = u[m + 1][n]
                        uE = u[m][n + 1]
                        uW = u[m][n - 1]

                        ugs[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                                     uN * aN[m][n] + SUx[m][n]) / aP[m][n]
                        sumRes = sumRes + (abs((uold[m][n] - ugs[m][n])) / max(1e-16, abs(uold[m][n])))
            iter += 1
            u = ugs  # update u
            residual = sumRes / np.size(u)  # update residual
            resHist.append((residual))
            diffs = abs(resHist[-1] - resHist[-2])

            if iter > 5000:
                break

        if flag in ['U', 'u']:
            print "Solved G-S for U with %i iterations and a final residual of %e" % (iter, diffs)
        elif flag in ['v', 'V']:
            print "Solved G-S for V with %i iterations and a final residual of %e" % (iter, diffs)
        elif flag in ['k', 'K']:
            print "Solved G-S for K with %i iterations and a final residual of %e" % (iter, diffs)
        else:
            print "Solved G-S for omega with %i iterations and a final residual of %e" % (iter, diffs)

        return u

    def gaussSeidel6u(self, u, aW, aE, aN, aS, aP, SUx, resinp, flag):
        """Solves gauss seidel for a fixed number of iterations."""

        iter = 0 # Mock up error parameter only to run the While
        # Initialize all relevant variables:u, v, p etc.
        i = np.size(u, 0)
        j = np.size(u, 1)
        diffs = 1
        residual = 1
        resHist  = [diffs]

        while residual > resinp:
            if flag in ['omega']:
                for m in range(2, i - 2):  # loop through rows
                    for n in range(2, j - 2):  # loop through columns

                        uN = u[m - 1][n]
                        uS = u[m + 1][n]
                        uE = u[m][n + 1]
                        uW = u[m][n - 1]

                        u[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                                     uN * aN[m][n] + SUx[m][n]) / aP[m][n]

            else:
                for m in range(i):  # loop through rows
                    for n in range(j):  # loop through columns
                        if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes:

                            uN = u[m - 1][n]
                            uS = u[m + 1][n]
                            uE = u[m][n + 1]
                            uW = u[m][n - 1]

                            u[m][n] = (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                                       uN * aN[m][n] + SUx[m][n]) / aP[m][n]

            if iter > 5000 :
                break

            iter += 1
        print "Solved G-S for %s with %i iterations and a final residual of %e"%(flag, iter, diffs)

        return u

