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
################################################################################MODULE FOR PLOTTING RESULTS################################################################################################

import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint
import matplotlib.pyplot as plt



#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class plotResults(object):

    """A class to plot results

    """

    def __init__(self):
        """Return object"""

    def plotContours(self,X, Y, U, V, P, Pp, it):
        y=[]
        Uplot=[]
        Vplot=[]
        for j in range(np.shape(X)[0]-1):
            for i in range(np.shape(X)[1]-1):
                if X[j,i]>0.45 and X[j,i]<0.55 and Y[j,i]<0.99:
                    y.append(Y[j,i])
                    Uplot.append(U[j,i])
                    Vplot.append(V[j,i])


        LarsData = np.loadtxt('LarsData.txt', skiprows=1)

        plt.figure(figsize=(30, 10))
        plt.subplot(231)
        plt.contourf(X, Y, U)
        plt.title('U at it: %d' % (it))
        plt.colorbar()
        plt.subplot(232)
        plt.contourf(X, Y, V)
        plt.title('V at it: %d' % (it))
        plt.colorbar()
        plt.subplot(233)
        plt.contourf(X, Y, P)
        plt.title('P at it: %d' % (it))
        plt.colorbar()
        plt.subplot(234)
        plt.contourf(X, Y, Pp)
        plt.title('Pp at it: %d' % (it))
        plt.colorbar()
        plt.subplot(235)
        plt.quiver(X, Y, U, V, scale=1)
        plt.title('Velocity vectors at it: %d' % (it))
        plt.subplot(236)
        plt.plot(LarsData[:, 0], LarsData[:, 2], 'ro-', label='lada U')
        plt.plot(LarsData[:, 1], LarsData[:, 2], 'go-', label='lada V')
        plt.plot(Uplot, y, 'bo-', label='my U')
        plt.plot(Vplot, y, 'ko-', label='my V')
        plt.legend()
        plt.ylabel('y [m]')
        plt.xlabel('velocity [m/s]')
        plt.title('Comparison with Lars\'s data it: %d' % (it))
        plt.show()

    def genGrid(self, matU):

        from IO import IO
        IO_obj = IO("random")
        grid_x = IO_obj.grid_x
        grid_y = IO_obj.grid_y
        return grid_x, grid_y

