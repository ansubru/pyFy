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

    def plotContours(self,X, Y, U, V, P, resU, resV, resP, resB, it):
        y=[]
        Uplot=[]
        Vplot=[]
        iter = list(range(1, it + 2))

        i = np.size(X, 0)
        j = np.size(X, 1)

        UToPlot = np.zeros([i,j])
        VToPlot = np.zeros([i,j])
        PToPlot = np.zeros([i,j])
        YToPlot = np.zeros([i,j])

        for m in range(i):
            for n in range(j):
                UToPlot[m][n] = U[-m-1][n]
                VToPlot[m][n] = V[-m - 1][n]
                PToPlot[m][n] = P[-m - 1][n]
                YToPlot[m][n] = Y[-m - 1][n]

        for m in range(i):
            for n in range(j):
                if m>0 and n == 5:
                    y.append(Y[m,n])
                    Uplot.append(0.5*(U[m,n]+U[m,n+1]))
                    Vplot.append(0.5*(V[m,n]+V[m,n+1]))

        LarsData = np.loadtxt('LarsData.txt', skiprows=1)

        plt.figure(figsize=(30, 10))
        plt.subplot(231)
        plt.contourf(X, YToPlot, UToPlot)
        plt.title('U at it: %d' % (it))
        plt.colorbar()
        plt.subplot(232)
        plt.contourf(X, YToPlot, VToPlot)
        plt.title('V at it: %d' % (it))
        plt.colorbar()
        plt.subplot(233)
        plt.contourf(X, YToPlot, PToPlot)
        plt.title('P at it: %d' % (it))
        plt.colorbar()
        plt.subplot(234)
        plt.quiver(X, YToPlot, UToPlot, VToPlot, scale=1)
        plt.title('Velocity vectors at it: %d' % (it))
        plt.subplot(235)
        plt.semilogy(iter[1:-1], resU[1:-1], 'blue', label="U-velocity")
        plt.semilogy(iter[1:-1], resV[1:-1], 'black', label="V-velocity")
        plt.semilogy(iter[1:-1], resP[1:-1], 'red', label="Pressure")
        plt.semilogy(iter[1:-1], resB[1:-1], 'green', label="Continuity")
        plt.xlabel('Number of iterations')
        plt.ylabel('Residual')
        plt.grid(True)
        plt.grid(which='both')
        plt.subplot(236)
        plt.plot(LarsData[:, 0], LarsData[:, 2], 'ro', label='lada U')
        plt.plot(LarsData[:, 1], LarsData[:, 2], 'go', label='lada V')
        plt.plot(Uplot[0:-1], y[0:-1], 'b-', label='my U')
        plt.plot(Vplot[0:-1], y[0:-1], 'k-', label='my V')
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

    def plotdataTurb(self, U, V, k, omega, Y, mut, resU, resV, resK, resomega, it):

        i = np.size(U, 0)
        j = np.size(V, 1)

        y = []
        Uplot = []
        Vplot = []
        kplot = []
        omegaplot = []
        mutplot = []
        iter = list(range(1, it + 2))


        for m in range(i):
            for n in range(j):
                if m>0 and n == 5:
                    y.append(Y[m,n])
                    Uplot.append(0.5*(U[m,n]+U[m,n+1]))
                    Vplot.append(0.5*(V[m,n]+V[m,n+1]))
                    kplot.append(0.5 * (k[m, n] + k[m, n + 1]))
                    omegaplot.append(0.5 * (omega[m, n] + omega[m, n + 1]))
                    mutplot.append(0.5 * (mut[m, n] + mut[m, n + 1]))

        LarsData = np.loadtxt('LarsDataTurb.txt', skiprows=1)

        plt.figure
        plt.subplot(231)
        plt.plot(LarsData[:, 0], LarsData[:, 2], 'ro', label='lada U')
        plt.plot(Uplot[0:-1], y[0:-1], 'k-', label='my U')
        plt.ylabel('y [m]')
        plt.xlabel('velocity [m/s]')
        plt.subplot(232)
        plt.plot(LarsData[:, 1], LarsData[:, 2], 'ro', label='lada V')
        plt.plot(Vplot[0:-1], y[0:-1], 'k-', label='my V')
        plt.ylabel('y [m]')
        plt.xlabel('velocity [m/s]')
        plt.subplot(233)
        plt.plot(LarsData[:, 3], LarsData[:, 2], 'ro', label='lada k')
        plt.plot(kplot[0:-1], y[0:-1], 'k-', label='my k')
        plt.ylabel('y [m]')
        plt.xlabel('Turbulent kinetic energy [m2/s2]')
        plt.subplot(234)
        plt.plot(LarsData[:, 4], LarsData[:, 2], 'ro', label='lada omega')
        plt.plot(omegaplot[0:-1], y[0:-1], 'k-', label='my omega')
        plt.ylabel('y [m]')
        plt.xlabel('Turbulent dissipation rate [1/s]')
        plt.subplot(235)
        plt.plot(LarsData[:, 5], LarsData[:, 2], 'ro', label='lada mut')
        plt.plot(omegaplot[0:-1], y[0:-1], 'k-', label='my mut')
        plt.ylabel('y [m]')
        plt.xlabel('Turbulent viscosity [pa.s]')
        plt.subplot(236)
        plt.semilogy(iter[1:-1], resU[1:-1], 'blue', label="U-velocity")
        plt.semilogy(iter[1:-1], resV[1:-1], 'black', label="V-velocity")
        plt.semilogy(iter[1:-1], resK[1:-1], 'red', label="k")
        plt.semilogy(iter[1:-1], resomega[1:-1], 'green', label="omega")
        plt.xlabel('Number of iterations')
        plt.ylabel('Residual')
        plt.grid(True)
        plt.grid(which='both')
        plt.show()






