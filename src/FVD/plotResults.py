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
from matplotlib.pyplot import draw, figure, show
import pylab



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

    def genGridLam(self, matU):

        from IO import IO
        IO_obj = IO("random")
        dy = IO_obj.dy
        grid_x = 0.0*matU
        grid_y = 0.0 * matU
        # Create array with x,y co-ordinates from generated grid
        i = np.size(matU, 0)
        j = np.size(matU, 1)

        for m in range(i):
            for n in range(j):
                if m == 0:
                    x = 0
                elif m < i - 1:
                    x = dy[m][n] * 0.5 + dy[m][n] * (m - 1)
                else:
                    x = 0
                if n == 0:
                    y = 0
                elif n < j - 1:
                    y = dy[m][n] * 0.5 + dy[m][n] * (n - 1)
                else:
                    y = 1
                grid_x[n][m] = x
                grid_y[n][m] = 1 - y
        return grid_x, grid_y

    def plotdataTurb(self, U, V, k, omega, Y, mut, resU, resV, resK, resomega, resB, resPP, it):

        i = np.size(U, 0)
        j = np.size(V, 1)

        y = []
        Uplot = []
        Vplot = []
        kplot = []
        omegaplot = []
        mutplot = []
        iter = list(range(1, it + 2))

        UToPlot = np.zeros([i, j])
        VToPlot = np.zeros([i, j])
        kToPlot = np.zeros([i, j])
        omegaToPlot = np.zeros([i, j])
        mutToPlot = np.zeros([i, j])

        from IO import IO
        IO_obj = IO("random")
        fxe = IO_obj.fxe

        from Interpolate import Interpolate
        Interp_obj = Interpolate()
        for m in range(i):
            for n in range(j):
                UToPlot[m][n] = U[-m-1][n]
                VToPlot[m][n] = V[-m - 1][n]
                kToPlot[m][n] = k[-m - 1][n]
                omegaToPlot[m][n] = omega[-m - 1][n]
                mutToPlot[m][n] = mut[-m - 1][n]


        for m in range(i):
            for n in range(j):
                if m>0 and n == 31:
                    y.append(Y[m, n])
                    Uplot.append(Interp_obj.weighted_interp(UToPlot[m, n], UToPlot[m, n + 1], fxe[m][n]))
                    Vplot.append(Interp_obj.weighted_interp(VToPlot[m, n], VToPlot[m, n + 1], fxe[m][n]))
                    kplot.append(Interp_obj.weighted_interp(kToPlot[m, n], kToPlot[m, n + 1], fxe[m][n]))
                    mutplot.append(Interp_obj.weighted_interp(mutToPlot[m, n], mutToPlot[m, n + 1], fxe[m][n]))
                    omegaplot.append(Interp_obj.weighted_interp(omegaToPlot[m, n], omegaToPlot[m, n + 1], fxe[m][n]))

        LarsData = np.loadtxt('LarsDataTurb.txt', skiprows=1)

        pylab.figure(figsize=(30, 10))
        pylab.subplot(231)
        pylab.plot(LarsData[:, 0], LarsData[:, 2], 'ro', label='lada U')
        pylab.plot(Uplot[0:-1], y[0:-1], 'k.', label='my U')
        pylab.ylabel('y [m]')
        pylab.xlabel('velocity [m/s]')
        pylab.legend(loc='upper right')
        pylab.subplot(232)
        pylab.plot(LarsData[:, 1], LarsData[:, 2], 'ro', label='lada V')
        pylab.plot(Vplot[0:-1], y[0:-1], 'k.', label='my V')
        pylab.ylabel('y [m]')
        pylab.xlabel('velocity [m/s]')
        pylab.legend(loc='upper right')
        pylab.subplot(233)
        pylab.plot(LarsData[:, 3], LarsData[:, 2], 'ro', label='lada k')
        pylab.plot(kplot[0:-1], y[0:-1], 'k.', label='my k')
        pylab.ylabel('y [m]')
        pylab.xlabel('Turbulent kinetic energy [m2/s2]')
        pylab.legend(loc='upper right')
        pylab.subplot(234)
        pylab.plot(LarsData[:, 4], LarsData[:, 2], 'ro', label='lada omega')
        pylab.plot(omegaplot[0:-1], y[0:-1], 'k.', label='my omega')
        pylab.ylabel('y [m]')
        pylab.xlabel('Turbulent dissipation rate [1/s]')
        pylab.legend(loc='upper right')
        pylab.subplot(235)
        pylab.plot(LarsData[:, 5], LarsData[:, 2], 'ro', label='lada mut')
        pylab.plot(mutplot[0:-1], y[0:-1], 'k.', label='my mut')
        pylab.ylabel('y [m]')
        pylab.xlabel('Turbulent viscosity [pa.s]')
        pylab.legend(loc='upper right')
        pylab.subplot(236)
        pylab.semilogy(iter[1:-1], resU[1:-1], 'blue', label='U-velocity')
        pylab.semilogy(iter[1:-1], resV[1:-1], 'black', label='V-velocity')
        pylab.semilogy(iter[1:-1], resPP[1:-1], 'magenta', label='Pprime')
        pylab.semilogy(iter[1:-1], resK[1:-1], 'red', label='k')
        pylab.semilogy(iter[1:-1], resomega[1:-1], 'green', label='omega')
        pylab.semilogy(iter[1:-1], resB[1:-1], 'cyan', label='continuity')
        pylab.xlabel('Number of iterations')
        pylab.ylabel('Residual')
        pylab.grid(True)
        pylab.grid(which='both')
        pylab.show()

        pylab.figure
        pylab.plot(LarsData[:, 4], LarsData[:, 2], 'ro', label='lada omega')
        pylab.plot(omegaplot[0:-1], y[0:-1], 'k.', label='my omega')
        pylab.ylabel('y [m]')
        pylab.xlabel('Turbulent dissipation rate [1/s]')
        pylab.xlim(0, 10)
        pylab.legend(loc='upper right')
        pylab.show()

    def plotdata2(self, U, V, Y, resU, resV, resB, it):

        i = np.size(U, 0)
        j = np.size(V, 1)

        UToPlot = np.zeros([i, j])
        VToPlot = np.zeros([i, j])
        YToPlot = np.zeros([i, j])

        y = []
        Uplot = []
        Vplot = []
        iter = list(range(1, it + 1))
        from IO import IO
        IO_obj = IO("random")
        fxe = IO_obj.fxe

        from Interpolate import Interpolate
        Interp_obj = Interpolate()

        for m in range(i):
            for n in range(j):
                UToPlot[m][n] = U[-m-1][n]
                VToPlot[m][n] = V[-m - 1][n]
                YToPlot[m][n] = Y[-m - 1][n]

        for m in range(i):
            for n in range(j):
                if m > 0 and n == 5:
                    y.append(Y[m, n])
                    Uplot.append(Interp_obj.weighted_interp(U[m, n], U[m, n + 1], fxe[m][n]))
                    Vplot.append(Interp_obj.weighted_interp(V[m, n], V[m, n + 1], fxe[m][n]))

        LarsData = np.loadtxt('LarsData.txt', skiprows=1)

        plt.figure
        plt.subplot(131)
        plt.plot(LarsData[:, 0], LarsData[:, 2], 'ro', label='lada U')
        plt.plot(Uplot[0:-1], y[0:-1], 'k-', label='my U')
        plt.ylabel('y [m]')
        plt.xlabel('velocity [m/s]')
        plt.legend(loc='upper right')
        plt.subplot(132)
        plt.plot(LarsData[:, 1], LarsData[:, 2], 'ro', label='lada V')
        plt.plot(Vplot[0:-1], y[0:-1], 'k-', label='my V')
        plt.ylabel('y [m]')
        plt.xlabel('velocity [m/s]')
        plt.legend(loc='upper right')
        plt.subplot(133)
        plt.semilogy(iter[0:-1], resU[1:-1], 'blue', label="U-velocity")
        plt.semilogy(iter[0:-1], resV[1:-1], 'black', label="V-velocity")
        plt.semilogy(iter[0:-1], resB[1:-1], 'cyan', label="continuity")
        plt.xlabel('Number of iterations')
        plt.ylabel('Residual')
        plt.grid(True)
        plt.grid(which='both')
        plt.legend(loc='upper right')
        plt.show()






