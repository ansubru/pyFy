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
################################################################################MODULE FOR READING SETUP FILE#################################################################################################

import os
import numpy as np
import re
import json
from pprint import pprint

##############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
class IO(object):
    """A class for IO of the usr specs

    Attributes:
    None

    """

    def __init__(self,random):
        """Random instantiation of class"""
        self.random = None

#Path for the current file
    temp = "/usr"
    pathf = os.getcwd()
    pathf2 =os.path.dirname(os.path.dirname(pathf)) + temp

    def grid_gen(a, b):
        """Function that generates a x*y matrix (incl boundaries) based on grid I/P"""
        a += 2
        b += 2
        grid = np.zeros(shape=(a,b))
        return grid

    def cloneRowMat(F, G, l, p):
        """A function to clone pth row from a matrix G to into the lth row of a matrix F (p > l)"""
        j = np.size(F, 1)
        for m in range(p):  # loop through rows
            for n in range(j):  # loop through columns
                if (m == p - 1):
                    F[m][n] = G[l - 1][n]
        return F

#Open pyFy/usr/setup.txt
    dataStorage = {}
    with open(os.path.join(pathf2, "setup.json"),'r') as json_data:
        dataStorage = json.load(json_data)
        data = dataStorage
#Data extraction
        x = data["xNodes"]
        y = data["yNodes"]
        domainSize = data["domain"]
        Ubc = data["boundary-velocity"] #Boundary velocity at I/P
        nu = 1.0*data["viscosity"]
        rho = 1.0*data["rho"]
        alpha = 1.0*data["Under-relaxation"]
        x_dis = (domainSize[0]*1.0)/x #x-grid spacing
        y_dis = (domainSize[1]*1.0)/y #y-grid spacing
        growthRatio = 1.0 * data["growthRatio"] #Cell growth ratio (towards the wall)
        U_wall = 0.0 #No slip wall velocity is (0,0)
        Ubcx = 1.0*Ubc[0] #Boundary velocity x-comp
        Ubcy = 1.0*Ubc[1] #Boundary velocity y-comp

        dx = np.array(grid_gen(x, y))  # cell length matrix
        dy = np.array(grid_gen(x, y))  # cell width matrix
        xpt = np.array(grid_gen(x, y))  # x distances
        ypt = np.array(grid_gen(x, y))  # y distances
        delXE = np.array(grid_gen(x, y))
        delXW = np.array(grid_gen(x, y))
        delYN = np.array(grid_gen(x, y))
        delYS = np.array(grid_gen(x, y))
        grid_nodes = grid_gen(x, y)
        grid_x = grid_gen(x, y)
        grid_y = grid_gen(x, y)


        # Create delX and delY grids
        i = np.size(grid_nodes, 0)
        j = np.size(grid_nodes, 1)

        if data["gridType"] in ['Equidistant', 'equidistant']:

            for m in range(0, i):
                for n in range(0, j):

                    dx[m][n] = (domainSize[0]*1.0)/x  #cell size (x-direction)
                    dy[m][n] = (domainSize[1]*1.0)/y  #cell size (y-direction)

            for m in range(0, i):
                for n in range(0, j):

                    # Apply bcs

                    if (m >= 0 and n == 0):  # A Boundary
                        delXE[m][n] = 1.0  # pad bc distance
                        delXW[m][n] = 1.0
                    elif (m == 0 and n != j - 1 and n != 0):  # B Boundary
                        delXE[m][n] = 1.0  # pad bc distance
                        delXW[m][n] = 1.0
                    elif (m >= 0 and n == j - 1):  # C Boundary
                        delXE[m][n] = 1.0  # pad bc distance
                        delXW[m][n] = 1.0
                    elif (m == i - 1 and n != j - 1 and n != 0):  # D Boundary
                        delXE[m][n] = 1.0  # pad bc distance
                        delXW[m][n] = 1.0
                    elif (m >= 1 and n == 1):
                        delXE[m][n] = dx[m][n]  # distance to East neighbour
                        delXW[m][n] = 0.5* dx[m][n]  # distance to West neighbour
                    elif (m >= 1 and n == i-2):
                        delXW[m][n] = dx[m][n] # distance to West neighbour
                        delXE[m][n] = 0.5 * dx[m][n]  # distance to East neighbour
                    else:
                        delXW[m][n] = dx[m][n]
                        delXE[m][n] = dx[m][n]

            delYS = np.transpose(delXE) # distance to North neighbour
            delYN = np.transpose(delXW) # distance to South neighbour


        elif data["gridType"] in ['Non-Equidistant', 'non-equidistant']:

            for m in range(0, i):
                for n in range(0, j):

                    #Calculate size of wall cell using growth ratio y = a*(1-x^n)/(1-x)
                    nx = (x)/2 ; ny = (y)/2 #number of nodes till half way
                    ax = (0.5*domainSize[0]) / ((1 - growthRatio ** nx) / (1 - growthRatio))
                    ay = (0.5 *domainSize[1]) / ((1 - growthRatio ** ny) / (1 - growthRatio))

                    if m >= 0 and n < j/2:
                        xpt[m][n] = ax * ((1 - growthRatio ** n) / (1 - growthRatio))

                    if m >= 0 and n >= j/2:
                        mirrX = j / 2 - (n - j / 2) - 1
                        xpt[m][n] = 1 - xpt[m][mirrX]
                    if m == n == 0 :
                        xpt[m][n] = 0.0

            ypt = np.transpose(xpt)
            ypt = 1 - ypt

            for m in range(0, i):
                for n in range(0, j):
                    if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes
                        dx[m][n] = xpt[m][n+1] - xpt[m][n]
                        dy[m][n] = ypt[m+1][n] - xpt[m][n]
            print dx[1]
            print xpt[1]

#Boundary conditions
        if data["A-bound"] in ['walls', 'Walls', 'wall', 'Wall']:
            UA = VA = U_wall

        elif data["A-bound"] in ['fixed velocity']:
            UA = Ubcx
            VA = Ubcy

        if data["B-bound"] in ['walls', 'Walls', 'wall', 'Wall']:
            UB = VB = U_wall

        elif data["B-bound"] in ['fixed velocity']:
            UB = Ubcx
            VB = Ubcy

        if data["C-bound"] in ['walls', 'Walls', 'wall', 'Wall']:
            UC = VC =  U_wall

        elif data["C-bound"] in ['fixed velocity']:
            UC = Ubcx
            VC = Ubcy

        if data["D-bound"] in ['walls', 'Walls', 'wall', 'Wall']:
            UD = VD =  U_wall

        elif data["D-bound"] in ['fixed velocity']:
            UD = Ubcx
            VD = Ubcy

        # Create array with x,y co-ordinates from generated grid
        i = np.size(grid_nodes, 0)
        j = np.size(grid_nodes, 1)

        # for m in range(0, i):
        #     for n in range(0, j):
        #
        #         if (m != 0 and n != 0 and m != (i-1) and n != (j-1)):
        #             grid_x[m][n] = (n * dy[m][n] * 0.5) + ((n-1) * dy[m][n] * 0.5)
        #         if ( m == i-1 or n == j-1):
        #             grid_x[m][n] = 1.0
        #
        #
        # grid_x = cloneRowMat(grid_x,grid_x,2,1)
        # grid_x = cloneRowMat(grid_x, grid_x, 11, 12)
        # grid_y = np.transpose((grid_x))

        # for m in range(0, i):
        #     for n in range(0, j):
        #         grid_y[m][n] = 1 - grid_y[m][n]

        for m in range(i):
            for n in range(j):
                if m==0:
                    x = 0
                elif m<i-1:
                    x=dy[m][n]*0.5 + dy[m][n]*(m-1)
                else:
                    x= domainSize[0]
                if n==0:
                    y = 0
                elif n<j-1:
                    y=dy[m][n]*0.5 + dy[m][n]*(n-1)
                else:
                    y= domainSize[1]
                grid_x[n][m] = x
                grid_y[n][m] = 1-y


#############################################################################################################################################################################################################

