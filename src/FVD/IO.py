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
        crv = 1.0 * data["CRV"] #Cell refinement value
        x_dis = (domainSize[0]*1.0)/x #x-grid spacing
        y_dis = (domainSize[1]*1.0)/y #y-grid spacing
        U_wall = 0.0 #No slip wall velocity is (0,0)
        Ubcx = 1.0*Ubc[0] #Boundary velocity x-comp
        Ubcy = 1.0*Ubc[1] #Boundary velocity y-comp

        dx = np.array(grid_gen(x, y))  # cell length matrix
        dy = np.array(grid_gen(x, y))  # cell width matrix
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


        #elif data["gridType"] in ['Non-Equidistant', 'non-equidistant']:

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

        for m in range(0, i):
            for n in range(0, j):
                if(m > 0 and n == 1):
                    grid_x[m][n] = (n)*dy[m][n]*0.5
                if(m > 0 and n == j-2):
                    grid_x[m][n] = ((j-2)*dy[5][5]) - (dy[m][n] * 0.5)
                if (m > 0 and n == j - 1):
                    grid_x[m][n] = ((j - 2) * dy[m][n])
                if (m == 1 and n > 0):
                    grid_y[m][n] = (m) * dx[m][n] * 0.5
                if (m == i-2 and n > 0):
                    grid_y[m][n] = ((i-2) * dx[5][5]) - (dx[m][n] * 0.5)
                if (m == i-1 and n > 0):
                    grid_y[m][n] = ((i-2) * dx[5][5])
                if( m > 0 and n>1 and n < j-2):
                    grid_x[m][n] = n * dy[m][n]
                if (n > 0  and m > 1 and m < i - 2):
                    grid_y[m][n] = m * dx[m][n]





#############################################################################################################################################################################################################

