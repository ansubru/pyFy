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
        """Function that generates a x*y matrix based on grid I/P"""
        a += 2
        b += 2
        grid_nodes = np.zeros(shape=(a,b))
        return grid_nodes

#Open pyFy/usr/setup.txt
    dataStorage = {}
    with open(os.path.join(pathf2, "setup.json"),'r') as json_data:
        dataStorage = json.load(json_data)
        data = dataStorage
        x = 1.0*data["xNodes"]
        y = 1.0*data["yNodes"]
        domainSize = data["domain"]
        nu = data["viscosity"]
        grid_zeros = np.array(grid_gen(x, y))
        grid_x = np.array(grid_gen((2.0*x-1.0), (2.0*y-1.0)))
        grid_y = np.array(grid_gen((2.0*x-1.0), (2.0*y-1.0)))
        x_dis = (domainSize[0]*1.0)/x
        y_dis = (domainSize[1]*1.0)/y
        delta_x = x_dis/2.0
        delta_y = y_dis/2.0

#Create array with x,y co-ordinates from generated grid
        i = np.size(grid_x,0)
        j = np.size(grid_y,1)

        for m in range(0,j):
            for n in range(0,i):
                if(m==n==0):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = y
                elif(m == 0 and n < i):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = y
                elif(m > 0 and n > 0):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = (y - ((m)*delta_y))
                elif(n == 0 and m == j):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = 0.0
                elif(n == 0 and m < j):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = (y - ((m)*delta_y))
                elif(n == i and m == j):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = 0.0
                elif(n == i and m == 0):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = y

#############################################################################################################################################################################################################

