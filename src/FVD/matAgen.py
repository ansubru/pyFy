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
################################################################################MODULE FOR CREATING 'A' MATRIX################################################################################################
import os
import numpy as np
import re
import sys
import IO
import Interpolate
import pprint
#####################################################################################################################################
class matAgen(object):

    """A class for creating 'A' matrix needed for solving equations

    Attributes:
       U --> Reads in a given U vector
       Returns A matrix after applying Discretization to it
    """

    def __init__(self):
        """Return object"""
 ###----------------------------------------------------------------------------RELEVANT FUNCTION DEFINITIONS-----------------------------------------------------------------------------------------------###

    #############################################################################################################################################################################################################
    ###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
    def matA(self, matU,matP,flag):
        """A function that returns A matrix after discretization using first order upwind"""
        f = flag
        from IO import IO
        IO_obj = IO("random")
        grid = IO_obj.grid_nodes #generate a base grid node layout

        #Boundary conditions
        if f in ('U', 'u'):
            UA = IO_obj.UA
            UB = IO_obj.UB
            UC = IO_obj.UC
            UD = IO_obj.UD
        else:
            UA = IO_obj.VA
            UB = IO_obj.VB
            UC = IO_obj.VC
            UD = IO_obj.VD

        #Initialize all variables

        U = matU
        P = matP
        i = np.size(U,0) #get indices for rows
        j = np.size(U,1) #get indices for columns
        k = i*j #total number of grid nodes

        from Discretize import Discretize
        disc_obj = Discretize()


        # Discretize U velocity using FOU
        Px, dPx, dPy, Fe, Fw, Fn, Fs, ufe, ufw, ufn, ufs, aP, aW, aE, aN, aS, SUx, SPx, SUy, SPy = disc_obj.FOU_disc(U, P, UA, UB, UC, UD)
