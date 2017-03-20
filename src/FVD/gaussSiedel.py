#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = "Ananda S. Kannan"
__copyright__ = "Copyright 2017, pyFy project"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ananda S. Kannan"
__email__ = "ansubru@gmail.com"
__credits to__ Worasait Suwannik http://bit.ly/wannik

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

class gaussSiedel(object):

    def __init__(self):
        """Return object"""

    def gauss_seidel(m, x0=None, eps=1e-5, max_iteration=100):
        """
          Parameters
          ----------
          m  : list of list of floats : coefficient matrix
          x0 : list of floats : initial guess
          eps: float : error tolerance
          max_iteration: int

          Returns
          -------
          list of floats
              solution to the system of linear equation

          Raises
          ------
          ValueError
              Solution does not converge
          """
        n = len(m)
        x0 = [0] * n if x0 == None else x0
        x1 = x0[:]

        for __ in range(max_iteration):
            for i in range(n):
              s = sum(-m[i][j] * x1[j] for j in range(n) if i != j)
              x1[i] = (m[i][n] + s) / m[i][i]
            if all(abs(x1[i]-x0[i]) < eps for i in range(n)):
              return x1
            x0 = x1[:]
        raise ValueError('Solution does not converge')




