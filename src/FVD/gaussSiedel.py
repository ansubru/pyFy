#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = "Ananda S. Kannan"
__copyright__ = "Copyright 2017, pyFy project"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ananda S. Kannan"
__email__ = "ansubru@gmail.com"
__credits to__ http://austingwalters.com/gauss-seidel-method/

"""
################################################################################MODULE FOR GAUSS SEIDEL METHOD################################################################################################

import os
import numpy as np
import re
import sys
import pprint

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###
from scipy.linalg import solve

class gaussSiedel(object):

    """A class for Gauss seidel method

    Attributes:
       A matrix
       b matrix --> solution to Ax = b
       x matrix --> values of x
       number of iterations n
    """

    def __init__(self):
        """Return object"""



    def gauss(self, A, b, x, n):

        L = np.tril(A)
        U = A - L
        for i in range(n):
            x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
            print str(i).zfill(3),
            print(x)
        return x

