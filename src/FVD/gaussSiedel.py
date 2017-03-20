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

#############################################################################################################################################################################################################
###----------------------------------------------------------------------------CLASS DEFINITION-----------------------------------------------------------------------------------------------###

class gaussSiedel(object):

    def __init__(self):
        """Return object"""

    def gauss_seidel(self, a, b, tol_typ, tol):
        """#Takes an nxn matrix and a 1xn vector and solves for x by iterating until the
            given tolerance (tol) is met. tol_type is a single character to indicate absolute
            or relative convergence ('a' or 'r').Prints the solution and required iterations
            to meet the tolerance."""
        print("Trying to solve equations using the Gauss-Seidel method")
        print ("A matrix",a)
        print ("b matrix",b)
        shape = np.shape(a)
        tol_type = tol_typ
        m = shape[0]
        n = shape[1]

        if m != n:
            print("This solver only works for square matrices")
            print("This matrix is %dx%d." % (m, n))
            exit(1)

        if m != np.shape(b)[0]:
            print("b must be the same dimensions as A.")
            print("b appears to be %d elements long" % np.shape(b)[0])
            exit(1)

        x = np.zeros(np.shape(b))
        sum = np.zeros(np.shape(b))
        prev_x = np.zeros(np.shape(b))
        diff = np.zeros(np.shape(b))

        if tol_type in "a":
            numpy_solution = np.linalg.solve(a, b);

        num_iterations = 0
        error = tol + 1 # Mock up error parameter only to run the While
        while error > tol:
            for i in range(m):
                prev_x[i] = x[i]
                sum[i] = b[i]
                for j in range(n):
                    if i != j:
                        sum[i] = sum[i] - a[i][j] * x[j]
                sum[i] = sum[i] / a[i][i]
                x[i] = sum[i]
            num_iterations += 1
            if tol_type in "a":
                # print("Absolute tolerance")
                diff = np.subtract(x, numpy_solution)
                error = np.linalg.norm(diff) / np.linalg.norm(x)
            if tol_type in "r":
                diff = np.subtract(x, prev_x)
                error = np.linalg.norm(diff) / np.linalg.norm(x)

        if tol_type in "a":
            print("Using GaussSeidel to converge to an absolute error of %.8f requiring %d iterations." % (
            tol, num_iterations))
            print("The solution is:")
            print(x)
        if tol_type in "r":
            print(
            "Using GaussSeidel to converge to a relative error of %.8f requiring %d iterations." % (tol, num_iterations))
            print("The solution is:")
            print(x)
        return x












