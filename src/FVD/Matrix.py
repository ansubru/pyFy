#!/usr/bin/env python

from pprint import pprint

class A:
    def __init__(self):
        print "From init"


    def f(self,mat):
        tmpMatrix = mat
        for i,x in enumerate(mat):
            for j,y in enumerate(x):
                tmpMatrix[i][j] = y+2
        return tmpMatrix


a = A()

extMat = [[1,2],[3,4]]

pprint(extMat)
extMat = a.f(extMat)
pprint(extMat)




