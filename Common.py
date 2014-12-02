from numpy import *
from visual import *
from math import *
from random import *

#Gives a fourvector with the given components
def vectorR4(x, y, z, w):
    return array([x, y, z, w])
def magR4(v1):
    return sqrt(dot(v1, v1))

#Returns the norm of a vector in R4
def normR4(v1):
    if (magR4(v1) == 0):
        return vectorR4(1.0, 0.0, 0.0, 0.0)
    return v1 * (1.0 / magR4(v1))

#Given two orthonormal vectors and an angle, find the vector obtained by
#rotating the first by the angle through the plane defined by the two.
def orthonorm_rotate(h, v, th):
    return cos(th) * h + sin(th) * v
#Projects a vector from R4 onto R3 by dropping last component
def projectR3(v1):
    return vector(v1[0], v1[1], v1[2])


