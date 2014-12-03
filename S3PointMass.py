from PointMass import *
from Common import *

class S3PointMass(PointMass):
        #Will represent S3 with unit vectors in R4

    #Given two points in S3, return the distance between them
    def distance(self, v1, v2):
        if (dot(v1, v2) >= 0.9):
            return 0.0
        return acos(dot(v1, v2))

    #Given two points in S3, return a tangent vector at v1 pointing toward v2
    def tangentdirection(self, v1, v2):
        #We will use the Gram-Schmidt process
        #to find a vector orthogonal to v1 in the plane of v1 and v2
        #Project v2 onto v1
        v1norm = normR4(v1)
        projv1v2 = (dot(v1norm, v2)) * v1norm
        return normR4(v2 - projv1v2)

    #Given a vector on the sphere, a tangent vector, and an angular displacement
    #Move the vector in the direction of the tangent vector by the displacement
    def translate(self, r, v, d):
        return (normR4(r*cos(d) + v*sin(d)), normR4(v*cos(d)-r*sin(d)))

    #Gives the force of gravity
    def gravityforce(self, m1, m2, d):
        if (d < dth):
            d = dth
        if (d > 3.1415 - dth):
            d = dth
        return G*m1*m2*(1.0/(sin(d) * sin(d)))
    #Given a unit vector, a magnitude, and a vector to add, return unit, magnitude of result
    def addtounitvec(self, v1, m1, v2):
        vec1 = v1 * m1
        result = vec1 + v2
        return (normR4(result), magR4(result))
    #Return a random point
    @staticmethod
    def randompoint():
        return randomnormedR4()
    #Return a random tangent vector
    @staticmethod
    def randomvel():
        return randomnormedR4()
    def updateappearance(self):
        self.sprite.pos = projectR3(self.pos)
    def zerovec(self):
        return vectorR4(0.0,0.0,0.0,0.0)
