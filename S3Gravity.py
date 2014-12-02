from numpy import *
from visual import *
from math import *
from random import *

#Gravity simulation in S3

G=0.1
dth = 0.1
N = 10

#Will represent S3 with unit vectors in R4

#Given two points in S3, return the distance between them
def distance(v1, v2):
    if (dot(v1, v2) >= 0.9):
        return 0.0
    return acos(dot(v1, v2))

def magR4(v1):
    return sqrt(dot(v1, v1))

#Returns the norm of a vector in R4
def normR4(v1):
    if (magR4(v1) == 0):
        return vectorR4(1.0, 0.0, 0.0, 0.0)
    return v1 * (1.0 / magR4(v1))

#Given two points in S3, return a tangent vector at v1 pointing toward v2
def tangentdirection(v1, v2):
    #We will use the Gram-Schmidt process
    #to find a vector orthogonal to v1 in the plane of v1 and v2
    #Project v2 onto v1
    v1norm = normR4(v1)
    projv1v2 = (dot(v1norm, v2)) * v1norm
    return normR4(v2 - projv1v2)

#Given a vector on the sphere, a tangent vector, and an angular displacement
#Move the vector in the direction of the tangent vector by the displacement
def translate(r, v, d):
    return (normR4(r*cos(d) + v*sin(d)), normR4(v*cos(d)-r*sin(d)))

#Gives the force of gravity
def gravityforce(m1, m2, d):
    if (d < dth):
        d = dth
    if (d > 3.1415 - dth):
        d = dth
    return G*m1*m2*(1.0/(sin(d) * sin(d)))

#Given a unit vector, a magnitude, and a vector to add, return unit, magnitude of result
def addtounitvec(v1, m1, v2):
    result = (v1 * m1) + v2
    return (normR4(result), magR4(result))

#Projects a vector from R4 onto R3 by dropping last component
def projectR3(v1):
    return vector(v1[0], v1[1], v1[2])

#Gives a fourvector with the given components
def vectorR4(x, y, z, w):
    return array([x, y, z, w])


#Assumption for now: The sphere is of radius 1, and all velocity vectors are normed

class SphericalPointMass:
    def __init__(self, p, mass, velocitydir=vectorR4(1,0,0,0), speed=0.0):
        rand_color = (uniform(0.0, 1.0), uniform(0.0, 1.0), uniform(0.0, 1.0))
        self.sprite = sphere(pos=projectR3(p), radius=0.05, color=rand_color, make_trail=True, retain=200)
        self.mass = mass
        self.velocitydir = velocitydir
        self.speed=speed
        self.pos = p
    #Move along the sphere according to our velocity and a timestep
    def move(self, t):
        r = self.pos
        (rp, vp) = translate(r, self.velocitydir, self.speed * t)
        self.pos = rp
        self.sprite.pos=projectR3(rp)
        self.velocitydir=vp
    #Given the self and a list of others, compute gravitational attraction
    def gravitate(self, others, t):
        accum = vectorR4(0.0,0.0,0.0,0.0)
        for other in others:
            if (not array_equal(other.pos, self.pos)):
                #Compute tangent direction to other
                dir = tangentdirection(self.pos, other.pos)
                #Compute distance
                dist = distance(self.pos, other.pos)
                #Add to the accumulator
                accum += dir * gravityforce(self.mass, other.mass, dist) * t
        (self.velocitydir, self.speed) = addtounitvec(self.velocitydir, self.speed, accum)                  

points = []

seed()
for i in range(N):
    pos = vectorR4(uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0))
    pos = normR4(pos)
    vel = vectorR4(uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0))
    vel = normR4(vel)
    points.append(SphericalPointMass(pos, uniform(0.1, 2.0), vel, 0.0))
    



GlobalSphere = sphere(pos=(0,0,0), radius=1, color=color.blue, opacity=0.5)

while True:
    dt = 1.0/480.0
    rate(1.0/dt)
    for point in points:
        point.gravitate(points, dt)
    for point in points:
        point.move(dt)
        




