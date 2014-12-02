from numpy import *
from visual import *
from math import *
from random import *
from Common import *

G = 0.3
dth = 0.1


class PointMass:
    def __init__(self, p, mass, velocitydir, speed):
        rand_color = (uniform(0.0, 1.0), uniform(0.0, 1.0), uniform(0.0, 1.0))
        self.sprite = sphere(pos=(0,0,0), radius=0.05, color=rand_color, make_trail=True, retain=200)
        self.mass = mass
        self.velocitydir = velocitydir
        self.speed=speed
        self.pos = p
   
    #Move in the space according to the velocity and a timestep
    def move(self, t):
        r = self.pos
        (rp, vp) = self.translate(r, self.velocitydir, self.speed * t)
        self.pos = rp
        self.velocitydir=vp
        self.updateappearance() #TODO: Implement this!

    #TODO: should implement zero vector as a thing
    def gravitate(self, others, t):
        accum = self.zerovec()
        for other in others:
            if (not array_equal(other.pos, self.pos)):
                #Compute tangent direction to other
                dir = self.tangentdirection(self.pos, other.pos)
                #Compute distance
                dist = self.distance(self.pos, other.pos)
                #Add to the accumulator
                accum += dir * self.gravityforce(self.mass, other.mass, dist) * t
        (self.velocitydir, self.speed) = self.addtounitvec(self.velocitydir, self.speed, accum)                  
