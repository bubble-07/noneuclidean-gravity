from numpy import *
from visual import *
from math import *
from random import *

#Gravity simulation in poincare sphere

#Convenience func to convert numpy column vec to vector
def make_vector(x):
    return vector(x[0][0], x[1][0], x[2][0])

#Given a point and a vector representing a heading
#Determine the center of the circle representing the geodesic
#That the given point would travel along
def get_geodesic_center(x, v):
    #If x and v are in the same direction
    if (mag(cross(x, v)) == 0.0):
        #Center is a point at infinity (on diameter)
        inf = float('inf')
        return vector(inf, inf, inf)
        
    A = array([x.astuple(), v.astuple(), cross(x, v).astuple()])
    B = array([[mag2(x)+1], [dot(x, v)], [0]])
    return make_vector(linalg.inv(A).dot(B))

#Given the center of a sphere that a geodesic lies on
#in the poincare model, find its radius
def get_geodesic_radius(c):
    return sqrt(mag2(c) - 1.0)

#Number of particles
N = 20

#Universal gravitational constant
G = 0.9

#Minimum threshold distance beyond which gravitation becomes eased
dth = 0.1

#Given masses and a distance, determine gravitational force
def gravityforce(m1, m2, d):
    if (d < dth):
        d = dth
    return G*m1*m2*(1.0/(sinh(d)*sinh(d)))

#Given masses and a distance, determine gravitational PE
def gravityenergy(m1, m2, d):
    if (d < dth):
        return gravityforce(m1, m2, d) * d
    return (-1*G*m1*m2)*((1/tanh(d)) - (1/tanh(dth))) + gravityforce(m1, m2, d) * d

#Given two points in the poincare sphere, return their
#distance in hyperbolic space
def get_hyperbolic_distance(x1, x2):
    return arccosh(1.0+2.0*(mag2(x2-x1)/((1.0-mag2(x1))*(1.0-mag2(x2)))))

#Given two orthonormal vectors and an angle, find the vector obtained by
#rotating the first by the angle through the plane defined by the two.
def orthonorm_rotate(h, v, th):
    return cos(th) * h + sin(th) * v

#Given a position in hyperbolic space and the center of the circle
#representing a geodesic, return the maximum angle along the circle
#that the point may be rotated before going to infinity
def get_max_rotation(x, c):
    r = get_geodesic_radius(c)
    phi = arccos(1.0/mag(c))
    alpha = arccos((r*r+mag2(c)-mag2(x))/(2*r*mag(c)))
    return 0.5 * pi - phi - alpha

#Given a radius, give the scaling factor for distances and lengths
#for all points of that distance from the origin
def get_scale_factor(r):
    return 1 - r*r


#Given two points in the disk, compute a unit tangent vector from x1 to x2 at x1
def tangentdirection(x1, x2):
    #Compute a normal to the plane
    n = cross(x1, x2)
    #if both are along some diameter
    if (mag(n) == 0.0):
        #Return a direction along the diameter
        return norm(x2 - x1)
    #otherwise, we should be able to solve for the center
    A = array([x1.astuple(), x2.astuple(), n.astuple()])
    B = array([[mag2(x1)+1], [mag2(x2)+1], [0]])
    #Center of the geodesic
    c = make_vector(linalg.inv(A).dot(B))
    
    #Tangent vector at x1
    v = norm(cross((c-x1), n))
    #Ensure that v is pointing in the right direction
    if (dot(v, x2-x1) < 0.0):
        return -v
    return v

#Precision for translations (to be proportional to distances)
precision = 0.01

#Maximum radius for a geodesic before we use a quasi-Euclidean approx
max_radius = 100000

#Maximum magnitude of a position before we renormalize it
max_position = 0.999999999

#Helper method that translates using a simple, locally-Euclidean approx
def translate_euclidean_approx(x, v, d):
    return ((d * get_scale_factor(mag(x))) * v + x, v)

#TODO: Handle case where we have a diameter
#Given a point, a heading, and a displacement, move the point
#Return the point and a new heading
def translate(x, v, d):
    if (d == 0):
        return x
    #Determine the center of the sphere the geodesic lies on
    c = get_geodesic_center(x, v)
    
    #Determine its radius
    r = mag(x - c)
    if (r > max_radius):
        return translate_euclidean_approx(x, v, d)
    #FIXME: Calculation for r is wrong if use get_geodesic_radius
    
    #Get a normal to the plane our geodesic lies in
    n = norm(cross(x, v))
    #Get an orthonormal basis within the plane
    h = norm(x - c)
    vh = norm(cross(h, n))
    
    #Use binary search to determine how far is far enough.
    th_min = 0.0
    th_max = get_max_rotation(x, c)
    th_current = (th_min + th_max) / 2.0
    d_est = 0.0 #Estimated distance
    x2 = x #Resultant vector
    while (abs(d_est - d) > abs(d) * precision):
        #TODO: handle case that what we're asking for is too far away!
        x2 = c + r * orthonorm_rotate(h, vh, th_current)

        d_est = get_hyperbolic_distance(x, x2)
            
        if (d_est > d):
            #We are too far away. Make the current theta the max
            th_max = th_current
        else:
            #We are too close. Make the current theta the min
            th_min = th_current
        th_current = (th_min + th_max) / 2.0 #Try again at the middle of the interval
    #Force the magnitude of x2 to not meet or exceed 1
    if (mag(x2) >= max_position):
        x2 = (2 * max_position - 1.0) * norm(x2)
    return (x2, orthonorm_rotate(vh, -h, th_current))


#Given a unit vector, a magnitude, and a vector to add, return unit, magnitude of result
def addtounitvec(v1, m1, v2):
    vec1 = v1 * m1
    result = vec1 + v2
    return (norm(result), mag(result))

E = 0

class HyperbolicPointMass:
    def __init__(self, p, mass, velocitydir=vector(1,0,0), speed=0.0):
        rand_color = (uniform(0.0, 1.0), uniform(0.0, 1.0), uniform(0.0, 1.0))
        self.sprite = sphere(pos=p, radius=0.05, color=rand_color, make_trail=True, retain=200)
        self.mass = mass
        self.velocitydir = velocitydir
        self.speed=speed
    #Move inside the disk according to our velocity and a timestep
    def move(self, t):
        r = self.sprite.pos
        #Get our new direction and heading
        (rp, vp) = translate(r, self.velocitydir, self.speed * t)
        self.sprite.pos=rp
        self.velocitydir=vp
        #Adjust apparent size
        self.sprite.radius = 0.05 * get_scale_factor(mag(self.sprite.pos))
        
    #Given the self and a list of others, compute gravitational attraction
    def gravitate(self, others, t):
        global E
        accum = vector(0,0,0)
        for other in others:
            if (not (other.sprite.pos == self.sprite.pos)):
                #Compute tangent direction to other
                dir = tangentdirection(self.sprite.pos, other.sprite.pos)
                #Compute distance
                dist = get_hyperbolic_distance(self.sprite.pos, other.sprite.pos)
                #Add to the accumulator
                accum += dir * gravityforce(self.mass, other.mass, dist) * t
                E += gravityenergy(self.mass, other.mass, dist)
        (self.velocitydir, self.speed) = addtounitvec(self.velocitydir, self.speed, accum)
        E += (1/2) * self.mass * self.speed * self.speed

def main():

    points = []

    seed()
    for i in range(N):
        pos = vector(uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0))
        pos = norm(pos) * uniform(-0.8, 0.8)
        points.append(HyperbolicPointMass(pos, uniform(0.1, 0.2), vector(0.0,1.0,0.0), 0.0))
        



    GlobalSphere = sphere(pos=(0,0,0), radius=1, color=color.blue, opacity=0.5)

    while (True):
        dt = 1.0/30.0
        rate(1.0/dt)
        global E
        E = 0
        for point in points:
            point.gravitate(points, dt)
        print E
        sum = vector(0.0, 0.0, 0.0)
        for point in points:
            point.move(dt)
            sum += point.sprite.pos
main()
    
    
    

    
    


