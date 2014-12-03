from PointMass import *
from random import *

class H3PointMass(PointMass):

    #Gravity simulation in poincare sphere

    #Given a point and a vector representing a heading
    #Determine the center of the circle representing the geodesic
    #That the given point would travel along
    def get_geodesic_center(self, x, v):
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
    def get_geodesic_radius(self, c):
        return sqrt(mag2(c) - 1.0)

    #Given masses and a distance, determine gravitational force
    def gravityforce(self, m1, m2, d):
        if (d < dth):
            d = dth
        return G*m1*m2*(1.0/(sinh(d)*sinh(d)))

    #Given masses and a distance, determine gravitational PE
    def gravityenergy(self, m1, m2, d):
        if (d < dth):
            return self.gravityforce(m1, m2, d) * d
        return (-1*G*m1*m2)*((1/tanh(d)) - (1/tanh(dth))) + self.gravityforce(m1, m2, d) * d

    #Given two points in the poincare sphere, return their
    #distance in hyperbolic space
    def distance(self, x1, x2):
        return arccosh(1.0+2.0*(mag2(x2-x1)/((1.0-mag2(x1))*(1.0-mag2(x2)))))

    #Given a position in hyperbolic space and the center of the circle
    #representing a geodesic, return the maximum angle along the circle
    #that the point may be rotated before going to infinity
    def get_max_rotation(self, x, c):
        r = self.get_geodesic_radius(c)
        phi = arccos(1.0/mag(c))
        alpha = arccos((r*r+mag2(c)-mag2(x))/(2*r*mag(c)))
        return 0.5 * pi - phi - alpha

    #Given a radius, give the scaling factor for distances and lengths
    #for all points of that distance from the origin
    def get_scale_factor(self, r):
        return 1 - r*r


    #Given two points in the disk, compute a unit tangent vector from x1 to x2 at x1
    def tangentdirection(self, x1, x2):
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
    def translate_euclidean_approx(self, x, v, d):
        return ((d * self.get_scale_factor(mag(x))) * v + x, v)

    def addtounitvec(self, v1, m1, v2):
        vec1 = v1 * m1
        result = vec1 + v2
        return (norm(result), mag(result))

    #TODO: Handle case where we have a diameter
    #Given a point, a heading, and a displacement, move the point
    #Return the point and a new heading
    def translate(self, x, v, d):
        if (d == 0):
            return x
        #Determine the center of the sphere the geodesic lies on
        c = self.get_geodesic_center(x, v)
        
        #Determine its radius
        r = mag(x - c)
        if (r > self.max_radius):
            return self.translate_euclidean_approx(x, v, d)
        #FIXME: Calculation for r is wrong if use get_geodesic_radius
        
        #Get a normal to the plane our geodesic lies in
        n = norm(cross(x, v))
        #Get an orthonormal basis within the plane
        h = norm(x - c)
        vh = norm(cross(h, n))
        
        #Use binary search to determine how far is far enough.
        th_min = 0.0
        th_max = self.get_max_rotation(x, c)
        th_current = (th_min + th_max) / 2.0
        d_est = 0.0 #Estimated distance
        x2 = x #Resultant vector
        while (abs(d_est - d) > abs(d) * self.precision):
            #TODO: handle case that what we're asking for is too far away!
            x2 = c + r * orthonorm_rotate(h, vh, th_current)

            d_est = self.distance(x, x2)
                
            if (d_est > d):
                #We are too far away. Make the current theta the max
                th_max = th_current
            else:
                #We are too close. Make the current theta the min
                th_min = th_current
            th_current = (th_min + th_max) / 2.0 #Try again at the middle of the interval
        #Force the magnitude of x2 to not meet or exceed 1
        if (mag(x2) >= self.max_position):
            x2 = (2 * self.max_position - 1.0) * norm(x2)
        return (x2, orthonorm_rotate(vh, -h, th_current))
    @staticmethod
    def randompoint():
        return randomnormedR3() * uniform(0.0, 0.8)
    @staticmethod
    def randomvel():
        return randomnormedR3()
    def zerovec(self):
        return vector(0.0,0.0,0.0)
    def updateappearance(self):
        self.sprite.pos = self.pos
        self.sprite.radius = 0.05 * self.get_scale_factor(mag(self.pos))
