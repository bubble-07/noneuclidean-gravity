from PointMass import *
from H3PointMass import *
from S3PointMass import *
import wx

N = 40

points = []

seed()
for i in range(N):
    pos = S3PointMass.randompoint()
    vel = S3PointMass.randomvel()
    points.append(S3PointMass(pos, uniform(0.1, 2.0), vel, 0.0))


GlobalSphere = sphere(pos=(0,0,0), radius=1, color=color.blue, opacity=0.5)

while True:
    dt = 1.0/480.0
    rate(1.0/dt)
    for point in points:
        point.gravitate(points, dt)
    for point in points:
        point.move(dt)
        


