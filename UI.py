import PointMass
from H3PointMass import *
from S3PointMass import *
import wx

#Create the main window
L = 320
Hregion = 400
w = window(width=2*(L+window.dwidth), height=L+window.dheight+Hregion,
            menus = False, title='Noneuclidean Gravity')
d = 20
disp = display(window=w, x=d, y=d, width = L-2*d, height=L-2*d)

p = w.panel


minspeedrangeselector = wx.TextCtrl(p, pos=(1.1*L, 0.1*L), value='0', size=(50, 20))
maxspeedrangeselector = wx.TextCtrl(p, pos=(1.4*L, 0.1*L), value='0', size=(50, 20))
wx.StaticText(p, pos=(1.2*L, 0.05*L), label='Speed range')

minmassrangeselector = wx.TextCtrl(p, pos=(1.1*L, 0.2*L), value='0', size=(50, 20))
maxmassrangeselector = wx.TextCtrl(p, pos=(1.4*L, 0.2*L), value='0', size=(50, 20))
wx.StaticText(p, pos=(1.2*L, 0.15*L), label='Mass range')

gravselector = wx.TextCtrl(p, pos=(1.25*L, 0.3*L), value='0.3', size=(50,20))
wx.StaticText(p, pos=(1.25*L, 0.25*L), label='G')

numselector = wx.TextCtrl(p, pos=(1.25*L, 0.4*L), value='10', size=(50,20))
wx.StaticText(p, pos=(1.25*L, 0.35*L), label='Particles')

modeselector = wx.RadioBox(p, pos=(1.25*L, 0.5*L), size=(0.25*L, 0.2*L), 
                           choices = ['S3', 'H3'])
minmass = 0.1
maxmass = 2.0
mode = 'S3'
N = 40
PointMass.G = 0.3

#loads a floating point value from a given string
#throws exception on malformed input
#TODO: Throw proper exceptions!
def loadfloat(string):
    return float(string)
def loadint(string):
    return int(string)
def loadmode(string):
    return string

def loadconsts():
    global minmass
    global maxmass
    global PointMass
    global N
    global mode
    
    minmass = loadfloat(minmassrangeselector.GetValue())
    maxmass = loadfloat(minmassrangeselector.GetValue())
    PointMass.G = loadfloat(gravselector.GetValue())
    N = loadint(numselector.GetValue())
    mode = modeselector.GetStringSelection()


points = []
def regenerate(evt):
    loadconsts()

    global points
    for point in points:
        point.cleanup()
    points = []
    for i in range(N):
        pos = S3PointMass.randompoint()
        vel = S3PointMass.randomvel()
        points.append(S3PointMass(pos, uniform(0.1, 2.0), vel, 0.0))




generate = wx.Button(p, label='Generate!', pos=(1.2*L, 0.8*L))
generate.Bind(wx.EVT_BUTTON, regenerate)



seed()
for i in range(N):
    pos = S3PointMass.randompoint()
    vel = S3PointMass.randomvel()
    points.append(S3PointMass(pos, uniform(0.1, 2.0), vel, 0.0))

GlobalSphere = sphere(pos=(0,0,0), radius=1, color=color.blue, opacity=0.5)

while True:
    dt = 1.0/60.0
    rate(1.0/dt)
    for point in points:
        point.gravitate(points, dt)
    for point in points:
        point.move(dt)
        


