#!/usr/bin/env python
from __future__ import division
from pylab import *

from filtering import FilterWeight, conv
from ResolutionDF import *
from ResolutionIR import *



def CheckData(weight):
    impulse = zeros((2400))
    rez = NDAAC_ResolIr(100, weight, impulse)
    print "dzIR", rez[1]

    impulse = ones((2400))

    rez = NDAAC_ResolDf(100, weight, impulse)
    print "dzDF", rez[1]



print "We First check the degree 0"

print "3 points"
weight = FilterWeight(0,1,0,0)
print weight
CheckData(weight)

print "5 points"
weight = FilterWeight(0,2,0,0)
print weight
CheckData(weight)

print "9 points"
weight = FilterWeight(0,4,0,0)
print weight
CheckData(weight)

print "17 points"
weight = FilterWeight(0,8,0,0)
print weight
CheckData(weight)

print "We  check the degree 0, hanning"

print "3 points"
weight = [0.,1.,0.]
print weight
CheckData(weight)

print "5 points"

weight = [0.00000000, 0.25000000, 0.50000000, 0.25000000, 0.00000000]
print weight
CheckData(weight)

print "9 points"

weight =[0.00000000,0.036611652,0.12500000,0.21338835,0.25000000,0.21338835,0.12500000,0.036611652,0.00000000]

print weight
CheckData(weight)

print "17 points"
weight= [0.00000000,0.0047575292,0.018305826,0.038582285,0.062500000,0.086417715,0.10669417,0.12024247,0.12500000,0.12024247,0.10669417,0.086417715,0.062500000,0.038582285,0.018305826,0.0047575292,0.00000000]
print weight
CheckData(weight)




print "We  check the degree 2"

print "3 points"
weight = FilterWeight(0,1,2,0)
print weight
CheckData(weight)

print "5 points"
weight = FilterWeight(0,2,2,0)
print weight
CheckData(weight)

print "9 points"
weight = FilterWeight(0,4,2,0)
print weight
CheckData(weight)

print "17 points"
weight = FilterWeight(0,8,2,0)
print weight
CheckData(weight)




print "We check the degree 2, hanning"

print "3 points"
weight = [0,1,0]
print weight
CheckData(weight)

print "5 points"

weight = [-0.00000000,0.20689655,0.58620690,0.20689655,-0.00000000]
print weight
CheckData(weight)

print "9 points"

weight =[-0.00000000,0.010552849,0.10036839,0.23723940,0.30367873,0.23723940,0.10036839,0.010552849,-0.00000000]
print weight
CheckData(weight)

print "17 points"
weight= [-0.00000000,-0.00082412208,0.0036995271,0.020050227,0.048719477,0.084828255,0.12013351,0.14580285,0.15518056,0.14580285,0.12013351,0.084828255,0.048719477,0.020050227,0.0036995271,-0.00082412208,-0.00000000]
print weight
CheckData(weight)



print "We check the degree 1 deriv 1"

print "3 points"
weight = FilterWeight(0,1,1,1)
print weight
CheckData(weight)

print "5 points"
weight = FilterWeight(0,2,1,1)
print weight
CheckData(weight)

print "9 points"
weight = FilterWeight(0,4,1,1)
print weight
CheckData(weight)

print "17 points"
weight = FilterWeight(0,8,1,1)
print weight
CheckData(weight)


print "We check the degree 3 deriv 1"

#print "3 points"
#weight = FilterWeight(0,1,3,1)
#print weight

print "5 points"
weight = FilterWeight(0,2,3,1)
print weight
CheckData(weight)

print "9 points"
weight = FilterWeight(0,4,3,1)
print weight
CheckData(weight)

print "17 points"
weight = FilterWeight(0,8,3,1)
print weight
CheckData(weight)




