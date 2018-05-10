#!/usr/bin/env python
from __future__ import division
from pylab import *
import random


def conv(irconvol,vCoef):
    nc = (len(vCoef) - 1) // 2
    ir = zeros(len(irconvol))
    for j in range(nc, len(irconvol) - nc ):
       ir[j] = sum(array(irconvol[j-nc :j+nc + 1]) * array(vCoef)) # We have to do the convolution like this
    return ir

def GenFact(a, b):
    # integers
    gf = 1
    for i in range(a-b+1, a + 1):
        gf *= i
    return gf

def GramPoly(i, m , k, s):
    if k >  0:
        return (4 * k - 2) / (k * (2 * m - k + 1)) * (i * GramPoly(i, m, k - 1, s) + s * GramPoly(i, m, k - 1, s - 1)) -((k - 1) * (2 * m + k)) / (k * (2 * m - k + 1)) * GramPoly(i, m, k - 2, s)
        print "probleme",k
        return 0
    if k == 0 and s == 0:
            return 1
    return 0

def FilterWeight(t, m, n, s):
    """ filter weight for the t Least Square (t = 0: Savitzky Golay)
    s'th derivative
    m points: 2 * m + 1
    n: order"""

    halfweight = []
    for i in range(-m, m + 1):
        um = 0
        for k in range(0, n + 1):
            um += (2 * k + 1) * (GenFact(2 * m, k) / GenFact(2 * m + k + 1, k + 1)) * GramPoly(i, m, k, 0) * GramPoly(t, m, k, s)
        halfweight.append(um)
    return halfweight


if '__main__' == __name__:

#    print "Filters"
#    print array(FilterWeight(-2, 2, 2, 0)) * 35
#    print array(FilterWeight(-1, 2, 2, 0)) * 35
#    print array(FilterWeight(0, 2, 2, 0)) * 35
#    print array(FilterWeight(1, 2, 2, 0)) * 35
#    print array(FilterWeight(2, 2, 2, 0)) * 35
#    print array(FilterWeight(0, 3, 2, 0)) * 21
#    print array(FilterWeight(0, 2, 2, 1)) * 10
#    print array(FilterWeight(0, 3, 2, 1)) * 28
#    print array(FilterWeight(0, 3, 3, 1)) * 252


#    print "blablablalb"
#    print array(FilterWeight(0, 2, 2, 0)) *35
#    print array(FilterWeight(0, 3, 2, 0)) *21
#    print array(FilterWeight(0, 4, 2, 0)) *231
#    print array(FilterWeight(0, 5, 2, 0)) * 429
#    print array(FilterWeight(0, 3, 4, 0)) * 231
#    print array(FilterWeight(0, 9, 2, 1)) * 570
#    print array(FilterWeight(0, 10, 3, 1)) * 3634092


    x =linspace(-10000,10000, 10000*2)
    
    y = (x / 100) ** 2
    y2 = y.copy()
    for i in range(len(y2)):
        y2[i] = random.gauss(y2[i], sqrt(y2[i]) / 10)
    plot(x, y2)
    #plot(x, y)
# we do a smoothinFilterWeight(0, 2, 2, 0)g over 21 points
    weight = FilterWeight(0, 10, 2, 0)

    y3 = conv(y2, weight)
    plot(x, y3)


# We do a first order derivative over 5 points
    weight = FilterWeight(0, 40, 4, 1)
    print weight
    y4 = conv(y2, weight)
    plot(x, y4 * 10000/2)


# We do a second order derivative over 5 points
    weight = FilterWeight(0, 40, 4, 2)
    print weight
    y4 = conv(y2, weight)
    plot(x, y4 * 10000/2)

    show()

    clf()
# We do a second order derivative over 5 points
    weight = FilterWeight(0, 2, 4, 2)
    print weight
    y3 = conv(y, weight)
    plot(x, y3)


    show()



