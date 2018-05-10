#!/usr/bin/env python
from __future__ import division
from pylab import *
import sys 
from pyhdf.SD import *

# Example file to read 

try:
	fi = sys.argv[1]
	print "We are reading the file ",fi
except:
	print "Please provide a valid hdf file"

d = SD(fi, SDC.READ)

print d.attributes()
print d.datasets()



for da in d.datasets().keys():
	g = d.select(da)
	print ""
	print ""
	print "DATASET : ",da 
	print "========="
	print g.attributes()
	if da =="DATETIME.START":
		print "datetime"
		print g[:]
	if da =="ALTITUDE":
		print "Altitude"
		print g[:100]
	if da == "O3.MIXING.RATIO.VOLUME_DERIVED":
		print g[0, :100]
		print g[0, -101:-1]


