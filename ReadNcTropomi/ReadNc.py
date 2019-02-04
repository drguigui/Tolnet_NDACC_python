#!/usr/bin/env python
from __future__ import division
from pylab import *
import h5py
# Needs to be python3 in the near future
import sys

da = h5py.File(sys.argv[1])

print da.items()
for item in da['PRODUCT'].items():
	print item
print da['PRODUCT']['latitude'].value
print da['PRODUCT']['longitude'].value
op = da['PRODUCT']['ozone_profile'].value
#print op
#print op.shape
#aozone_profile_subcolumns
#ozone_profile_precision
#print da['METADATA'].items()
for item in da['PRODUCT']['SUPPORT_DATA'].items():
	print item
print da['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']

for item in da['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].items():
	print item
print da['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']["altitude"].value
