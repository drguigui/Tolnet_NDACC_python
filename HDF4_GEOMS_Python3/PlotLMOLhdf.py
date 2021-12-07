#!/usr/bin/env python
from __future__ import division
from pylab import *
#import h5py
import pyhdf
from pyhdf import SD
from pyhdf.SD import *
#from pyhdf import SDC
import sys
import pytz
import matplotlib.dates as mdates
print("For quick view only")
print("Disclaimer: known bugs on maximum value of 100ppbv and if data are cut in time")
print("Contact Guillaume.P.Gronoff@nasa.gov for better plots")



file_name = sys.argv[1]
da = SD(file_name, SDC.READ)

print(da.info())

datasets_dic = da.datasets()

for idx,sds in enumerate(datasets_dic.keys()):
	print(idx,sds)

def MJD2KToDateTimeString(mjd):
	JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
	J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates

	dt = JULIAN_EPOCH + datetime.timedelta(mjd)
	print("DT=", dt)
	return dt


def DatetimeToDatenum(dt):
#	ordi = dt.toordinal()
	mdn = dt + datetime.timedelta(days = 366)
	frac = (dt - datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
	return mdn.toordinal() + frac



def DatenumToDatetime(nb):
	return datetime.datetime(1,1,1) + datetime.timedelta(nb -1 - 366, 0,0,0,0,0,0)




data2D = da.select("O3.MIXING.RATIO.VOLUME_DERIVED")
data = data2D[:,:] *1e3
data = data.transpose()
data[data>100] = 100
data[data<0] = 0
data = data[::-1]
time1 = datetime.datetime.fromtimestamp(da.select("DATETIME.START")[0] / 1e3, pytz.utc)

time1 = MJD2KToDateTimeString(da.select("DATETIME.START")[0])
dastop = da.select("DATETIME.STOP")
dastop = dastop[:]
time2 = MJD2KToDateTimeString(dastop[-1])
xlims = [time1, time2]
x_lims = matplotlib.dates.date2num(xlims)

fig=figure()
ax = subplot(111)


ncolors = [array([255,  140,  255]) / 255.,  array([221,  111,  242]) / 255., array([187,  82,  229]) / 255.,  array([153,  53,  216]) / 255.,  array([119,  24,  203]) / 255.,  array([0,  0,  187]) / 255.,  array([0,  44,  204]) / 255.,  array([0,  88,  221]) / 255.,  array([0,  132,  238]) / 255.,  array([0,  175,  255]) / 255.,  array([0,  235,  255]) / 255.,  array([39,  255,  215]) / 255.,  array([99,  255,  155]) / 255.,  array([163,  255,  91]) / 255.,  array([211,  255,  43]) / 255.,  array([255,  255,  0]) / 255.,  array([255,  207,  0]) / 255.,  array([255,  159,  0]) / 255.,  array([255,  111,  0]) / 255.,  array([255,  63,  0]) / 255.,  array([255,  0,  0]) / 255.,  array([216,  0,  15]) / 255.,  array([178,  0,  31]) / 255.,  array([140,  0,  47]) / 255.,  array([102,  0,  63]) / 255.,  array([52,  52,  52]) / 255.,  array([96,  96,  96]) / 255.,  array([140,  140,  140]) / 255.,  array([184,  184,  184]) / 255.,  array([228,  228,  228]) / 255., [1.,1.,1.] ]

import matplotlib as mpl
import matplotlib.pyplot as plt
ncmap = mpl.colors.ListedColormap(ncolors)
ncmap.set_under([1,1,1])
ncmap.set_over([0,0,0])
bounds = [0.001, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 125, 150, 200, 300,600]
#print len(bounds), len(ncolors)
nnorm = mpl.colors.BoundaryNorm(bounds, ncmap.N)


title("O3 Mixing Ratio [PPBV]")
ylabel("Altitude [m]")
xlabel("Time [UTC]")
alt=  da.select("ALTITUDE") 
#alt = alt[:] 
#alt = alt[::-1]
y_lims = [alt[0][0], alt[-1][0]]
#print(type(data))
#print(data)

print(x_lims, y_lims)
ax.imshow((abs(data)), interpolation='nearest', aspect='auto', extent= [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]], cmap = ncmap, norm=nnorm)
#ax.imshow((abs(data)), interpolation='nearest', aspect='auto', cmap = ncmap, norm=nnorm)

ax.xaxis_date()

ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H"))
ax.xaxis.set_minor_formatter(mdates.DateFormatter("%H"))

fig.autofmt_xdate()

colorbar(ax.get_children()[-2], ax=ax)
print(data)

show()

sys.exit()



