#!/usr/bin/env python3
#from __future__ import division
from pylab import *
import h5py
import sys
import pytz
print("For quick view only")
print("Disclaimer: known bugs on maximum value of 100ppbv and if data are cut in time")
print("Contact Guillaume.P.Gronoff@nasa.gov for better plots")


da =h5py.File(sys.argv[1])

def DatetimeToDatenum(dt):
#   ordi = dt.toordinal()
    mdn = dt + datetime.timedelta(days = 366)
    frac = (dt - datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac



def DatenumToDatetime(nb):
    return datetime.datetime(1,1,1) + datetime.timedelta(nb -1 - 366, 0,0,0,0,0,0)


print(da.items())

i = 0
for (a,b) in  da["INSTRUMENT_ATTRIBUTES"].attrs.items(): #.iteritems():
    print( i, a, b)
    i += 1
print( da["DATA"].items())

data = da["DATA"]["O3MR"][()]
data[data>100] = 100
data[data<0] = 0
data = data[::-1]
print(da["DATA"]["TIME_START_UT_UNIX"][()][0])
time1 = datetime.datetime.fromtimestamp(da["DATA"]["TIME_START_UT_UNIX"][()][0][0] / 1e3, pytz.utc)
time2 = datetime.datetime.fromtimestamp(da["DATA"]["TIME_STOP_UT_UNIX"][()][-1][0] / 1e3, pytz.utc)
xlims = [time1, time2]
x_lims = matplotlib.dates.date2num(xlims)


fig=figure()
ax = subplot(221)


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
alt=  da["DATA"]["ALT"][()]
y_lims = [alt[0][0], alt[-1][0]]
print(x_lims, y_lims)
ax.imshow((abs(data)), interpolation='nearest', aspect='auto', extent= [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]], cmap = ncmap, norm=nnorm)

ax.xaxis_date()
fig.autofmt_xdate()

colorbar(ax.get_children()[-2], ax=ax)




#cb2 = mpl.colorbar.ColorbarBase(ax =ax, cmap=cmap, norm=norm, boundaries=[0] + bounds + [600], extend='both', ticks=bounds, spacing='proportional', orientation='vertical')
#cb2.set_label("O3 [PPBV]") #( Discrete intervals, some other units')
#colorbar(ax.get_children()[2], ax=ax, cmap=cmap, norm=norm)
#colorbar(ax.get_children()[-2], ticks= linspace(0,100,11)) #,boundaries=[0] + bounds + [200])






try:
    ax = subplot(222)
    title("Aerosols Beta")
    ylabel("Altitude [m]")
    xlabel("Time [UTC]")
    data = da["DATA"]["AerosolsBeta"][()]
    data[data>1] = 1
    data[data<0] = 0
    data = data[::-1]
#x_lims = list(map(datetime.datetime.fromtimestamp, [da["DATA"]["TIME_START_UT_UNIX"].value[0] / 1e3,  da["DATA"]["TIME_STOP_UT_UNIX"].value[-1] /1e3]))
#x_lims = matplotlib.dates.date2num(x_lims)

    alt=  da["DATA"]["ALT"][()]
    y_lims = [alt[0][0], alt[-1][0]]

    ax.imshow((abs(data)), interpolation='nearest', aspect='auto', extent= [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
    ax.xaxis_date()
    fig.autofmt_xdate()
#try:
    colorbar(ax.get_children()[-2], ax=ax)
#

#try:
except:
    print("No aerosols")
#colorbar(ax.get_children()[2], ax=ax)



ax = subplot(224)
title("Vertical Resolution [m]")
ylabel("Altitude [m]")
xlabel("Time [UTC]")
data = da["DATA"]["O3NDResol"][()]
#data[data>500] = 500
data[data<0] = 0
data = data[::-1]
#x_lims = list(map(datetime.datetime.fromtimestamp, [da["DATA"]["TIME_START_UT_UNIX"].value[0] / 1e3,  da["DATA"]["TIME_STOP_UT_UNIX"].value[-1] /1e3]))
#x_lims = matplotlib.dates.date2num(x_lims)

alt=  da["DATA"]["ALT"][()]
y_lims = [alt[0][0], alt[-1][0]]

ax.imshow((abs(data)), interpolation='nearest', aspect='auto', extent= [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
ax.xaxis_date()
fig.autofmt_xdate()
#try:
colorbar(ax.get_children()[-2], ax=ax)
#except:
#colorbar(ax.get_children()[2], ax=ax)

ax = subplot(223)

title("O3 Mixing Ratio Uncertainty [PPBV]")
ylabel("Altitude [m]")
xlabel("Time [UTC]")
data = da["DATA"]["O3MRUncert"][()]
data[data>100] = 100
data[data<0] = 0
data = data[::-1]
alt=  da["DATA"]["ALT"][()]
y_lims = [alt[0][0], alt[-1][0]]
ax.imshow((abs(data)), interpolation='nearest', aspect='auto', extent= [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
ax.xaxis_date()
fig.autofmt_xdate()
#atry:
#   colorbar(ax.get_children()[-2], ax=ax)
#except:
colorbar(ax.get_children()[-2], ax=ax)
show()
