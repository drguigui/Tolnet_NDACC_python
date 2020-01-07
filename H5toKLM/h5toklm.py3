#!/usr/bin/env python
from __future__ import division
import sys
import pytz
print( "For quick view only")
print( "Disclaimer: known bugs on maximum value of 100ppbv and if data are cut in time")
print( "Contact Guillaume.P.Gronoff@nasa.gov for better plots")
import simplekml
import h5py
from pylab import *
def DatetimeToDatenum(dt):
#	ordi = dt.toordinal()
    mdn = dt + datetime.timedelta(days = 366)
    frac = (dt - datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac



def DatenumToDatetime(nb):
    return datetime.datetime(1,1,1) + datetime.timedelta(nb -1 - 366, 0,0,0,0,0,0)



redarr = [255,221,187,153,119,0,0,0,0,0,0,39,99,163,211,255,255,255,255,255,255,216,178,140,102,52,96,140,184,228]
greenarr =[140,111,82,53,24,0,44,88,132,175,235,255,255,255,255,255,207,159,111,63,0,0,0,0,0,52,96,140,184,228]
bluearr = [255,242,229,216,203,187,204,221,238,255,255,215,155,91,43,0,0,0,0,0,0,15,31,47,63,52,96,140,184,228]
valarr = [4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,125,150,200,300,600]


def ValueToColor(val):

    for i in range(len(valarr)):
        if val < valarr[i]:
            opacity = 200
            va = hex(opacity*16**6 + bluearr[i] * 16 **4 +  greenarr[i] * 16 **2 + redarr[i])
            return va[2:]
    #if(val < 40):
    #    return simplekml.Color.green
    #if(val < 70):
    #    return simplekml.Color.orange
    #return simplekml.Color.red




da =h5py.File(sys.argv[1])
data = da["DATA"]["O3MR"].value # Mixing ratio in PPBV
alt=  da["DATA"]["ALT"].value # Altitude in m

deltaalt = abs(alt[1] - alt[0]) / 2.

try:
    l1 = da["INSTRUMENT_ATTRIBUTES"].attrs["Location_Latitude"][0].decode().strip()
    l2 = da["INSTRUMENT_ATTRIBUTES"].attrs["Location_Longitude"][0].decode().strip()
    #print((l1))
    lat= abs(float(l1[:-1]))
    lon = abs(float(l2[:-1]))
    if "S" in l1:
        lat = -lat
    if "W" in l2:
        lon = -lon
    print( lat, lon, l1, l2)
except:
    l1 = da["INSTRUMENT_ATTRIBUTES"].attrs["Location_Latitude"].decode().strip()
    l2 = da["INSTRUMENT_ATTRIBUTES"].attrs["Location_Longitude"].decode().strip()
    lat= float(l1[:-1])
    lon = float(l2[:-1])
    if "S" in l1:
        lat = -lat
    if "W" in l2:
        lon = -lon
    print( lat, lon, l1, l2)
kml = simplekml.Kml()

#data = data[::-1]
for i in range(len(da["DATA"]["TIME_START_UT_UNIX"])):
    print( i)
    time1 = datetime.datetime.fromtimestamp(da["DATA"]["TIME_START_UT_UNIX"].value[i] / 1e3, pytz.utc).isoformat()
    time2 = datetime.datetime.fromtimestamp(da["DATA"]["TIME_STOP_UT_UNIX"].value[i] / 1e3, pytz.utc).isoformat()

    MR = data[:, i] # Hopefully it is this one


    for a in range(len(alt)):
        if MR[a] > 0 and MR[a] < 9999:
            print(lat,lon, (alt[a] - deltaalt)[0])
            lin = kml.newlinestring(name="LMOL O3 [ppbv]", description="altitude: "+str(alt[a][0])+" m <br/> value: "+ str(int(MR[a])) +" ppbv", coords=[(lon,lat, (alt[a] - deltaalt)[0]), (lon,lat, (alt[a] + deltaalt)[0])])
            lin.linestyle.width= 20
            lin.linestyle.color= ValueToColor(MR[a])
            lin.timespan.begin = time1
            lin.timespan.end = time2
            lin.altitudemode = simplekml.AltitudeMode.relativetoground
            lin.extrude=0



kml.save(sys.argv[1][:-2]+ "kml")

