#!/usr/bin/env python3
from __future__ import division
from pylab import *
import os
import h5py
import sys
import datetime, matplotlib
from glob import glob
from write_ndacc_HDF import *
import pytz
from dateutil.relativedelta import relativedelta


def DatetimeToDatenum(dt):
#   ordi = dt.toordinal()
    mdn = dt + datetime.timedelta(days = 366)
    frac = (dt - datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac



def DatenumToDatetime(nb):
    return datetime.datetime(1,1,1) + datetime.timedelta(nb -1 - 366, 0,0,0,0,0,0)



def DatetimeToMJD2K(dt):
#   DAY = datetime.timedelta(1)
    JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
    #J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates
#   day_of_year = (dt - relativedelta(years=dt.year)) + datetime.timedelta(1)  # Jan the 1st is day 1
#   day_of_year = dt.timetuple().tm_yday
    #julian_day2k = (dt - JULIAN_EPOCH + J2000_JD).total_seconds() / 86400 -2400000.5
    julian_day2k = (dt - JULIAN_EPOCH).total_seconds() / 86400 #-2400000.5
    return julian_day2k



def MJD2KToDateTime(mjd):
    JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
    J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates

    dt = JULIAN_EPOCH + datetime.timedelta(mjd)
    return dt


def TimeUnixToMJD2K(li):
    li = li.reshape(li.shape[0])
    print(li,type(li))
    print(li.shape)
    
    out = []
    for i in li:
        out.append(DatetimeToMJD2K(datetime.datetime.fromtimestamp(i * 1e-3, tz=pytz.utc)))

    #   dt =  MJD2KToDateTime(out[-1])
    #   print "We take the timestamp date:", i, out[-1], MJD2KToDateTime(out[-1]), dt.strftime("%Y%m%dt%H%M%Sz")
    return array(out)
    #tmp1 = list(map(lambda x: x*1e-3, li))
    #tmp2 = list(map(datetime.datetime.fromtimestamp, tmp1))
    #tmp2 = list(map(datetime.datetime.fromtimestamp, tmp1))




def LoadH5File(hfil):

    mandatory_variables = {}
    optional_variables = {}
    da =h5py.File(hfil)
    for (a,b) in  da["INSTRUMENT_ATTRIBUTES"].attrs.items():
        print("We are starting with the next attribute",a, b)
        try:
            b = b.decode("ASCII")
        except:
            try:
               # b = b.decode("UTF-8")
                b = b[0].decode("ASCII")
            except:
                pass
        print(a)
        if "Location_Latitude" in a:
            #tmp = b.strip().split(" ")
            try:
                tmp = b.strip().split(" ")
            except:
                tmp = b[0].strip().split(" ")
            print(tmp)
            val = float(tmp[0])
            if "S" == tmp[1][0]:
                val = -val
            print("Latitude", val)
            mandatory_variables["lat"] = array([val])
        if "Location_Longitude" in a:
            #tmp = b.strip().split(" ")
            try:
                tmp = b.strip().split(" ")
            except:
                tmp = b[0].strip().split(" ")
            val = float(tmp[0])
            if "W" == tmp[1][0]:
                val = -val
            print("Longitude", val)
            mandatory_variables["lon"] = array([val])
        if "Location_ASL" in a:
            try:
                tmp = b.strip().split(" ")
            except:
                tmp = b[0].strip().split(" ")
            val = float(tmp[0])
            print("Altitude", val)
            mandatory_variables["elev"] = array([val])
    print("DATETIME/datetime")
    print("INTEGRATION.TIME/integhrs")
    print(da["DATA"]["TIME_MID_UT"][()][0], da["DATA"]["TIME_MID_UT_UNIX"][()][0] ,da["DATA"]["TIME_MID_UT_DATEVEC"][()][0])
    print(da["DATA"]["TIME_MID_UT"][()][0] - 201707200000)
    mandatory_variables["datetime"] = TimeUnixToMJD2K(da["DATA"]["TIME_MID_UT_UNIX"][()])
    mandatory_variables["datetimestart"] = TimeUnixToMJD2K(da["DATA"]["TIME_START_UT_UNIX"][()])
    mandatory_variables["datetimestop"] = TimeUnixToMJD2K(da["DATA"]["TIME_STOP_UT_UNIX"][()])
    dttt = (mandatory_variables["datetimestop"] -  mandatory_variables["datetimestart"]) * 24.

    mandatory_variables["integhrs"] = array(dttt)
    mandatory_variables["z"] = array(da["DATA"]["ALT"][()]).flatten()
    mandatory_variables["o3nd"] = transpose(da["DATA"]["O3ND"][()])
    mandatory_variables["uo3nd"] = transpose(da["DATA"]["O3NDUncert"][()])
    mandatory_variables["uo3ndrand"] = transpose(da["DATA"]["Precision"][()]) * transpose(da["DATA"]["O3ND"][()]) / 100.
    mandatory_variables["uo3ndsyst"] = sqrt( mandatory_variables["uo3nd"] ** 2 -  mandatory_variables["uo3ndrand"] ** 2)
    mandatory_variables["dz"] = transpose(da["DATA"]["O3NDResol"][()])
    mandatory_variables["o3mr"] = transpose(da["DATA"]["O3MR"][()]) * 1e-3
    mandatory_variables["uo3mr"] = transpose(da["DATA"]["O3MRUncert"][()]) * 1e-3
    mandatory_variables["uo3mrrand"] = transpose(da["DATA"]["Precision"][()]) * transpose(da["DATA"]["O3MR"][()]) / 100. * 1e-3
    mandatory_variables["uo3mrsyst"] = sqrt( mandatory_variables["uo3mr"] ** 2 -  mandatory_variables["uo3mrrand"] ** 2)
    mandatory_variables["xp"] =da["DATA"]["Press"][()].flatten()
    mandatory_variables["xt"] =da["DATA"]["Temp"][()].flatten()

    
    # because this file format is made by maniacs with 0 brain cells
#    src = "GEOS-5"
#    import numpy as np
#a    source = np.zeros((len(src)), dtype="S1")
#    for k in range(len(src)):
#        source[k] = src[k]
    source = "ECC Sonde"
    #mandatory_variables["xpsce"] = source #[[source]*len(mandatory_variables["datetime"])]*len(mandatory_variables["z"]) # array([b"GEOS-5"]*len(mandatory_variables["z"]))  # HERE ADD THE ACTUAL SOURCE
    #mandatory_variables["xtsce"] = source  # [[source]*len(mandatory_variables["datetime"])]*len(mandatory_variables["z"]) # array([b"GEOS-5"]*len(mandatory_variables["z"]))  # HERE ADD THE ACTUAL SOURCE
    mandatory_variables["xpsce"] = [source]*len(mandatory_variables["z"]) # array([b"GEOS-5"]*len(mandatory_variables["z"]))  # HERE ADD THE ACTUAL SOURCE
    mandatory_variables["xtsce"] =  [source]*len(mandatory_variables["z"]) # array([b"GEOS-5"]*len(mandatory_variables["z"]))  # HERE ADD THE ACTUAL SOURCE
    return mandatory_variables, optional_variables





def WriteNDAACFile(hfi, meta):
    mand, opt = LoadH5File(hfi)
    write_NDACC_HDF4_O3(meta, mand, opt)




if __name__ == "__main__":
    print( "Usage: ./commonH5tondaac.py H5file.h5 metafile")
    if True:

#   try:
        h5fi = sys.argv[1]
        meta = sys.argv[2]
        print(h5fi, meta)
        WriteNDAACFile(h5fi, meta)
#   except:
#       print "Error"






