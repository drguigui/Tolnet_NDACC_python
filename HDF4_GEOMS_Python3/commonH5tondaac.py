#!/usr/bin/env python
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
#	ordi = dt.toordinal()
	mdn = dt + datetime.timedelta(days = 366)
	frac = (dt - datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
	return mdn.toordinal() + frac



def DatenumToDatetime(nb):
	return datetime.datetime(1,1,1) + datetime.timedelta(nb -1 - 366, 0,0,0,0,0,0)



def DatetimeToMJD2K(dt):
#	DAY = datetime.timedelta(1)
	JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
	#J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates
#	day_of_year = (dt - relativedelta(years=dt.year)) + datetime.timedelta(1)  # Jan the 1st is day 1
#	day_of_year = dt.timetuple().tm_yday
	#julian_day2k = (dt - JULIAN_EPOCH + J2000_JD).total_seconds() / 86400 -2400000.5
	julian_day2k = (dt - JULIAN_EPOCH).total_seconds() / 86400 #-2400000.5
	return julian_day2k



def MJD2KToDateTime(mjd):
	JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
	J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates

	dt = JULIAN_EPOCH + datetime.timedelta(mjd)
	return dt


def TimeUnixToMJD2K(li):
	out = []
	for i in li:
		out.append(DatetimeToMJD2K(datetime.datetime.fromtimestamp(i * 1e-3, tz=pytz.utc)))

	#	dt =  MJD2KToDateTime(out[-1])
	#	print "We take the timestamp date:", i, out[-1], MJD2KToDateTime(out[-1]), dt.strftime("%Y%m%dt%H%M%Sz")
	return array(out)


def LoadH5File(hfil):

	mandatory_variables = {}
	optional_variables = {}
	da =h5py.File(hfil)
	for (a,b) in  da["INSTRUMENT_ATTRIBUTES"].attrs.items():
		b[0] = b[0].decode("ASCII")
		if "Location_Latitude" in a:

			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			if "S" == tmp[1][0]:
				val = -val
			print("Latitude", val)
			mandatory_variables["lat"] = array([val])
		if "Location_Longitude" in a:
			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			if "W" == tmp[1][0]:
				val = -val
			print("Longitude", val)
			mandatory_variables["lon"] = array([val])
		if "Location_ASL" in a:
			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			print("Altitude", val)
			mandatory_variables["elev"] = array([val])
	print("DATETIME/datetime")
	print("INTEGRATION.TIME/integhrs")
	print(da["DATA"]["TIME_MID_UT"].value[0], da["DATA"]["TIME_MID_UT_UNIX"].value[0] ,da["DATA"]["TIME_MID_UT_DATEVEC"].value[0])
	print(da["DATA"]["TIME_MID_UT"].value[0] - 201707200000)
	mandatory_variables["datetime"] = TimeUnixToMJD2K(da["DATA"]["TIME_MID_UT_UNIX"].value)
	mandatory_variables["datetimestart"] = TimeUnixToMJD2K(da["DATA"]["TIME_START_UT_UNIX"].value)
	mandatory_variables["datetimestop"] = TimeUnixToMJD2K(da["DATA"]["TIME_STOP_UT_UNIX"].value) 
	dttt = (mandatory_variables["datetimestop"] -  mandatory_variables["datetimestart"]) * 24.

	mandatory_variables["integhrs"] = array(dttt)
	mandatory_variables["z"] = da["DATA"]["ALT"].value
	mandatory_variables["o3nd"] = transpose(da["DATA"]["O3ND"].value)
	mandatory_variables["uo3nd"] = transpose(da["DATA"]["O3NDUncert"].value)
	mandatory_variables["uo3ndrand"] = transpose(da["DATA"]["Precision"].value) * transpose(da["DATA"]["O3ND"].value) / 100.
	mandatory_variables["uo3ndsyst"] = sqrt( mandatory_variables["uo3nd"] ** 2 -  mandatory_variables["uo3ndrand"] ** 2)
	mandatory_variables["dz"] = transpose(da["DATA"]["O3NDResol"].value)
	mandatory_variables["o3mr"] = transpose(da["DATA"]["O3MR"].value) * 1e-3
	mandatory_variables["uo3mr"] = transpose(da["DATA"]["O3MRUncert"].value) * 1e-3
	mandatory_variables["uo3mrrand"] = transpose(da["DATA"]["Precision"].value) * transpose(da["DATA"]["O3MR"].value) / 100. * 1e-3
	mandatory_variables["uo3mrsyst"] = sqrt( mandatory_variables["uo3mr"] ** 2 -  mandatory_variables["uo3mrrand"] ** 2)
	mandatory_variables["xp"] =da["DATA"]["Press"].value
	mandatory_variables["xt"] =da["DATA"]["Temp"].value
	mandatory_variables["xpsce"] = "GEOS-5"  # HERE ADD THE ACTUAL SOURCE
	mandatory_variables["xtsce"] =  "GEOS-5"  # HERE ADD THE ACTUAL SOURCE
	return mandatory_variables, optional_variables





def WriteNDAACFile(hfi, meta):
	mand, opt = LoadH5File(hfi)
	write_NDACC_HDF4_O3(meta, mand, opt)




if __name__ == "__main__":
	print( "Usage: ./commonH5tondaac.py H5file.h5 metafile")
	if True:

#	try:
		h5fi = sys.argv[1]
		meta = sys.argv[2]
		print(h5fi, meta)
		WriteNDAACFile(h5fi, meta)
#	except:
#		print "Error"






