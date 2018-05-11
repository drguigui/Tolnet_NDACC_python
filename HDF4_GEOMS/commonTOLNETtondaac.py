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
	#tmp1 = list(map(lambda x: x*1e-3, li))
	#tmp2 = list(map(datetime.datetime.fromtimestamp, tmp1))
	#tmp2 = list(map(datetime.datetime.fromtimestamp, tmp1))




def LoadH5File(hfil):

	mandatory_variables = {}
	optional_variables = {}

	da =h5py.File(hfil)

	
	for (a,b) in  da["INSTRUMENT_ATTRIBUTES"].attrs.iteritems():
		if "Location_Latitude" in a:

			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			if "S" == tmp[1][0]:
				val = -val
			print "Latitude", val
			mandatory_variables["lat"] = array([val])
		if "Location_Longitude" in a:
			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			if "W" == tmp[1][0]:
				val = -val
			print "Longitude", val
			mandatory_variables["lon"] = array([val])
		if "Location_ASL" in a:
			tmp = b[0].strip().split(" ")
			val = float(tmp[0])
			print "Altitude", val
			mandatory_variables["elev"] = array([val])
		#if "Start Time Data Acquisition" in a:
	#		mandatory_variables["datetimestart"] = b
	#	if "End Time Data Acquisition" in a:
	#		mandatory_variables["datetimestop"] = b
	#	i += 1
	#import sys
	#sys.exit()

	print "DATETIME/datetime"
	print "INTEGRATION.TIME/integhrs"
	print da["DATA"]["TIME_MID_UT"].value[0], da["DATA"]["TIME_MID_UT_UNIX"].value[0] ,da["DATA"]["TIME_MID_UT_DATEVEC"].value[0]
	print da["DATA"]["TIME_MID_UT"].value[0] - 201707200000
	mandatory_variables["datetime"] = TimeUnixToMJD2K(da["DATA"]["TIME_MID_UT_UNIX"].value)
	mandatory_variables["datetimestart"] = TimeUnixToMJD2K(da["DATA"]["TIME_START_UT_UNIX"].value)
	mandatory_variables["datetimestop"] = TimeUnixToMJD2K(da["DATA"]["TIME_STOP_UT_UNIX"].value) 


	integ = []
	for i in range(len(mandatory_variables["datetimestart"] )):
		dttt = (mandatory_variables["datetimestop"] -  mandatory_variables["datetimestart"]) * 24.
		integ.append(dttt)



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
	mandatory_variables["xpsce"] = "GEOS-5"
	mandatory_variables["xtsce"] =  "GEOS-5"
	return mandatory_variables, optional_variables

def ReadTolnetFile(fi):
	""" This function reads a TOLnet file that has the same number of altitudes/comment lines for each profiles. Otherwise, it is not possible to create a unique GEOMS curtain plot."""
	f = open(fi,'r')
	da = []
	mandatory_variables = {}
	optional_variables = {}
	for l in f:
	    da.append(l.strip())
	ngh = int(da[0].split(" ")[0])
#	print ngh
	nprof = int(da[2].split(" ")[0])
#	print nprof
	ncol = int(da[3].split(" ")[0])
#	print "ncol", ncol
	ngc = int(da[2 + ngh - 1].split(" ")[0])
#	print "NGC", ngc
	nph = int(da[4 + ngc + ngh - 1].split(" ")[0])
#	print "nph", nph
	nalt = int(da[5 + ngc + ngh - 1].split(" ")[0])
#	print "Nalt", nalt
#	print "date", da[11 - 1 + ngh + ngc]


	locationline =  da[2 + ngh - 1 + 4].strip().split(" ")

	longitude = float(locationline[0].strip(","))
	latitude = float(locationline[1].strip(","))
	altASL = float(locationline[2].strip())

	# 1) Take latitude, longitude, and altitude
	mandatory_variables["lat"] = array([latitude])
	mandatory_variables["lon"] = array([longitude])
	mandatory_variables["elev"] = array([altASL])
	
	altitudes = []
	starttime = []
	endtime = []
	middletime = []
	integrhours = []
	O3s = []   # O3 VMR in PPBV
	O3uncs = [] # O3 VMR uncertainty
	O3nds = []
	O3nduncerts = []
	Resolutions = []
	Precisions = []
	Pressures = []
	Temperatures = []



	# 2) Datetime for all profiles and integration


	for nb in range((nprof)):
    		bprof = 2 + ngh + ngc + (nph + nalt + 2) * nb
    		firstl = bprof + 2 + nph 
    		lastl = bprof + 1 + nph  + nalt 
    		date = bprof + 6
      		tstartdate = DatetimeToMJD2K(datetime.datetime.strptime(da[date][:40].strip(), "%Y-%m-%d, %H:%M:%S").replace(tzinfo=pytz.utc))
      		tenddate = DatetimeToMJD2K(datetime.datetime.strptime(da[date + 1][:40].strip(), "%Y-%m-%d, %H:%M:%S").replace(tzinfo=pytz.utc))
      		tmiddate = DatetimeToMJD2K(datetime.datetime.strptime(da[date + 2][:40].strip(), "%Y-%m-%d, %H:%M:%S").replace(tzinfo=pytz.utc))
		integr = (tenddate - tstartdate) * 24.
		starttime.append(tstartdate)
		endtime.append(tenddate)
		middletime.append(tmiddate)
		integrhours.append(integr)
		altitude = []
		O3 = []   # O3 VMR in PPBV
		O3unc = [] # O3 VMR uncertainty
		O3nd = []
		O3nduncert = []
		Resolution = []
		Precision = []
		Pressure = []
		Temperature = []

		for i in range(firstl, lastl):
			vals = da[i].split(",")
			altitude.append(float(vals[0]))
			O3nd.append(float(vals[1]))
			O3nduncert.append(float(vals[2]))
			Resolution.append(float(vals[3]))
			Precision.append(float(vals[4]))
			O3.append(float(vals[6]))
			O3unc.append(float(vals[7]))
			Pressure.append(float(vals[8]))
			Temperature.append(float(vals[10]))

		O3s.append(O3)
		O3uncs.append(O3unc)
		O3nds.append(O3nd)
		O3nduncerts.append(O3nduncert)
		Resolutions.append(Resolution)
		Precisions.append(Precision)
		Pressures.append(Pressure)
		Temperatures.append(Temperature)




		altitudes.append(altitude)
	#	return array(altitude), array(O3), array(O3unc)

	integrhours = array(integrhours)
	O3s = array(O3s)
	O3uncs = array(O3uncs)
	O3nds = array(O3nds)
	O3nduncerts = array(O3nduncerts)
	Resolutions =array(Resolutions)
	Precisions = array(Precisions)
	Pressures = array(Pressures)
	Temperatures = array(Temperatures)
	mandatory_variables["datetime"] = array(middletime)
	mandatory_variables["datetimestart"] = array(starttime)
	mandatory_variables["datetimestop"] = array(endtime)
	mandatory_variables["integhrs"] = integrhours
	# 3) altitude grid
	#print altitude
	mandatory_variables["z"] = altitudes[0]
	# 4) o3nd, uo3nd....
	mandatory_variables["o3nd"] = (O3nds)
	mandatory_variables["uo3nd"] = (O3nduncerts)
	mandatory_variables["uo3ndrand"] = Precisions * O3nds / 100. 
	#transpose(da["DATA"]["Precision"].value) * transpose(da["DATA"]["O3ND"].value) / 100.
	mandatory_variables["uo3ndsyst"] = sqrt( mandatory_variables["uo3nd"] ** 2 -  mandatory_variables["uo3ndrand"] ** 2)
	mandatory_variables["dz"] = Resolutions
	mandatory_variables["o3mr"] = O3s * 1e-3
	mandatory_variables["uo3mr"] = O3uncs * 1e-3
	mandatory_variables["uo3mrrand"] = Precisions * O3s / 100. * 1e-3
	mandatory_variables["uo3mrsyst"] = sqrt( mandatory_variables["uo3mr"] ** 2 -  mandatory_variables["uo3mrrand"] ** 2)
	mandatory_variables["xp"] = Pressures[0]
	mandatory_variables["xt"] = Temperatures[0]


	# 5) dz
	# 6) Mixing ratio
	# 7) Pressure
	# 8) temperature




	mandatory_variables["xpsce"] = "GEOS-5"  # HERE ADD THE ACTUAL SOURCE
	mandatory_variables["xtsce"] =  "GEOS-5"  # HERE ADD THE ACTUAL SOURCE
	print mandatory_variables
	return mandatory_variables, optional_variables




def WriteNDAACFile(hfi, meta):
	mand, opt = LoadH5File(hfi)
	write_NDACC_HDF4_O3(meta, mand, opt)

def WriteTolnetNDAACFile(tofi, meta):
	mand, opt = ReadTolnetFile(tofi)
	write_NDACC_HDF4_O3(meta, mand, opt)



if __name__ == "__main__":
	print "Usage: ./commonTOLNETtondaac.py TOLNETfile.txt metafile"
	if True:

#	try:
		tofi = sys.argv[1]
		meta = sys.argv[2]
		print tofi, meta
		#WriteNDAACFile(h5fi, meta)
		WriteTolnetNDAACFile(tofi, meta)
#	except:
#		print "Error"






