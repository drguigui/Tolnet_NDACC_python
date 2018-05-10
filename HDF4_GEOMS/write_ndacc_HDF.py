#!/usr/bin/env python
from __future__ import division
from pylab import *

import os
from pyhdf.SD import *
import re
import datetime
import pytz


def MJD2KToDateTimeString(mjd):
	JULIAN_EPOCH = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) # noon (the epoch name is unrelated)
	J2000_JD = datetime.timedelta(2451545) # julian epoch in julian dates

	dt = JULIAN_EPOCH + datetime.timedelta(mjd)
	return dt.strftime("%Y%m%dT%H%M%SZ")


def BoundarySave(da, value):
    maxi = da["VAR_VALID_MAX"]
    mini = da["VAR_VALID_MIN"]
    fill = da["VAR_FILL_VALUE"]


    typ = da["VAR_DATA_TYPE"]
    if "STRING" == typ:
        return value
    
    if "REAL" == typ:
        maxi = float(maxi)
        mini = float(mini)
        fill = float(fill)
   
    if "DOUBLE" == typ:
        maxi = float(maxi)
        mini = float(mini)
        fill = float(fill)

    if "SHORT" == typ:
        maxi = int(maxi)
        mini = int(mini)
        fill = int(fill)

    if "INTEGER" == typ:
        maxi = long(maxi)
        mini = long(mini)
        fill = long(fill)



    try:
        arr = array(value)
	arr[isnan(arr)] = fill
        arr[arr > maxi] = fill
        arr[arr < mini] = fill
        return arr
    except:
        if (arr > maxi) or (arr < mini):
            return fill
        return arr



def write_NDACC_HDF4_O3(metafile, mandatory_variables, optional_variables = {},pathout = "./"):
    

    # We read the meta-data file

    fi = open(metafile,'r')

    li = fi.readline()
    gattname = []
    gattval = []
    attributes = {}
    while not (re.match("^! Variable Attributes", li)):
        li = li.strip()
        l = li.split("=")
        if(2 == len(l)):
            gattname.append(l[0].strip())
            val = l[1].strip()
            if '' == val:
                val = "  " # To prevent any problem with null strings
            gattval.append(val)
            attributes[l[0].strip()] = val
        li = fi.readline()

    varlist = {}
    varname = ""
    while li.strip() != "!END":
        li = li.strip()
	#print "LINE:" , li
	if li[0] == "!":
		li = fi.readline()
		continue
        l = li.split("=")
        if l[0].strip() == "VAR_NAME":
            varname = l[1].strip()
            varlist[varname] = {}
        else:
            val = l[1].strip()
            if '' == val:
                val = "  " # To prevent any problem with null strings
            varlist[varname][l[0].strip()] = val
	li = fi.readline()
    # We create the output file name
    output = pathout 
    fileout = ""
    attributes["DATA_START_DATE"] = MJD2KToDateTimeString(mandatory_variables["datetimestart"][0])
    attributes["DATA_STOP_DATE"] = MJD2KToDateTimeString(mandatory_variables["datetimestop"][-1])
    attributes["FILE_GENERATION_DATE"]= datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")

    tmp = attributes["DATA_DISCIPLINE"].split(";")
    fileout += tmp[-1].lower() +"_"
    fileout += attributes["DATA_SOURCE"].lower() + "_"
    fileout += attributes["DATA_LOCATION"].lower().split(";")[0] + "_"
    fileout += attributes["DATA_START_DATE"].lower() +"_"
    fileout += attributes["DATA_STOP_DATE"].lower() + "_"
    fileout += attributes["DATA_FILE_VERSION"].lower()
    fileout += ".hdf"
    output += fileout
    attributes["FILE_NAME"] = fileout
    print "Output file: ", output
    
    # We create the main output file
    d = SD(output, SDC.WRITE|SDC.CREATE)

#    for varname in varlist.keys():
#       maxi = varlist[varname]["VAR_VALID_MAX"]
#       mini = varlist[varname]["VAR_VALID_MIN"]
#       fill = varlist[varname]["VAR_FILL_VALUE"]
#       typ = varlist[varname]["VAR_DATA_TYPE"]
      # if "STRING" == typ:
      #    return value
    
#       if "REAL" == typ:
#           varlist[varname]["VAR_VALID_MAX"] = float(maxi)
#           varlist[varname]["VAR_VALID_MIN"] = float(mini)
#           varlist[varname]["VAR_FILL_VALUE"] = float(fill)
    
#       if "DOUBLE" == typ:
#           varlist[varname]["VAR_VALID_MAX"] = float(maxi)
#           varlist[varname]["VAR_VALID_MIN"] = float(mini)
#           varlist[varname]["VAR_FILL_VALUE"] = float(fill)

#       if "SHORT" == typ:
#           varlist[varname]["VAR_VALID_MAX"] = int(maxi)
#           varlist[varname]["VAR_VALID_MIN"] = int(mini)
#           varlist[varname]["VAR_FILL_VALUE"] = int(fill)

#       if "INTEGER" == typ:
#           varlist[varname]["VAR_VALID_MAX"] = long(maxi)
#           varlist[varname]["VAR_VALID_MIN"] = long(mini)
#           varlist[varname]["VAR_FILL_VALUE"] = long(fill)






    for varname in varlist.keys():
        value = []
        # We add the mandatory data into the varlist values  


        if "LATITUDE.INSTRUMENT" == varname:
            value =  mandatory_variables["lat"]
        if "LONGITUDE.INSTRUMENT" == varname:
            value =  mandatory_variables["lon"]
        if "ALTITUDE.INSTRUMENT" == varname:
            value =  mandatory_variables["elev"]
        if "DATETIME" == varname:
            value =  mandatory_variables["datetime"]
        if "DATETIME.START" == varname:
            value =  mandatory_variables["datetimestart"]
        if "DATETIME.STOP" == varname:
            value =  mandatory_variables["datetimestop"]
        if "INTEGRATION.TIME" == varname:
            value =  mandatory_variables["integhrs"]
        if "ALTITUDE" == varname:
            value =  mandatory_variables["z"]
        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL" == varname:
            value =  mandatory_variables["o3nd"]
        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD" == varname:
            value =  mandatory_variables["uo3nd"]
        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD" == varname:
            value =  mandatory_variables["uo3ndrand"]
        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD" == varname:
            value =  mandatory_variables["uo3ndsyst"]
	#    varname = "O3NDADUSS"
        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM" == varname:
            value =  mandatory_variables["dz"]
        if "O3.MIXING.RATIO.VOLUME_DERIVED" == varname:
            value =  mandatory_variables["o3mr"]
        if "O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD" == varname:
            value =  mandatory_variables["uo3mr"]
        if "O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD" == varname:
            value =  mandatory_variables["uo3mrrand"]
        if "O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD" == varname:
            value =  mandatory_variables["uo3mrsyst"]
        if "PRESSURE_INDEPENDENT" == varname:
            value =  mandatory_variables["xp"]
        if "TEMPERATURE_INDEPENDENT" == varname:
            value =  mandatory_variables["xt"]
        if "PRESSURE_INDEPENDENT_SOURCE" == varname:
            value =  mandatory_variables["xpsce"]
        if "TEMPERATURE_INDEPENDENT_SOURCE" == varname:
            value =  mandatory_variables["xtsce"]

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DF.CUTOFF" == varname:
            try:
                value =  optional_variables["dzdf"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.ORIGINATOR" == varname:
            try:
                value =  optional_variables["dzorig"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.ORIGINATOR" == varname:
            try:
                value =  optional_variables["uo3ndorig"]
            except: # Skip to the next variable if not present
                continue

        if "O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.ORIGINATOR" == varname:
            try:
                value =  optional_variables["uo3mrorig"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DISTANCE.FROM.IMPULSE" == varname:
            try:
                value =  optional_variables["irdist"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE" == varname:
            try:
                value =  optional_variables["ir"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DF.NORMALIZED.FREQUENCY" == varname:
            try:
                value =  optional_variables["tffreq"]
            except: # Skip to the next variable if not present
                continue

        if "O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DF.TRANSFER.FUNCTION" == varname:
            try:
                value =  optional_variables["tf"]
            except: # Skip to the next variable if not present
                continue

        if "SOURCE.PRODUCT" == varname:
            try:
                value =  optional_variables["sceprod"]
            except: # Skip to the next variable if not present
                print " We skip SOURCE PRODUCT"
                continue
        print varname


        # We add the optional data into the varlist values


        # We check the bounds for the values

    
        # We save the data
        typ = varlist[varname]["VAR_DATA_TYPE"]
	print "We write the variable", varname
#	print type(value)
        value = BoundarySave(varlist[varname], value)
#	print value

        maxi = varlist[varname]["VAR_VALID_MAX"]
        mini = varlist[varname]["VAR_VALID_MIN"]
        fill = varlist[varname]["VAR_FILL_VALUE"]
        #maxi = maxi.strip()
        #mini = mini.strip()
        #fill = fill.strip()

        if "STRING" == typ:
            valsize = 1
	    if (varname == "TEMPERATURE_INDEPENDENT_SOURCE") or (varname == "PRESSURE_INDEPENDENT_SOURCE"):
		    valsize = "1" #len(mandatory_variables["z"])
	    try:
               v = d.create(varname, SDC.CHAR, value.shape)
	       setattr(v, "VAR_SIZE", str(valsize))
	    except:
               v = d.create(varname, SDC.CHAR, (1))
	       setattr(v, "VAR_SIZE",str(valsize))
            g1 = v.attr("VAR_VALID_MAX")

            maxi = maxi.strip()
            mini = mini.strip()
            fill = fill.strip()
            if maxi == "":
                maxi = " "
            if mini == "":
                mini = " "
            if fill == "":
                fill = " "

            g1.set(SDC.CHAR, maxi)
            g2 = v.attr("VAR_VALID_MIN")
            g2.set(SDC.CHAR, mini)
            g3 = v.attr("VAR_FILL_VALUE")
            g3.set(SDC.CHAR, fill)


	valsize = ""
	try:
	   valsize = str(value.shape[0]) #";".join([str(i) for i in value.shape])
	   if value.shape[1] >1:
		   valsize += ";" + str(value.shape[1])

	   print valsize
	   print value.shape
	except:
	   try:
	      valsize = str(len(value))
	   except:
	      valsize = "1"




        if "REAL" == typ:
	    try:
               v = d.create(varname, SDC.FLOAT32, value.shape)
	    except:
               v = d.create(varname, SDC.FLOAT32, (1))
	    try:
	    	value = value.astype(float32)
	    except:
		value = float32(value)
		pass
	    setattr(v, "VAR_SIZE", valsize)
            g1 = v.attr("VAR_VALID_MAX")
            g1.set(SDC.FLOAT32, float(maxi))
            g2 = v.attr("VAR_VALID_MIN")
            g2.set(SDC.FLOAT32, float(mini))
            g3 = v.attr("VAR_FILL_VALUE")
            g3.set(SDC.FLOAT32, float(fill))


        if "DOUBLE" == typ:
	    try:
               v = d.create(varname, SDC.FLOAT64, value.shape)
	    except:
               v = d.create(varname, SDC.FLOAT64, (1))
	    setattr(v, "VAR_SIZE", valsize)
#            v = d.create(varname, SDC.FLOAT64, value.shape)
 #           v = d.create(varname, SDC.CHAR, value.shape)
            g1 = v.attr("VAR_VALID_MAX")
            g1.set(SDC.FLOAT64, float(maxi))
            g2 = v.attr("VAR_VALID_MIN")
            g2.set(SDC.FLOAT64, float(mini))
            g3 = v.attr("VAR_FILL_VALUE")
            g3.set(SDC.FLOAT64, float(fill))
        if "SHORT" == typ:
	    try:
               v = d.create(varname, SDC.INT16, value.shape)
	    except:
               v = d.create(varname, SDC.INT16, (1))
	    setattr(v, "VAR_SIZE", valsize)
#            v = d.create(varname, SDC.INT16, value.shape)
 #           v = d.create(varname, SDC.CHAR, value.shape)
            g1 = v.attr("VAR_VALID_MAX")
            g1.set(SDC.INT16, int(maxi))
            g2 = v.attr("VAR_VALID_MIN")
            g2.set(SDC.INT16, int(mini))
            g3 = v.attr("VAR_FILL_VALUE")
            g3.set(SDC.INT16, int(fill))
        if "INTEGER" == typ:
	    try:
               v = d.create(varname, SDC.INT32, value.shape)
	    except:
               v = d.create(varname, SDC.INT32, (1))
	    setattr(v, "VAR_SIZE", valsize)

            #v = d.create(varname, SDC.INT32, value.shape)
            g1 = v.attr("VAR_VALID_MAX")
            g1.set(SDC.INT32, int(maxi))
            g2 = v.attr("VAR_VALID_MIN")
            g2.set(SDC.INT32, int(mini))
            g3 = v.attr("VAR_FILL_VALUE")
            g3.set(SDC.INT32, int(fill))        
        print "We set the attributes"
            
	setattr(v, "VAR_NAME", varname)

 #       try:
#            if len(v.shape) < 2:
	#        setattr(v, "VAR_SIZE", int(v.size()))
  #          else:
#	        setattr(v, "VAR_SIZE", )))

	for k in varlist[varname]:
            if not(k == "VAR_VALID_MAX" or k == "VAR_VALID_MIN" or  k == "VAR_FILL_VALUE" or k == "VAR_SIZE"):
                tmp =  varlist[varname][k]
                tmp = tmp.strip()
                if tmp == "":
                    tmp = " "

                setattr(v, k, tmp)
#	try:
#		print "We set the value", value.min(), value.max(), value.dtype, value.shape
#		
#	except:
#		print "We set the value", value
	v[:] = value
#		v[k] = varlist[varname][k]
#        print "We continue"
    # We populate the output file from the meta-data file
    
    for i in attributes.keys():
	print "ATTRIBUTES", i
        tmp = d.attr(i)
        tmp2 = attributes[i]
        tmp2.strip()
        if tmp2 == "":
            tmp2 = " "
        tmp.set(SDC.CHAR, tmp2) #attributes[i])
    
    return



