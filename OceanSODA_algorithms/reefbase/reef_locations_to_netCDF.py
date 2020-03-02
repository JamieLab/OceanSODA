#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:41:50 2020

Processes the reef base reef location data into netCDF format

@author: tom holding
"""

import pandas as pd;
import numpy as np;
from os import path;
from netCDF4 import Dataset, stringtochar;

inputPath = "ReefLocations.csv";
outputPath = "reef_locations.nc";

df = pd.read_csv(inputPath, sep=",", encoding="ISO-8859-1");

#find max string sizes
maxSizeSystem = 0;
for s in df["REEF_SYSTEM"]:
    if isinstance(s, str):
        if len(s) > maxSizeSystem:
            maxSizeSystem = len(s);
maxSizeName = 0;
for s in df["REEF_NAME"]:
    if isinstance(s, str):
        if len(s) > maxSizeName:
            maxSizeName = len(s);
stringLength = max(maxSizeSystem, maxSizeName);
stringFormat = "S"+str(stringLength);

#systemStrs = [];
#for s in df["REEF_SYSTEM"]:
#    arr = stringtochar(np.array(s, stringFormat));
#    
#systemStrs = [stringtochar(np.array(s, stringFormat)) for s in df["REEF_SYSTEM"] if isinstance(s, str)];
#
#s = np.array(["test"], "S100");
#stringtochar(s)



#Write output to netCDF
nc = Dataset(outputPath, 'w');
nc.createDimension("index", len(df));


#dimension variables
var = nc.createVariable("index", int, ("index",));
var.units = "integer";
var[:] = range(0, len(df));


#data variables
var = nc.createVariable("lat", float, ("index",));
var.units = "lat (degrees North)";
var[:] = np.array(df["LAT"]);

var = nc.createVariable("lon", float, ("index",));
var.units = "lon (degrees East)";
var[:] = np.array(df["LON"]);

var = nc.createVariable("id", float, ("index",));
var.units = "integer";
var[:] = np.array(df["ID"]);

#var = nc.createVariable("reef_system", str, ("index",));
#var.units = "string";
#var[:] = np.array([stringtochar[s] for s in df["REEF_SYSTEM"]]);
#
#var = nc.createVariable("reef_name", str, ("index",));
#var.units = "string";
#var[:] = np.array([stringtochar[s] for s in df["REEF_NAME"]]);


nc.close();
