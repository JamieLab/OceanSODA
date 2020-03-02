#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:56:01 2020

Zone definitions from table 1 of Lee et al 2006.
See: Lee, K., Tong, L.T., Millero, F.J., Sabine, C.L., Dickson, A.G., Goyet, C., Park, G.H., Wanninkhof, R., Feely, R.A. and Key, R.M., 2006. Global relationships of total alkalinity with salinity and temperature in surface waters of the world's oceans. Geophysical research letters, 33(19).

@author: tom holding
"""

import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;

#tmp
#ncin = Dataset("Peters_LeeEtAl2006_masks.nc", 'r');
#pregions = ncin.variables["regionNo"][:]; #Peter's region data
#regionNames = ncin.variables["regionName"][:];


outputPath = "../algorithms/algo_data/lee2006_masks_tmh.nc";

oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);

#ocean mask values
LAND=0;
ATLANTIC=30;
PACIFIC=70;
INDIAN=50;
SOUTHERN=90;

#plt.figure();
#plt.imshow(oceans);


#create NC file and dimensions
nc = Dataset(outputPath, 'w');
nc.createDimension("lat", 180);
nc.createDimension("lon", 360);

#dimension variables
var = nc.createVariable("lat", float, ("lat",));
var.units = "lat (degrees North)";
var[:] = np.arange(-90, 90)+0.5;
var = nc.createVariable("lon", float, ("lon",));
var.units = "lon (degrees East)";
var[:] = np.arange(-180, 180)+0.5;


#Zone 1: Subtropics
zone1 = np.zeros((180, 360), dtype=float);
zone1[90-30:90+30, :] = 1; #30S - 30N
zone1[90-20:90+20, 180-110:180-75] = 0; #Exclude Equatorial Pacific upwelling region
zone1[90-10:90+10, 180-140:180-110] = 0; #Exclude Equatorial Pacific upwelling region
var = nc.createVariable("zone1_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1;

#Equatorial Pacific upwelling region
zone2 = np.zeros((180, 360), dtype=float);
zone2[90-20:90+20, 180-110:180-75] = 1; #Equatorial Pacific upwelling region
zone2[90-10:90+10, 180-140:180-110] = 1; #Equatorial Pacific upwelling region
zone2[oceans!=PACIFIC] = 0; #Remove non-pacific regions
var = nc.createVariable("zone2_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone2;

#North Atlantic
zone3 = np.zeros((180, 360), dtype=float);
zone3[90+30:90+80, 180-100:180+106] = 1; #Using Peter's mask's boundaries
var = nc.createVariable("zone3_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone3;

#North Pacific
zone4 = np.zeros((180, 360), dtype=float);
zone4[90+30:90+90+1, 180-180:180-100] = 1; #Using Peter's mask's boundaries
zone4[90+30:90+90+1, 180+106:180+180] = 1; #Using Peter's mask's boundaries
var = nc.createVariable("zone4_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone4;

#Southern Ocean
zone5 = np.zeros((180, 360), dtype=float);
zone5[90-80:90-30, :] = 1;
var = nc.createVariable("zone5_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone5;

nc.close();



#sanity check for overlaps
import matplotlib.pyplot as plt;
zoneSum = zone1+zone2+zone3+zone4+zone5;
plt.figure();
plt.imshow(zoneSum);
plt.colorbar();
plt.ylim(0,180)

