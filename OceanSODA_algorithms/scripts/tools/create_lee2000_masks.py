#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:56:01 2020

Zone definitions from table 3, 4 and 5,  of Lee et al 2000.
See: Lee, K., Wanninkhof, R., Feely, R.A., Millero, F.J. and Peng, T.H., 2000. Global relationships of total inorganic carbon with temperature and nitrate in surface seawater. Global Biogeochemical Cycles, 14(3), pp.979-994.

@author: tom holding
"""

import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;

#tmp
#ncin = Dataset("Peters_LeeEtAl2006_masks.nc", 'r');
#pregions = ncin.variables["regionNo"][:]; #Peter's region data
#regionNames = ncin.variables["regionName"][:];


outputPath = "../algorithms/algo_data/lee2000_masks_tmh.nc";


oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);

#ocean mask values
LAND=0;
ATLANTIC=30;
PACIFIC=70;
INDIAN=50;
SOUTHERN=90;

petersMaskPath = "Peters_LeeEtAl2000_masks.nc";
pmasks = np.flipud(Dataset(petersMaskPath, 'r').variables["regionNo"][:]);
P_ZONE1_MINUS_NORTH_PACIFIC = 1;
P_ZONE1_NOTH_PACIFIC_ONLY = 2;
P_ZONE1_BAY_OF_BENGAL = 3;
P_ZONE2_SUBTROPICS_MINUS_WEST_ATLANTIC_AND_NORTH_PACIFIC = 4;
P_ZONE2_SUBTROPICS_NORTH_PACIFIC = 5;
P_ZONE2_SUBTROPICAL_WESTERN_ATLANTIC = 6;
P_ZONE3_NORTH_WESTERN_ATLANTIC = 7;
P_ZONE3_NORTH_EASTERN_ATLANTIC = 8;
P_ZONE4_NORTH_PACIFIC = 9;
P_ZONE5_SOUTHER_OCEAN = 10;

pmasks[oceans==LAND]= 0;
#plt.figure(); plt.imshow(np.flipud(pmasks));



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


#Zone 1, equation 1
zone1eq1 = (pmasks==P_ZONE1_MINUS_NORTH_PACIFIC).astype(int);
var = nc.createVariable("zone1_equation1_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq1;

#Zone 1, equation 2
zone1eq2 = (pmasks==P_ZONE1_NOTH_PACIFIC_ONLY).astype(int);
var = nc.createVariable("zone1_equation2_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq2;

#Zone 1, equation 3
zone1eq3 = (pmasks==P_ZONE1_MINUS_NORTH_PACIFIC).astype(int);
var = nc.createVariable("zone1_equation3_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq3;

#Zone 1, equation 4
zone1eq4 = (pmasks==P_ZONE1_NOTH_PACIFIC_ONLY).astype(int);
var = nc.createVariable("zone1_equation4_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq4;

#Zone 1, equation 5
zone1eq5 = (pmasks==P_ZONE1_BAY_OF_BENGAL).astype(int);
var = nc.createVariable("zone1_equation5_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq5;


#Zone 2, equation 6S and 6W
zone2eq6 = (pmasks==P_ZONE2_SUBTROPICS_MINUS_WEST_ATLANTIC_AND_NORTH_PACIFIC).astype(int);
var = nc.createVariable("zone2_equation6_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone2eq6

#Zone 2, equation 7S and 7W
zone2eq7 = (pmasks==P_ZONE2_SUBTROPICS_NORTH_PACIFIC).astype(int);
var = nc.createVariable("zone2_equation7_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone2eq7

#Zone 2, equation 8
zone2eq8 = (pmasks==P_ZONE2_SUBTROPICAL_WESTERN_ATLANTIC).astype(int);
var = nc.createVariable("zone2_equation8_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone2eq8

#Zone3, equation 9S and 9W
zone3eq9 = (pmasks==P_ZONE3_NORTH_WESTERN_ATLANTIC).astype(int);
var = nc.createVariable("zone3_equation9_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone3eq9

#Zone3, equation 10S and 10W
zone3eq10 = (pmasks==P_ZONE3_NORTH_EASTERN_ATLANTIC).astype(int);
var = nc.createVariable("zone3_equation10_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone3eq10


#Zone4, equation 11S and 11W
zone4eq11 = (pmasks==P_ZONE4_NORTH_PACIFIC).astype(int);
var = nc.createVariable("zone4_equation11_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone4eq11

#Zone5, equation 12S and 12W
zone5eq12 = (pmasks==P_ZONE5_SOUTHER_OCEAN).astype(int);
var = nc.createVariable("zone5_equation12_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone5eq12


nc.close();


#sanity check for overlaps
zoneSum = zone1eq1+zone1eq2+zone1eq3+zone1eq4+zone1eq5+zone2eq6+zone2eq7+zone2eq8+zone3eq9+zone3eq10+zone4eq11+zone5eq12;
plt.figure();
plt.imshow(zoneSum);
plt.colorbar();
plt.ylim(0,180)
plt.title("note: equatorial region should equal 2");














