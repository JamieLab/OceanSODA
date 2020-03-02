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


outputPath = "../algorithms/algo_data/lee2000_masks_tmh.nc";

petersMaskPath = "Peters_LeeEtAl2000_masks.nc";
pmasks = petersMaskPath.variables["regionNo"][:];
P_BAY_OF_BENGAL = 3;
P_ZONE1_MINUS_NORTH_PACIFIC = 1;
P_ZONE1_NOTH_PACIFIC_ONLY = 2;

oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);

#ocean mask values
LAND=0;
ATLANTIC=30;
PACIFIC=70;
INDIAN=50;
SOUTHERN=90;

def get_zone1_mask():
    pacificRegion = np.zeros((180, 360), dtype=int);
    pacificRegion[90-20:90+20, 180-110:180-75] = 1;
    pacificRegion[90-10:90+10, 180-160:180-110] = 1;
    pacificRegion[oceans!=PACIFIC] = 0;
    atlanticIndianRegion = np.zeros((180, 360), dtype=int);
    atlanticIndianRegion[90-5:90+5, :] = 1;
    atlanticIndianRegion[(oceans==PACIFIC)] = 0;
    zone1 = (pacificRegion | atlanticIndianRegion).astype(int);
    return zone1;


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
zone1eq1 = get_zone1_mask();
northPacific = (oceans==PACIFIC);
northPacific[:90,:] = False;
zone1eq1[northPacific] = 0; #Eq 1 is not applicable to the north pacific.
var = nc.createVariable("zone1_equation1_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq1;

#Zone 1, equation 2
zone1eq2 = (get_zone1_mask() & northPacific).astype(int);
var = nc.createVariable("zone1_equation2_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1eq2;

#Zone 1, equation 3



nc.close();


a = zone1eq1;
a += zone1eq2*2
a[oceans==LAND] = -1;
plt.figure();
plt.imshow(np.flipud(a))

##sanity check for overlaps
#import matplotlib.pyplot as plt;
#zoneSum = zone1+zone2+zone3+zone4+zone5;
#plt.figure();
#plt.imshow(zone2-oceanMask);
#plt.colorbar();
#plt.ylim(0,180)



peter = Dataset("/home/verwirrt/Projects/Work/20190816_OceanSODA/scripts/algorithms/algo_data/unused_from_peter/LeeEtAl2000.nc", 'r').variables["regionNo"][:];
plt.figure();
plt.imshow(peter)














