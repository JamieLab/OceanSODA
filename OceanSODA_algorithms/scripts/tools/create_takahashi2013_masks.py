#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:56:01 2020

Selected zone definitions from Takahashi and Sutherland 2013
See: Takahashi, T., Sutherland, S.C., Chipman, D.W., Goddard, J.G., Ho, C., Newberger, T., Sweeney, C. and Munro, D.R., 2014. Climatological distributions of pH, pCO2, total CO2, alkalinity, and CaCO3 saturation in the global surface ocean, and temporal changes at selected locations. Marine Chemistry, 164, pp.95-125.

Extents taken from table 1

@author: tom holding
"""

import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;

outputPath = "../algorithms/algo_data/takahashi2013_masks_tmh.nc";

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


#Zone 7: North Atlantic Drift
zone7 = np.zeros((180, 360), dtype=float);
zone7[90+40:90+55, 180-60:180+10] = 1; #40N-55N, 60W-10E
zone7[oceans!=ATLANTIC] = 0; #Remove non-atlantic regions
zone7[90+46:90+54, 180-70:180-50] = 1; #Gulf of St Lawrence should be included (also means the algorithm can be a single algorithm - barring Mediterranean)
var = nc.createVariable("zone7_north_atlantic_drift", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone7;

#Zone 8: Central Atlantic
zone8 = np.zeros((180, 360), dtype=float);
zone8[90-40:90+40, :] = 1; #40S-40N, coast-to-coast
zone8[oceans!=ATLANTIC] = 0; #Remove non-atlantic regions
var = nc.createVariable("zone8_central_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone8;

nc.close();



#sanity check for overlaps
zoneSum = zone7+zone8
plt.figure();
plt.imshow(zoneSum);
plt.colorbar();
plt.ylim(0,180)

