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


outputPath = "../algorithms/algo_data/millero1998_masks_tmh.nc";
oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);
#plt.figure(); plt.imshow(oceans); plt.colorbar();

#ocean mask values
LAND=0;
ATLANTIC=30;
PACIFIC=70;
INDIAN=50;
SOUTHERN=90;


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


#Zone 1 atlantic and indian: Equitorial and mid latitudes
atlantic = np.zeros((180, 360), dtype=float);
atlantic[90-30:90+30, :] = 1; #30S - 30N
atlantic[oceans!=ATLANTIC] = 0; #Exclude non-atlantic regions
indian = np.zeros((180, 360), dtype=float);
indian[90-30:90+25, :] = 1; #30S - 30N
indian[oceans!=INDIAN] = 0; #Exclude non-atlantic regions
zone1 = atlantic+indian;
var = nc.createVariable("zone1_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone1;

#Zone 2: North Atlantic
zone2 = np.zeros((180, 360), dtype=float);
zone2[90+30:90+80, :] = 1; #20N - 80N
zone2[oceans!=ATLANTIC] = 0; #Exclude non-atlantic regions
#add an extra bit that covers the Gulf of St Lawrens (as this isn't included as 'Atlantic' by the ATLANTIC mask)
zone2[90+46:90+54, 180-70:180-50] = 1; #Gulf of St Lawrence should be included (also means Millero1998 can be a single algorithm)
#plt.figure(); plt.imshow(np.flipud(zone2));
var = nc.createVariable("zone2_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone2;

#Zone 3 pacific: Equitorial pacific upwelling region
zone3 = np.zeros((180, 360), dtype=float);
zone3[90-20:90+20, 180-110:180-75] = 1; #110W-75W, 20S-20N
zone3[90-10:90+10, 180-140:180-110] = 1; #140W-110W, 10N-10S
zone3[oceans!=PACIFIC] = 0; #Exclude non-atlantic regions
var = nc.createVariable("zone3_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone3;

#Zone 4: pacific Gyres (excluding equatorial upwelling)
zone4 = np.zeros((180, 360), dtype=float);
zone4[90-20:90+30, :] = 1;
zone4[90-20:90+20, 180-110:180-75] = 0; #110W-75W, 20S-20N
zone4[90-10:90+10, 180-140:180-110] = 0; #140W-110W, 10N-10S
zone4[oceans!=PACIFIC] = 0; #Exclude non-atlantic regions
var = nc.createVariable("zone4_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone4;

#Zone 5 pacific: North
zone5 = np.zeros((180, 360), dtype=float);
zone5[90+30:90+52, :] = 1;
zone5[oceans!=PACIFIC] = 0; #Exclude non-atlantic regions
var = nc.createVariable("zone5_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone5;

#Zone 6 Southern Ocean: Atlantic and Indian
zone6a = np.zeros((180, 360), dtype=float);
zone6a[90-70:90-30, :] = 1;
zone6a[(oceans==ATLANTIC) | (oceans==INDIAN)==False] = 0; #Exclude non-atlantic and non-indian regions
tmp = oceans==SOUTHERN; #Southern ocean mask
tmp[:, 0:113] = 0; #Just the section below atlantic and indian oceans
tmp[:, 328:] = 0; #Just the section below atlantic and indian oceans
tmp[:90-70,:] = 0; #Don't include below 70S
zone6a += tmp;
var = nc.createVariable("zone6_atlantic_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone6a;

#Zone 6 Southern Ocean: Pacific
zone6p = np.zeros((180, 360), dtype=float);
zone6p[90-70:90-30, :] = 1;
zone6p[oceans!=PACIFIC] = 0;
tmp = oceans==SOUTHERN;
tmp[:, 113:328] = 0; #exclude southern ocean between atlantic and indian
tmp[:90-70,:] = 0; #Don't include below 70S
zone6p += tmp;
#plt.figure(); plt.imshow(zone6p)
var = nc.createVariable("zone6_pacific_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zone6p;


#sanity check for overlaps
zoneSum = zone1+zone2+zone3+zone4+zone5+zone6a+zone6p;
var = nc.createVariable("mask_sum_all_zones", int, ("lat", "lon"));
var.units = "Integer";
var[:] = zoneSum;

nc.close();


plt.figure();
plt.imshow(zoneSum);
plt.colorbar();
plt.ylim(0,180)



